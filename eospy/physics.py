
import numpy as np
import shen
import fermion
import units
import quantities as pq


class EquationOfStateTerms(object):
    """
    Base class for EOS terms, which in general are responsible for returning the
    number and mass density, pressure, internal energy, entropy, and chemical
    potential of some component of the gas.

    Each instantiation of classes inheriting from this represents one point in
    the space of independent thermodynamic variables.
    """
    kB = pq.constants.Boltzmann_constant.rescale('MeV/K')
    hc = (pq.constants.hbar * pq.c).rescale('MeV*fm')
    me = (pq.constants.electron_mass * pq.c**2).rescale('MeV')

    def number_density(self):
        return self._gencall('n')

    def mass_density(self):
        raise NotImplementedError("mass density not defined for this term")

    def enthalpy(self):
        """
        Generates the enthalpy parameter rho*h, where h = 1 + e + p/rho.
        """
        D = self.mass_density()
        u = self.internal_energy()
        p = self.pressure()
        return D + (u + p)/pq.c**2

    def pressure(self):
        return self._gencall('p')

    def internal_energy(self):
        return self._gencall('u')

    def entropy(self):
        return self._gencall('s')

    def specific_internal_energy(self):
        return self._terms['u'] / self.mass_density()

    def internal_energy_per_particle(self):
        return self._terms['u'] / self._terms['n']

    def _gencall(self, key):
        return self._terms[key]

    def temperature_in_MeV(self, val):
        """
        Convenience function to convert temperature to energy.
        """
        if isinstance(val.dimensionality.keys()[0], pq.UnitTemperature):
            return (val * self.kB).rescale('MeV')
        else:
            return val.rescale('MeV')

    def temperature_in_Kelvin(self, val):
        """
        Convenience function to convert energy to temperature.
        """
        if isinstance(val.dimensionality.keys()[0], pq.UnitTemperature):
            return val.rescale('K')
        else:
            return (val / self.kB).rescale('K')



class EquationOfStateEvaluator(object):
    """
    Instances of this class are composite equations of state, consisting of
    several different terms, or components. Each term is a class inheriting from
    EquationOfStateTerms, not an instance of a class. It knows how to create
    instances when it needs and query them for values or create derivatives.
    """
    _units = {
        'number_density': 1/pq.m**3,
        'pressure': pq.pascal,
        'internal_energy': pq.J/pq.m**3,
        'specific_internal_energy': pq.J/pq.kg,
        'internal_energy_per_particle': pq.J,
        'entropy': pq.J/pq.K,
        'enthalpy': pq.kg/pq.m**3
        }
    _num_deriv_dx = 1e-8

    def build_terms(self, args):
        """
        This method should be over-ridden to return a list of
        EquationOfStateTerms instances which are built from *args*.
        """
        raise NotImplementedError()

    def set_numerical_derivative_step(self, dx):
        """
        Sets the size of the step used for the numerical derivative
        estimates. The value is a fraction of the independent variable being
        stepped.
        """
        self._num_deriv_step = dx

    def _call_attr(self, args, attr):
        return sum([getattr(term, attr)() for term in
                    self.build_terms(args)], 0.0 * self._units[attr])

    def number_density(self, *args):
        return self._call_attr(args, 'number_density')

    def pressure(self, *args):
        return self._call_attr(args, 'pressure')

    def internal_energy(self, *args):
        return self._call_attr(args, 'internal_energy')

    def specific_internal_energy(self, *args):
        return self._call_attr(args, 'specific_internal_energy')

    def internal_energy_per_particle(self, *args):
        return self._call_attr(args, 'internal_energy_per_particle')

    def enthalpy(self, *args):
        return self._call_attr(args, 'enthalpy')

    def entropy(self, *args):
        return self._call_attr(args, 'entropy')

    def derivative(self, attr, var, *args):
        dx = self._num_deriv_dx
        X0, X1 = list(args), list(args)

        n = self._vars[var]

        X1[n] = X1[n]*(1.0 + dx)
        X0[n] = X0[n]*(1.0 - dx)

        f = lambda x: self._call_attr(x, attr)
        return (f(X1) - f(X0)) / (X1[n] - X0[n])

    def gamma_effective(self, *args, **kwargs):
        """
        Compute the effective Gamma := (dlogp/dlogD)|_s to be used in the sound
        speed: c_s^2 = (\Gamma*p) / (rho*h) using one of three different
        methods.
        """
        method = kwargs.get('method', 2)

        # Some terms may use number density 'n', while others may use the mass
        # density 'D'.
        den = self._density_var
        ekey = ('specific_internal_energy' if den is 'D' else
                'internal_energy_per_particle')

        if method == 1:
            n = args[self._vars[den]]
            T = args[self._vars['T']]
            p = self.pressure(*args)
            s = self.entropy(*args)
            dpdn = self.derivative('pressure', den, *args)
            dpdT = self.derivative('pressure', 'T', *args)
            dsdn = self.derivative('entropy', den, *args)
            dsdT = self.derivative('entropy', 'T', *args)
            return (n/p) * (dpdn*dsdT - dpdT*dsdn) / dsdT
        elif method == 2:
            n = args[self._vars[den]]
            T = args[self._vars['T']]
            p = self.pressure(*args)
            dpdn = self.derivative('pressure', den, *args)
            dpdT = self.derivative('pressure', 'T', *args)
            dedn = self.derivative(ekey, den, *args)
            dedT = self.derivative(ekey, 'T', *args)
            return (n/p) * (dpdn*dedT - dpdT*dedn + p/n**2 * dpdT) / dedT
        elif method == 3:
            n = args[self._vars[den]]
            T = args[self._vars['T']]
            p = self.pressure(*args)
            dpdn = self.derivative('pressure', den, *args)
            dpdT = self.derivative('pressure', 'T', *args)
            dedn = self.derivative(ekey, den, *args)
            dedT = self.derivative(ekey, 'T', *args)
            return (n/p) * (dpdn*dedT + T/n**2 * dpdT**2) / dedT
        else:
            raise ValueError("method must be 1,2, or 3")

    def sound_speed(self, *args, **kwargs):
        p = self.pressure(*args)
        rhoh = self.enthalpy(*args)
        gam = self.gamma_effective(*args, **kwargs)
        return np.sqrt(gam * p / rhoh).rescale('m/s')


class IdealAdiabatic(EquationOfStateTerms):
    """
    Represents an adiabatic equation of state, with adiabatic constant 'gamma'.
    """
    def __init__(self, n, T, gamma=1.4, particle_mass=pq.constants.proton_mass):
        self.n = n
        self.T = self.temperature_in_Kelvin(T)
        self.gamma = gamma
        self.particle_mass = particle_mass
        self._terms = { }
        self._set_terms()

    def mass_density(self):
        return self._terms['n'] * self.particle_mass

    def _set_terms(self):
        g1 = self.gamma - 1.0

        n = self.n
        p = self.n * self.kB * self.T
        s = np.log(p.magnitude/n.magnitude**self.gamma) / g1

        f = self._terms
        f['n'] = n
        f['p'] = p
        f['u'] = p / g1
        f['s'] = s * self.kB



class BlackbodyPhotons(EquationOfStateTerms):
    """
    Evaluates all terms for a photon gas. Only needs the temperature as a
    parameter.
    """
    def __init__(self, T):
        self.kT = self.temperature_in_MeV(T)
        self._terms = { }
        self._set_terms()

    def mass_density(self):
        return 0.0

    def chemical_potential(self):
        return 0.0

    def _set_terms(self):
        """
        http://en.wikipedia.org/wiki/Photon_gas
        """
        z3 = 1.202056903159594285 # RiemannZeta(3), Mathematica
        T3 = np.power(self.kT, 3)
        T4 = np.power(self.kT, 4)
        a = pow(np.pi, 2) / (15*np.power(self.hc, 3))

        f = self._terms
        f['n'] = T3 * 2 * z3 / (np.power(np.pi, 2.0) * np.power(self.hc, 3.0))
        f['p'] = T4 * a / 3.0
        f['u'] = T4 * a
        f['s'] = (4./3.) * (T4 * a) / (self.kT / self.kB)



class NeutrinoComponent(EquationOfStateTerms):
    """
    Represents an EOS component of neutrinos and anti-neutrinos. Must be built
    with the desired chemical potential, which would typically be recovered by
    through coupling to a baryon component.
    """
    def __init__(self, mu, T):
        """
        Parameters:
        --------------------------------------------------------

        mu   : chemical potential (MeV)
        T    : temperature (MeV)
        """
        self.mu = mu
        self.kT = self.temperature_in_MeV(T)
        self._terms = { }
        self._set_terms()

    def mass_density(self):
        """
        Neutrinos are massless.
        """
        return 0.0

    def chemical_potential(self):
        return self._terms['mu']

    def _set_terms(self):
        Volume = np.power(np.pi, 2) * np.power(self.hc/self.kT, 3)
        Energy = self.kT

        f = fermion.neutrino_everything(self.mu/self.kT)

        f['n'] *= (1.0 / Volume)
        f['p'] *= (Energy / Volume)
        f['u'] *= (Energy / Volume)
        f['s']  = (f['u'] + f['p']) / self.kT - f['n'] * f['eta']

        for k in "npus":
            self._terms[k] = f[k]

        self._terms['eta'] = f['eta']
        self._terms['mu'] = self._terms['eta'] * self.kT


class FermionComponent(EquationOfStateTerms):
    """
    Represents an EOS component of electrons or positrons.
    """

    _volume = (np.power(np.pi, 2) *
               np.power(EquationOfStateTerms.hc/EquationOfStateTerms.me, 3))
    _energy = EquationOfStateTerms.me

    def __init__(self, sgn, mu, T):
        """
        Parameters:
        --------------------------------------------------------

        mu   : chemical potential, without rest-mass (MeV)
        T    : temperature (MeV)
        """
        self.sgn = sgn
        self.mu = mu
        self.kT = self.temperature_in_MeV(T)
        self._terms = { }
        self._set_terms()

    def mass_density(self):
        """
        Returns the mass density of electrons.
        """
        return self._terms['n'] * pq.constants.electron_mass

    def chemical_potential(self):
        return self._terms['mu']

    def _set_terms(self):
        """
        Sets the number density, pressure, internal energy, and specific entropy
        (n,p,u,s) for both electrons and positrons.
        """
        eta = self.mu / self.kT
        beta = self.kT / self.me
        f = fermion.fermion_everything(self.sgn, eta, beta)

        f['n'] *= (1.0 / self._volume)
        f['p'] *= (self._energy / self._volume)
        f['u'] *= (self._energy / self._volume)
        f['s']  = (f['u'] + f['p']) / self.kT - f['n'] * eta

        for k in "npus":
            self._terms[k] = f[k]

        self._terms['eta'] = f['eta'] # without rest-mass
        self._terms['mu'] = self._terms['eta'] * self.kT + self.me # with rest-mass


class FermiDiracElectrons(FermionComponent):
    """
    Evaluates the electron thermodynamics using exact Fermi-Dirac integrals.
    """
    def __init__(self, np, T):
        """
        Parameters:
        --------------------------------------------------------
        T    : temperature (MeV)
        np   : the number density of positively charged baryons (1/fm^3)
        """
        kT = self.temperature_in_MeV(T)
        eta = fermion.solve_eta_pairs(kT / self.me, self._volume * np)
        super(FermiDiracElectrons, self).__init__(+1, eta*kT, kT)



class FermiDiracPositrons(FermionComponent):
    """
    Evaluates the positron thermodynamics using exact Fermi-Dirac integrals.
    """
    def __init__(self, np, T):
        """
        Parameters:
        --------------------------------------------------------
        T    : temperature (MeV)
        np   : the number density of positively charged baryons (1/fm^3)
        """
        kT = self.temperature_in_MeV(T)
        eta = fermion.solve_eta_pairs(kT / self.me, self._volume * np)
        super(FermiDiracPositrons, self).__init__(-1, eta*kT, kT)



class ColdElectrons(FermionComponent):
    """
    Evaluates the thermodynamic variables for a cold (fully degenerate) electron
    gas. Accurate as long as the Fermi momentum is not relativistic.
    """
    def __init__(self, ne):
        """
        Parameters:
        --------------------------------------------------------
        ne   : the number density of electrons (1/fm^3)
        """
        self.ne = ne
        self._terms = { }
        self._set_terms()

    def _set_terms(self):
        """
        http://scienceworld.wolfram.com/physics/ElectronDegeneracyPressure.html
        """
        num = np.power(np.pi, 2) * np.power(self.hc, 2)
        den = 5.0 * self.me
        las = np.power(3.0/np.pi, 2./3.) * np.power(self.ne, 5./3.)

        self._terms['n'] = self.ne
        self._terms['p'] = (num/den) * las
        # 'u' not implemented yet
        # 's' not implemented yet



class DenseElectrons(FermionComponent):
    """
    Evaluates the thermodynamic variables for a dense and fully degenerate
    electron gas, where the Fermi momentum is ultra-relativistic.
    """
    def __init__(self, ne):
        """
        Parameters:
        --------------------------------------------------------
        ne   : the number density of electrons (1/fm^3)
        """
        self.ne = ne
        self._terms = { }
        self._set_terms()

    def _set_terms(self):
        """
        The Physical Universe: An Introduction to Astronomy, Frank H. Shu (1982)
        """
        self._terms['n'] = self.ne
        self._terms['p'] = 0.123 * 2*np.pi * self.hc * np.power(self.ne, 4./3.)
        # 'u' not implemented yet
        # 's' not implemented yet



class NucleonsShenEos3(EquationOfStateTerms):
    """
    Evaluates the thermodynamic variables for dense baryons using the (updated
    2011) lookup table 'eos3' of Shen et. al. (1998).

    Notes:
    ----------------------------------------------------------------------------

    user guide: http://user.numazu-ct.ac.jp/~sumi/eos/table2/guide_EOS3.pdf
    full table: http://user.numazu-ct.ac.jp/~sumi/eos/table2/eos3.tab.gz
    """

    _table = None # caches the hdf5 table once it is loaded

    def __init__(self, D, T, Y):
        self.D = D.rescale('g/cm^3')
        self.kT = self.temperature_in_MeV(T)
        self.Y = Y
        self._terms = { }
        self._set_terms()

    def chemical_potential(self, type):
        if type == 'neutrons': return self._terms['mu_n']
        if type == 'protons' : return self._terms['mu_p']
        raise ValueError("'type' must be either 'neutrons' or 'protons'")

    def mass_density(self):
        return self._terms['rhoB']

    def _set_terms(self):
        if type(self)._table is None:
            cols = ['log10_rhoB', 'logT', 'Yp', 'p', 'nB', 'Eint',
                    'S', 'un', 'up']
            type(self)._table = shen.read_hdf5("data/eos3.h5", cols=cols)

        D, kT, Ye = self.D.magnitude, self.kT.magnitude, self.Y
        t = self._terms
        e = type(self)._table

        # user guide explains that chemical potentials are with respect to the
        # free-nucleon mass
        Eref = 938.0 * pq.MeV

        t['n'] = shen.sample(e, 'nB'   , D, kT, Ye) / pq.fm**3
        t['p'] = shen.sample(e, 'p'    , D, kT, Ye) * pq.MeV / pq.fm**3
        t['u'] = shen.sample(e, 'Eint' , D, kT, Ye) * t['n'] * pq.MeV
        t['s'] = shen.sample(e, 'S'    , D, kT, Ye) * self.kB # entropy per baryon
        t['mu_n'] = shen.sample(e, 'un', D, kT, Ye) * pq.MeV + Eref
        t['mu_p'] = shen.sample(e, 'up', D, kT, Ye) * pq.MeV + Eref
        t['rhoB'] = 10**shen.sample(e, 'log10_rhoB', D, kT, Ye) * pq.g/pq.cm**3

