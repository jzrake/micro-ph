
import numpy as np
import shen
import fermion
import units


LIGHT_SPEED        = 2.997924580e+10 # cm/s
HBAR_C             = 1.973269718e+02 # MeV-fm
ELECTRON_MASS      = 5.110998928e-01 # MeV
ATOMIC_MASS_UNIT   = 9.314940612e+02 # MeV
PROTON_MASS        = 9.382720462e+02 # MeV
MEV_TO_ERG         = 1.602176487e-06
FM3_TO_CM3         = 1.000000000e-39
BOLTZMANN_CONSTANT = 8.617332400e-11 # MeV/K



class EquationOfStateTerms(object):
    """
    Base class for EOS terms, which in general are responsible for returning the
    number and mass density, pressure, internal energy, entropy, and chemical
    potential of some component of the gas.

    Each instantiation of classes inheriting from this represents one point in
    the space of independent thermodynamic variables.
    """
    def number_density(self):
        return self._terms['n'](unit) if unit else self._terms['n']

    def pressure(self, unit=None):
        return self._terms['p'](unit) if unit else self._terms['p']

    def internal_energy(self, unit=None):
        return self._terms['u'](unit) if unit else self._terms['u']

    def specific_internal_energy(self, unit):
        return self._terms['u'](unit) / self._terms['n'](unit)

    def entropy(self, unit=None):
        return self._terms['s'](unit) if unit else self._terms['s']

    def chemical_potential(self):
        return self._terms['eta'] * self.T


class EquationOfStateEvaluator(object):
    """
    Instances of this class are composite equations of state, consisting of
    several different terms, or components. Each term is a class inheriting from
    EquationOfStateTerms, not an instance of a class. It knows how to create
    instances when it needs and query them for values or create derivatives.
    """
    def __init__(self, Terms=[]):
        self._num_deriv_dx = 1e-8
        self._Terms = [ ]
        for t in Terms: self.add_term(t)

    def set_numerical_derivative_step(self, dx):
        """
        Sets the size of the step used for the numerical derivative
        estimates. The value is a fraction of the independent variable being
        stepped.
        """
        self._num_deriv_step = dx

    def add_term(self, term):
        """
        Add a term to the composite EOS.
        """
        self._Terms.append(term)

    def number_density(self, D, T, Y, derivative=None):
        return self._sample('number_density', derivative, D, T, Y)

    def pressure(self, D, T, Y, derivative=None):
        return self._sample('pressure', derivative, D, T, Y)

    def internal_energy(self, D, T, Y, derivative=None):
        return self._sample('internal_energy', derivative, D, T, Y)

    def specific_internal_energy(self, D, T, Y, derivative=None):
        return self._sample('specific_internal_energy', derivative, D, T, Y)

    def entropy(self, D, T, Y, derivative=None):
        return self._sample('entropy', derivative, D, T, Y)

    def gamma_effective(self, D, T, Y, method=2):
        """
        Computes the effective Gamma := (dlogp/dlogD)|_s to be used in the sound
        speed: c_s^2 = (\Gamma*p) / (D*h) using one of three different methods.
        """
        if method == 1:
            p = self.pressure(D,T,Y)
            s = self.entropy(D,T,Y)
            dpdD = self.pressure(D, T, Y, derivative='D')
            dpdT = self.pressure(D, T, Y, derivative='T')
            dsdD = self.entropy(D, T, Y, derivative='D')
            dsdT = self.entropy(D, T, Y, derivative='T')
            return (D/p)*(dpdD*dsdT - dpdT*dsdD) / dsdT

        elif method == 2:
            p = self.pressure(D,T,Y)
            dpdD = self.pressure(D, T, Y, derivative='D')
            dpdT = self.pressure(D, T, Y, derivative='T')
            dedD = self.specific_internal_energy(D, T, Y, derivative='D')
            dedT = self.specific_internal_energy(D, T, Y, derivative='T')
            return (D/p)*(dpdD*dedT - dpdT*dedD + p/D**2 * dpdT)/dedT

        elif method == 3:
            p = self.pressure(D,T,Y)
            dpdD = self.pressure(D, T, Y, derivative='D')
            dpdT = self.pressure(D, T, Y, derivative='T')
            dedD = self.specific_internal_energy(D, T, Y, derivative='D')
            dedT = self.specific_internal_energy(D, T, Y, derivative='T')

            return (D/p) * (dpdD*dedT + T/D**2 * dpdT**2) / dedT

        else:
            raise ValueError("Method must be either 1, 2, or 3")


    def _build_terms(self, *args):
        """
        This method should be over-ridden for building a composite EOS whose
        terms are order-dependent, or are not instantiated with D, T, and Y.
        """
        raise NotImplementedError("Derived class must implement this method.")


    def _sample(self, component, derivative, *args):
        """
        Private function, general handler for samples of the EOS components and
        their derivatives.
        """
        if not derivative:
            return sum([getattr(t, component)() for t in self._build_terms(*args)])

        elif len(derivative) == 1:
            f = getattr(self, component)
            n = {'D': 0, 'T': 1, 'Y': 2}[derivative]

            dx = self._num_deriv_dx
            X0, X1 = list(args), list(args)

            X1[n] += X1[n]*dx
            X0[n] -= X0[n]*dx

            return (f(*X1) - f(*X0)) / (X1[n] - X0[n])

        elif len(derivative) == 2:
            """
            Evaluating second derivatives. Unfinished.
            """
            raise NotImplementedError("Second derivative calculation unfinished.")
            f = getattr(self, component)
            n0 = {'D': 0, 'T': 1, 'Y': 2}[derivative[0]]
            n1 = {'D': 0, 'T': 1, 'Y': 2}[derivative[1]]
            
            dx = 1e-8
            X00, X01, X10, X11 = list(args), list(args), list(args), list(args)

        else:
            raise ValueError("Derivative string must be 0, 1, or 2 characters")



class IdealAdiabatic(EquationOfStateTerms):
    """
    Represents an adiabatic equation of state, with adiabatic constant 'gamma'.
    """
    def __init__(self, n, T, gamma=1.4):
        self.n = units.NumberDensity(n).convert_to('1/fm^3')
        self.T = units.Temperature(T).convert_to('MeV')
        self.gamma = gamma
        self._terms = { }
        self._set_terms()


    def _set_terms(self):
        g1 = self.gamma - 1.0
        p = self.n * self.T
        s = np.log(p/self.n**self.gamma) / g1

        f = self._terms
        f['n'] = units.NumberDensity([self.n, '1/fm^3'])
        f['p'] = units.EnergyDensity([p, 'MeV/fm^3'])
        f['u'] = units.EnergyDensity([p / g1, 'MeV/fm^3'])
        f['s'] = units.Entropy([s, 'MeV/K'])



class BlackbodyPhotons(EquationOfStateTerms):
    """
    Evaluates all terms for a photon gas. Only needs the temperature as a
    parameter.
    """
    def __init__(self, T):
        self.T = T
        self._terms = { }
        self._set_terms()

    def mass_density(self):
        return 0.0

    def _set_terms(self):
        """
        http://en.wikipedia.org/wiki/Photon_gas
        """
        z3 = 1.202056903159594285 # RiemannZeta(3), Mathematica
        T3 = np.power(self.T, 3)
        T4 = np.power(self.T, 4)
        a = pow(np.pi, 2) / (15*np.power(HBAR_C, 3))

        t = self._terms

        t['n'] = T3 * 2 * z3 / (np.power(np.pi, 2.0) * np.power(HBAR_C, 3.0))
        t['p'] = T4 * a / 3.0
        t['u'] = T4 * a
        t['s'] = (4./3.) * (T4 * a) / (self.T / BOLTZMANN_CONSTANT)



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
        self.T = T
        self._terms = { }
        self._set_terms()

    def _set_terms(self):
        Volume = np.power(np.pi, 2) * np.power(HBAR_C/self.T, 3)
        Energy = self.T

        f = fermion.neutrino_everything(self.mu/self.T)

        f['n'] *= (1.0 / Volume)
        f['p'] *= (Energy / Volume)
        f['u'] *= (Energy / Volume)
        f['s']  = (f['u'] + f['p']) / self.T - f['n'] * f['eta']

        for k in "npus":
            self._terms[k] = f[k]

        self._terms['eta'] = f['eta']



class FermionComponent(EquationOfStateTerms):
    """
    Represents an EOS component of electrons or positrons.
    """
    Volume = np.power(np.pi, 2) * np.power(HBAR_C/ELECTRON_MASS, 3)
    Energy = ELECTRON_MASS

    def __init__(self, mu, T):
        """
        Parameters:
        --------------------------------------------------------

        mu   : chemical potential, without rest-mass (MeV)
        T    : temperature (MeV)
        """
        self.mu = mu
        self.T = T
        self._terms = { }
        self._set_terms()

    def mass_density(self):
        """
        Returns the mass density of electrons in g/cm^3.
        """
        return self._terms['n'] * (ELECTRON_MASS/(LIGHT_SPEED*LIGHT_SPEED))

    def _set_terms(self):
        """
        Sets the number density, pressure, internal energy, and specific entropy
        (n,p,u,s) for both electrons and positrons.
        """
        eta = self.mu / self.T
        beta = self.T / ELECTRON_MASS
        f = fermion.fermion_everything(eta, beta)

        f['n'] *= (1.0 / self.Volume)
        f['p'] *= (self.Energy / self.Volume)
        f['u'] *= (self.Energy / self.Volume)
        f['s']  = (f['u'] + f['p']) / self.T - f['n'] * eta

        for k in "npus":
            self._terms[k] = f[k]

        self._terms['eta'] = eta



class FermiDiracElectrons(FermionComponent):
    """
    Evaluates the electron thermodynamics using exact Fermi-Dirac integrals.
    """
    def __init__(self, T, np):
        """
        Parameters:
        --------------------------------------------------------
        T    : temperature (MeV)
        np   : the number density of positively charged baryons (1/fm^3)
        """
        nu = fermion.solve_eta_pairs(T / ELECTRON_MASS, self.Volume * np)
        eta = +(nu - ELECTRON_MASS/T)
        super(FermiDiracElectrons, self).__init__(eta*T, T)



class FermiDiracPositrons(FermionComponent):
    """
    Evaluates the positron thermodynamics using exact Fermi-Dirac integrals.
    """
    def __init__(self, T, np):
        """
        Parameters:
        --------------------------------------------------------
        T    : temperature (MeV)
        np   : the number density of positively charged baryons (1/fm^3)
        """
        nu = fermion.solve_eta_pairs(self.T / ELECTRON_MASS, self.Volume * np)
        eta = -(nu + ELECTRON_MASS/T)
        super(FermiDiracPositrons, self).__init__(eta*T, T)



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

    def _set_terms(self):
        """
        http://scienceworld.wolfram.com/physics/ElectronDegeneracyPressure.html
        """
        num = np.power(np.pi, 2) * np.power(HBAR_C, 2)
        den = 5.0 * ELECTRON_MASS
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

    def _set_terms(self):
        """
        The Physical Universe: An Introduction to Astronomy, Frank H. Shu (1982)
        """
        self._terms['n'] = self.ne
        self._terms['p'] = 0.123 * 2*np.pi * HBAR_C * np.power(self.ne, 4./3.)
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
        self.D = D
        self.T = T
        self.Y = Y
        self._terms = { }
        self._set_terms()

    def chemical_potential(self, type):
        if type == 'neutrons': return self._terms['mu_n']
        if type == 'protons' : return self._terms['mu_p']
        raise ValueError("'type' must be either 'neutrons' or 'protons'")

    def _set_terms(self):
        if type(self)._table is None:
            cols = ['log10_rhoB', 'logT', 'Yp', 'p', 'nB', 'Eint',
                    'S', 'un', 'up']
            type(self)._table = shen.read_hdf5("data/eos3.h5", cols=cols)

        D, kT, Ye = self.D, self.T, self.Y
        t = self._terms
        e = type(self)._table
 
        t['n'] = shen.sample(e, 'nB'   , D, kT, Ye)
        t['p'] = shen.sample(e, 'p'    , D, kT, Ye)
        t['u'] = shen.sample(e, 'Eint' , D, kT, Ye) * t['n']
        t['s'] = shen.sample(e, 'S'    , D, kT, Ye) # entropy per baryon
        t['mu_n'] = shen.sample(e, 'un', D, kT, Ye) + 938.0
        t['mu_p'] = shen.sample(e, 'up', D, kT, Ye) + 938.0



if __name__ == "__main__":
    """
    Runs some tests.
    """
    import unittest
    print "testing", __file__

    class TestEquationOfStateTerms(unittest.TestCase):

        def test_instantiate(self):
            with self.assertRaises(AttributeError):
                eos = EquationOfStateTerms(1.0, 1.0, 1.0)


    class TestBlackbodyPhotons(unittest.TestCase):

        def test_instantiate(self):
            with self.assertRaises(TypeError):
                eos = BlackbodyPhotons(1.0)

        def test_eos(self):
            eos = BlackbodyPhotons(1e13, 40.0, 0.08)
            self.assertIsInstance(eos.pressure(), float)
            self.assertIsInstance(eos.entropy(), float)
            self.assertEqual(eos.mass_density(), 0.0)


    class TestColdElectrons(unittest.TestCase):

        def test_eos(self):
            eos = ColdElectrons(1e13, 40.0, 0.08)

            self.assertIsInstance(eos.pressure(), float)
            self.assertIsInstance(eos.mass_density(), float)

            with self.assertRaises(KeyError):
                s = eos.entropy()


    class TestNucleonsShenEos3(unittest.TestCase):

        def test_instantiate(self):
            eos = NucleonsShenEos3(1e13, 40.0, 0.08)
            self.assertIsNotNone(NucleonsShenEos3._table)

        def test_eos(self):
            eos = NucleonsShenEos3(1e13, 40.0, 0.08)
            self.assertIsInstance(eos.pressure(), float)
            self.assertIsInstance(eos.entropy(), float)

        def test_chemical_potential(self):
            eos = NucleonsShenEos3(1e13, 40.0, 0.08)
            with self.assertRaises(ValueError):
                eos.chemical_potential('dr. seuss')
            self.assertIsInstance(eos.chemical_potential('protons'), float)


    class TestEquationOfStateEvaluator(unittest.TestCase):
        
        def test_evaluate(self):
            pho = BlackbodyPhotons(1e13, 40.0, 0.08)
            ele = FermiDiracElectrons(1e13, 40.0, 0.08)
            eos = EquationOfStateEvaluator([BlackbodyPhotons, FermiDiracElectrons])
            self.assertEqual(eos.pressure(1e13, 40.0, 0.08),
                             pho.pressure() + ele.pressure())

        def test_derivatives(self):
            eos = EquationOfStateEvaluator([BlackbodyPhotons, FermiDiracElectrons])
            self.assertIsInstance(eos.pressure(1e13, 40.0, 0.08), float)
            self.assertIsInstance(eos.pressure(1e13, 40.0, 0.08, 'D'), float)
            self.assertIsInstance(eos.pressure(1e13, 40.0, 0.08, 'T'), float)
            self.assertIsInstance(eos.pressure(1e13, 40.0, 0.08, 'Y'), float)

    unittest.main()

