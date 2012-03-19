
import numpy as np
import shen
import fermion



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
    Describes an EOS term.
    """
    def __init__(self, D, T, Y):
        self.D = D
        self.T = T
        self.Y = Y
        self._terms = { }
        self._set_terms()

    def number_density(self):
        return self._terms['n']

    def pressure(self):
        return self._terms['p']

    def internal_energy(self):
        return self._terms['u']

    def entropy_density(self):
        return self._terms['s']

    def chemical_potential(self):
        return self._terms['eta'] * self.T

    def mass_density(self):
        return self.D

    def consistency(self):
        pass


class BlackbodyPhotons(EquationOfStateTerms):
    """
    Evaluates all terms for a photon gas. Only needs the temperature as a
    parameter.
    """
    def __init__(self, T):
        self.D = 0.0
        self.T = T
        self._terms = { }
        self._set_terms()

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


class FermionComponent(EquationOfStateTerms):
    """
    Represents an EOS component of electrons and/or positrons.

    Notes:
    ----------------------------------------------------------------------------

    Objects inheriting from this class are instantiated with the baryon mass
    density and proton fraction as 'D'. However, the mass_density method returns
    the mass density of fermions.
    """
    def mass_density(self):
        """
        Returns the mass density of electrons in g/cm^3.
        """
        return self._terms['n'] * (ELECTRON_MASS/(LIGHT_SPEED*LIGHT_SPEED))


class ElectronPositronPairs(FermionComponent):
    """
    Evaluates the electron and positron pressure using exact Fermi-Dirac
    integrals.
    """
    def _set_terms(self):
        """
        Sets the number density, pressure, and internal (n,p,u) energy for both
        electrons and positrons.
        
        Parameters (taken from self):
        --------------------------------------------------------
        
        D  : density (g/cm^3)
        T  : temperature (MeV)
        Y  : proton/electron fraction

        """

        c2 = LIGHT_SPEED*LIGHT_SPEED
        Volume = np.power(np.pi, 2) * np.power(HBAR_C/ELECTRON_MASS, 3)
        Energy = ELECTRON_MASS

        Erest = self.D * c2 * FM3_TO_CM3 / MEV_TO_ERG
        C = Volume * self.Y * Erest / ATOMIC_MASS_UNIT

        beta = self.T / ELECTRON_MASS
        eta = solve_eta_pairs(beta, C)

        ele = fermion_everything(+1, eta, beta)
        pos = fermion_everything(-1, eta, beta)
    
        for f in [ele, pos]:
            f['n'] *= (1.0 / Volume)
            f['p'] *= (Energy / Volume)
            f['u'] *= (Energy / Volume)
            f['s'] *= (BOLTZMANN_CONSTANT / Volume)

        for k in "npus":
            self._terms[k] = ele[k] + pos[k]

        self._terms['eta'] = eta
        self.electrons = ele
        self.positrons = pos



class ColdElectrons(FermionComponent):
    """
    Evaluates the thermodynamic variables for a cold (fully degenerate) electron
    gas. Accurate as long as the Fermi momentum is not relativistic.
    """
    def _set_terms(self):
        """
        http://scienceworld.wolfram.com/physics/ElectronDegeneracyPressure.html
        """
        Erest = self.D * LIGHT_SPEED*LIGHT_SPEED * FM3_TO_CM3 / MEV_TO_ERG
        ne = self.Y * Erest / ATOMIC_MASS_UNIT

        num = np.power(np.pi, 2) * np.power(HBAR_C, 2)
        den = 5.0 * ELECTRON_MASS
        las = np.power(3.0/np.pi, 2./3.) * np.power(ne, 5./3.)

        self._terms['n'] = ne
        self._terms['p'] = (num/den) * las
        # 'u' not implemented yet
        # 's' not implemented yet



class DenseElectrons(FermionComponent):
    """
    Evaluates the thermodynamic variables for a dense and fully degenerate
    electron gas, where the Fermi momentum is ultra-relativistic.
    """
    def _set_terms(self):
        """
        The Physical Universe: An Introduction to Astronomy, Frank H. Shu (1982)
        """
        Erest = self.D * LIGHT_SPEED*LIGHT_SPEED * FM3_TO_CM3 / MEV_TO_ERG
        ne = self.Y * Erest / ATOMIC_MASS_UNIT

        self._terms['n'] = ne
        self._terms['p'] = 0.123 * 2*np.pi * HBAR_C * np.power(ne, 4./3.)
        # 'u' not implemented yet
        # 's' not implemented yet



class NucleonsShenEos3(EquationOfStateTerms):
    """
    Evaluates the thermodynamic variables for dense baryons using the (updated)
    lookup table 'eos3' of Shen et. al. (1998).

    user guide: http://user.numazu-ct.ac.jp/~sumi/eos/table2/guide_EOS3.pdf
    full table: http://user.numazu-ct.ac.jp/~sumi/eos/table2/eos3.tab.gz
    """

    _table = None # caches the hdf5 table once it is loaded

    def _set_terms(self):
        if type(self)._table is None:
            cols = ['log10_rhoB', 'logT', 'Yp', 'p', 'nB', 'Eint', 'S']
            type(self)._table = shen.read_hdf5("data/eos3.h5", cols=cols)

        D, kT, Ye = self.D, self.T, self.Y
        t = self._terms
        e = type(self)._table
 
        t['n'] = shen.sample(e, 'nB'  , D, kT, Ye)
        t['p'] = shen.sample(e, 'p'   , D, kT, Ye)
        t['u'] = shen.sample(e, 'Eint', D, kT, Ye)
        t['s'] = shen.sample(e, 'S'   , D, kT, Ye) * t['n']



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
                eos = BlackbodyPhotons(1.0, 1.0, 1.0)

        def test_eos(self):
            eos = BlackbodyPhotons(40.0)
            self.assertIsInstance(eos.pressure(), float)
            self.assertIsInstance(eos.entropy_density(), float)
            self.assertEqual(eos.mass_density(), 0.0)


    class TestColdElectrons(unittest.TestCase):

        def test_eos(self):
            eos = ColdElectrons(1e13, 40.0, 0.08)

            self.assertIsInstance(eos.pressure(), float)
            self.assertIsInstance(eos.mass_density(), float)

            with self.assertRaises(KeyError):
                s = eos.entropy_density()


    class TestNucleonsShenEos3(unittest.TestCase):

        def test_instantiate(self):
            eos = NucleonsShenEos3(1e13, 40.0, 0.08)
            self.assertIsNotNone(NucleonsShenEos3._table)

        def test_eos(self):
            eos = NucleonsShenEos3(1e13, 40.0, 0.08)
            self.assertIsInstance(eos.pressure(), float)


    unittest.main()



