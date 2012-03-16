
import numpy as np
import shen
from fermion import *

__all__ = ["eos", "convert_beta_to_T"]



LIGHT_SPEED        = 2.997924580e+10 # cm/s
HBAR_C             = 1.973269718e+02 # MeV-fm
ELECTRON_MASS      = 5.110998928e-01 # MeV
ATOMIC_MASS_UNIT   = 9.314940612e+02 # MeV
PROTON_MASS        = 9.382720462e+02 # MeV
MEV_TO_ERG         = 1.602176487e-06
FM3_TO_CM3         = 1.000000000e-39
BOLTZMANN_CONSTANT = 8.617332400e-11 # MeV/K

ShenNucleonTable = { "table": None }


class EquationOfStateTerms(object):
    """
    Describes an EOS term.
    """
    def __init__(self, D, T, Y):
        self.D = D
        self.T = T
        self.Y = Y
        self._terms = { }

    def density(self):
        return self._terms['rho']

    def number_density(self):
        return self._terms['n']

    def pressure(self):
        return self._terms['p']

    def internal_energy(self):
        return 0.0

    def entropy(self):
        return self._terms['s']

    def chemical_potential(self):
        return self._terms['eta'] * self.T

    def consistency(self):
        pass


class BlackbodyPhotons(EquationOfStateTerms):
    """
    http://en.wikipedia.org/wiki/Photon_gas
    """
    def _set_terms(self):
        z3 = 1.202056903159594285 # RiemannZeta(3), Mathematica
        T3 = np.power(self.T, 3)
        T4 = np.power(self.T, 4)

        self.terms = { }

        self.terms['n'] = T3 * 2 * z3 / (np.power(np.pi, 2.0) * np.power(HBAR_C, 3.0))
        self.terms['p'] = T4 * a / 3.0
        self.terms['u'] = T4 * a
        self.terms['s'] = (4./3.) * self.terms['u'] / self.T # entropy has no units here



class ElectronPositronPairs(EquationOfStateTerms):
    """
    
    """
    def _set_terms(self):
        fermion.solve_eta_pairs()


    def eval_pairs(self):
        """
        Sets the number density, pressure, and internal (n,p,u) energy for both
        electrons and positrons.
        
        Parameters:
        --------------------------------------------------------
        
        D   : density (g/cm^3)
        kT  : temperature (MeV)
        Ye  : proton/electron fraction

        """

        D, kT, Ye = self.D, self.T, self.Y

        c2 = LIGHT_SPEED*LIGHT_SPEED
        Volume = np.power(np.pi, 2) * np.power(HBAR_C/ELECTRON_MASS, 3)
        Energy = ELECTRON_MASS

        Erest = D * c2 * FM3_TO_CM3 / MEV_TO_ERG
        C = Volume * Ye * Erest / ATOMIC_MASS_UNIT

        beta = kT / ELECTRON_MASS
        eta = solve_eta_pairs(beta, C)

        self.electrons = fermion_everything(+1, eta, beta)
        self.positions = fermion_everything(-1, eta, beta)

        return self.electrons, self.positrons



def eos(D, kT, Ye, component):
    """
    Returns the number density, pressure, and internal energy density given the
    density, temperature, and internal energy.

    Parameters:
    --------------------------------------------------------

    component : A term or list of terms to be included. May contain:

    ["electrons", "positrons", "photons", "cold_electrons", "dense_electrons"]

    D   : density (g/cm^3)
    kT  : temperature (MeV)
    Ye  : proton/electron fraction

    """
    if type(component) is str: component = [component]

    n = 0.0
    p = 0.0
    u = 0.0
    mu = 0.0

    if "electrons" in component:
        ev = eval_pairs(D, kT, Ye, +1)
        n += ev[0]
        p += ev[1]
        u += ev[2]
        mu += ev[3]

    if "positrons" in component:
        ev = eval_pairs(D, kT, Ye, -1)
        n += ev[0]
        p += ev[1]
        u += ev[2]
        mu += ev[3]

    if "photons" in component:
        a = pow(np.pi, 2) / (15*np.power(HBAR_C, 3))
        n += 0.0
        p += a * np.power(kT, 4) / 3.0
        u += a * np.power(kT, 4)
        mu += 0.0

    if "cold_electrons" in component:
        """
        http://scienceworld.wolfram.com/physics/ElectronDegeneracyPressure.html
        """
        Erest = D * LIGHT_SPEED*LIGHT_SPEED * FM3_TO_CM3 / MEV_TO_ERG
        ne = Ye * Erest / ATOMIC_MASS_UNIT

        num = np.power(np.pi, 2) * np.power(HBAR_C, 2)
        den = 5.0 * ELECTRON_MASS
        las = np.power(3.0/np.pi, 2./3.) * np.power(ne, 5./3.)

        n += ne
        p += (num/den) * las
        u += 0.0 # not implemented yet
        mu += 0.0 # not implemented yet

    if "dense_electrons" in component:
        Erest = D * LIGHT_SPEED*LIGHT_SPEED * FM3_TO_CM3 / MEV_TO_ERG
        ne = Ye * Erest / ATOMIC_MASS_UNIT

        n += ne
        p += 0.123 * 2*np.pi * HBAR_C * np.power(ne, 4./3.)
        u += 0.0 # not implemented yet
        mu += 0.0 # not implemented yet

    if "nucleons" in component:

        if ShenNucleonTable["table"] is None:
            cols = ['log10_rhoB', 'logT', 'Yp', 'p', 'nB', 'Eint']
            ShenNucleonTable["table"] = shen.read_hdf5("data/eos3.h5", cols=cols)

        table = ShenNucleonTable["table"]

        n += shen.sample(table, 'nB'  , D, kT, Ye)
        p += shen.sample(table, 'p'   , D, kT, Ye)
        u += shen.sample(table, 'Eint', D, kT, Ye)
        mu += 0.0

    return n, p, u, mu


def convert_beta_to_kT(beta):
    """
    Convenience function, returns the temperature (in MeV) given beta := kT /
    mc^2.
    """
    return beta / ELECTRON_MASS


def convert_kT_to_Kelvin(kT):
    return kT / BOLTZMANN_CONSTANT
