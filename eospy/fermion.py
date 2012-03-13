
"""
 * -----------------------------------------------------------------------------
 * FILE: fermion.py
 *
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 *
 * DESCRIPTION: Calculates contributions to the EOS from fermions.
 *
 *
 * REFERENCES:
 *
 * (+) Sekiguchi (2010)
 *     http://arxiv.org/abs/1009.3320
 *     NOTE: typo in eqn for pressure, missing 1/3
 *
 * (+) Timmes & Arnett (1999)
 *     http://iopscience.iop.org/0067-0049/125/1/277
 *     NOTE: typo in eqn (5), missing \beta in N_ele
 *
 * (+) Beaudet & Tassoul (1971)
 *     http://adsabs.harvard.edu/abs/1971A%26A....13..209B
 *
 * (+) Tooper (1969)
 *     http://adsabs.harvard.edu/full/1969ApJ...156.1075T
 *
 * (+) Kothari (1941)
 *     http://www.jstor.org/stable/97729
 *
 * NOTES:
 *
 * (+) Throughout the code comments, the letter h means 'h bar'
 *
 * -----------------------------------------------------------------------------
 """

from scipy.optimize import brentq
import numpy as np
import timmes.fdfunc
import fdfunc

# don't worry about overflows in exp
# np.seterr(over='ignore')

__all__ = ["evaluate_term",
           "solve_eta_pairs"]

FnBackend = "timmes"

def Fn(n, eta, beta):
    if FnBackend == "timmes":
        return timmes.fdfunc.dfermi(n, eta, beta)[0]
    elif FnBackend == "scipy":
        return fdfunc.dfermi(n, eta, beta)


def Fn_all(n, eta, beta):
    """
    Generates the Fermi-Dirac integral, as well as its first and second
    derivatives in eta, beta.

    Returns:
    --------------------------------------------------------

    f
    df/d(eta)
    df/d(beta)
    d^2(f)/d(eta)^2
    d^2(f)/d(beta)^2
    d^2(f)/(d(eta)d(theta))

    """
    return timmes.fdfunc.dfermi(n, eta, beta)


"""
These definitions are taken from Beaudet & Tassoul (1971).
"""
def fermion_number_density(eta, beta):
    return (1./2.) * np.power(2, 1.5) * np.power(beta, 1.5) * \
        (Fn(0.5, eta, beta) + 1.0*beta*Fn(1.5, eta, beta))

def fermion_pressure(eta, beta):
    return (1./3.) * np.power(2, 1.5) * np.power(beta, 2.5) * \
        (Fn(1.5, eta, beta) + 0.5*beta*Fn(2.5, eta, beta))

def fermion_internal_energy(eta, beta):
    return (1./2.) * np.power(2, 1.5) * np.power(beta, 2.5) * \
        (Fn(1.5, eta, beta) + 1.0*beta*Fn(2.5, eta, beta))


def dPdeta(eta, beta):
    F0, Fe, Fb,  = Fn_all(1.5, eta, beta)[:3]
    G0, Ge, Gb,  = Fn_all(2.5, eta, beta)[:3]
    return (1./3.) * np.power(2, 1.5) * np.power(beta, 2.5) * (Fe + 0.5*beta*Ge)


def dPdbeta(eta, beta):
    F0, Fe, Fb,  = Fn_all(1.5, eta, beta)[:3]
    G0, Ge, Gb,  = Fn_all(2.5, eta, beta)[:3]

    K = (1./3.) * np.power(2, 1.5)
    x = (5./2.) * np.power(beta, 1.5) * (F0 + 0.5*(0 + beta*G0))
    y = (1./1.) * np.power(beta, 2.5) * (Fb + 0.5*(G + beta*Gb))

    return K*x*y



def evaluate_term(key, sgn, eta, beta):
    """
    Evaluates the dimensionless EOS variable 'key', which is one of
    'number_density', 'pressure', or 'internal energy'.
    
    Parameters:
    --------------------------------------------------------

    eta  : mu/kT   ... chemical potential (including rest-mass)
    beta : kT/mc^2 ... relativity parameter


    Normalized Units:
    --------------------------------------------------------

    Volume : mi^2 (h/mc)^3
    Energy : mc^2

    """

    # For positrons, see TA99 eqn (5)
    if sgn < 0:
        eta = -eta - 2/beta

    terms = {
        "number_density": fermion_number_density,
        "pressure": fermion_pressure,
        "internal_energy": fermion_internal_energy }

    res = terms[key](eta, beta)

    if key == "internal_energy" and sgn < 0:
        # Correct for the self-energy of positrons, TA99 eqn (9)
        res += 2 * ELECTRON_MASS * fermion_number_density(eta, beta)

    return res


def solve_eta_pairs(beta, C):
    """
    Solves the implicit equation ne(e,b) - np(e,b) = C for e := eta, where C is
    the total number of positively charged baryons in the characteristic volume
    V0 := pi^2 (h/mc)^3.

    Parameters:
    --------------------------------------------------------

    beta : kT / mc^2    ... relativity parameter
    C    : n * V0       ... dimensionless number density
    """
    def f(eta):
        """ ne(e,b) - np(e,b) = C """
        return \
            evaluate_term("number_density", +1, eta, beta) - \
            evaluate_term("number_density", -1, eta, beta) - C

    bracket = 1.0
    while f(bracket) * f(-bracket) > 0.0: bracket *= 2.0

    return brentq(f, -bracket, bracket, disp=True)
