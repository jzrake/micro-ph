
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

__all__ = ["solve_eta_pairs",
           "fermion_everything"]

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


def fermion_everything(sgn, eta, beta):
    """
    Evaluates the dimensionless number density, pressure, and internal energy,
    and entropy.

    Parameters:
    --------------------------------------------------------

    eta  : mu/kT   ... chemical potential (including rest-mass)
    beta : kT/mc^2 ... relativity parameter


    Returns:
    --------------------------------------------------------

    A dictionary whose keys are self-explanatory.


    Notes:
    --------------------------------------------------------

    Return values are normalized with the following units:

    Volume : mc^2 (h/mc)^3
    Energy : mc^2

    The formulae below are taken from Beaudet & Tassoul (1971). The idea of
    including the rest-mass in the chemical potential is inspired by TA99.

    """
    # For positrons, see TA99 eqn (5)
    if sgn < 0:
        eta = -eta - 2/beta

    F, Fe, Fb = Fn_all(0.5, eta, beta)[:3]
    G, Ge, Gb = Fn_all(1.5, eta, beta)[:3]
    H, He, Hb = Fn_all(2.5, eta, beta)[:3]

    t15 = np.power(2, 1.5)
    B05 = np.power(beta, 0.5)
    B15 = np.power(beta, 1.5)
    B25 = np.power(beta, 2.5)

    res = { }

    res['n'] = (1./2.) * t15 * B15 * (F + 1.0*beta*G)
    res['p'] = (1./3.) * t15 * B25 * (G + 0.5*beta*H)
    res['u'] = (1./2.) * t15 * B25 * (G + 1.0*beta*H)
    res['s'] = (res['u'] + res['p'])/beta - res['n'] * eta

    # Correct for the self-energy of positrons, TA99 eqn (9)
    if sgn < 0:
        res['u'] += 2*res['n']

    return res


def electron_everything(sgn, eta, beta):
    """
    Like fermion_everything, but also generates derivatives of n, p, and u. Only
    works for electrons right now.
    """
    assert(sgn > 0)

    F, Fe, Fb = Fn_all(0.5, eta, beta)[:3]
    G, Ge, Gb = Fn_all(1.5, eta, beta)[:3]
    H, He, Hb = Fn_all(2.5, eta, beta)[:3]

    t15 = np.power(2, 1.5)
    B05 = np.power(beta, 0.5)
    B15 = np.power(beta, 1.5)
    B25 = np.power(beta, 2.5)

    res = { }

    res['n'] = (1./2.) * t15 * B15 * (F + 1.0*beta*G)
    res['p'] = (1./3.) * t15 * B25 * (G + 0.5*beta*H)
    res['u'] = (1./2.) * t15 * B25 * (G + 1.0*beta*H)

    res['dndeta'] = (1./2.) * t15 * B15 * (Fe + 1.0*beta*Ge)
    res['dpdeta'] = (1./3.) * t15 * B25 * (Ge + 0.5*beta*He)
    res['dudeta'] = (1./2.) * t15 * B25 * (Ge + 1.0*beta*He)

    # number density
    x = (3./2.) * B05 * (F  + (0 + beta*G ))
    y = (1./1.) * B15 * (Fb + (G + beta*Gb))
    res['dndbeta'] = (1./2.) * t15 * (x + y)

    # pressure
    x = (5./2.) * B15 * (G  + 0.5*(0 + beta*H ))
    y = (1./1.) * B25 * (Gb + 0.5*(H + beta*Hb))
    res['dpdbeta'] = (1./3.) * t15 * (x + y)

    # inernal energy
    x = (5./2.) * B15 * (G  + (0 + beta*H ))
    y = (1./1.) * B25 * (Gb + (H + beta*Hb))
    res['dudbeta'] = (1./2.) * t15 * (x + y)

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
            fermion_everything(+1, eta, beta)['n'] - \
            fermion_everything(-1, eta, beta)['n'] - C

    bracket = 1.0
    while f(bracket) * f(-bracket) > 0.0: bracket *= 2.0

    return brentq(f, -bracket, bracket, disp=True)
