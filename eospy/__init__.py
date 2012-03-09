#!/usr/bin/env python

import numpy as np
from scipy.integrate import quadpack
from scipy.optimize import newton, brentq


LIGHT_SPEED      = 2.997924580e+10 # cm/s
HBAR_C           = 1.973269718e+02 # MeV-fm
ELECTRON_MASS    = 5.110998928e-01 # MeV
ATOMIC_MASS_UNIT = 9.314940612e+02 # MeV
PROTON_MASS      = 9.382720462e+02 # MeV
MEV_TO_ERG       = 1.602176487e-06
FM3_TO_CM3       = 1.000000000e-39


# don't worry about overflows in exp
np.seterr(over='ignore')

#import shen

def pdf_fermion(x, sgn, eta, beta):
    """
    Returns the integrand for fermion number density.

    Parameters:
    --------------------------------------------------------

    x     : pc/mc^2   ... independent variable
    eta   : mu/kT     ... degeneracy parameter
    beta  : mc^2/kT   ... unitless inverse temperature
    
    Notes:
    --------------------------------------------------------

    Volume             ... pi^2 (h/mc)^3
    Energy             ... m c^2
    """
    return x*x / (np.exp(beta*(np.sqrt(1 + x*x) - sgn*1.0) - sgn*eta) + 1.0)

def ferm_number_density (x,s,e,b): return 1.0
def ferm_pressure       (x,s,e,b): return np.power(x,2) / (np.sqrt(1 + x*x))
def ferm_internal_energy(x,s,e,b): return np.sqrt(1 + x*x) - s*1



def pdf_fermion_massless(x, sgn, eta, beta):
    """
    Returns the integrand for fermion number density, in the massless limit.

    Parameters:
    --------------------------------------------------------

    x     : pc/kT     ... independent variable
    eta   : mu/kT     ... degeneracy parameter
    beta  :           ... dummy argument

    Notes:
    --------------------------------------------------------

    Volume             ... pi^2 (hc/kT)^3
    Energy             ... kT
    """
    return 0.5 * sgn * (x*x) / (1 + np.cosh(x - sgn*eta))

def fmml_number_density (x,s,e,b): return 1.0
def fmml_pressure       (x,s,e,b): return x
def fmml_internal_energy(x,s,e,b): return x


def evaluate_term(key, sgn, eta, beta, massless=False):
    """
    Evaluates the EOS variable 'key', which is one of 'number_density',
    'pressure', or 'internal energy'.

    Notes:
    --------------------------------------------------------

    Use sgn = +1/-1 for particles and anti-particles respecticely. Eta is always
    the 'degeneracy parameter', i.e. the chemical potential over kT. Beta is
    inverse temperature, where the dimensions depend on whether the particle is
    massless or not.
    """
    massive_terms = {
        "pdf": pdf_fermion,
        "number_density": ferm_number_density,
        "pressure": ferm_pressure,
        "internal_energy": ferm_internal_energy }

    massless_terms = {
        "pdf": pdf_fermion_massless,
        "number_density": fmml_number_density,
        "pressure": fmml_pressure,
        "internal_energy": fmml_internal_energy }

    terms = massless_terms if massless else massive_terms
    pdf, val = terms["pdf"], terms[key]
    f = lambda x: val(x, sgn, eta, beta) * pdf(x, sgn, eta, beta)

    res = quadpack.quad(f, 0.0, quadpack.Inf)
    return res[0]


def solve_eta_pairs(beta, C):
    """
    Solves the implicit equation ne(e,b) - np(e,b) = C for e := eta, where C is
    the total number of positively charged baryons in the characteristic volume
    V0 := pi^2 (h/mc)^3.

    Parameters:
    --------------------------------------------------------

    beta : mc^2 / kT    ... dimensionless inverse temperature
    C    :              ... dimensionless number density
    """
    def f(eta):
        """
        ne(e,b) - np(e,b) = C
        """
        return \
            evaluate_term("number_density", +1, eta, beta) - \
            evaluate_term("number_density", -1, eta, beta) - C

    bracket = 1.0
    while f(bracket) * f(-bracket) > 1.0: bracket *= 2.0

    return brentq(f, -bracket, bracket, disp=True)



def solve_eta_neutrinos(C):
    """
    Solves the implicit equation n(e) = C for e := eta, where C the total
    neutrino number in the characteristic volume V0 := pi^2 (hc/kT)^3.

    Parameters:
    --------------------------------------------------------

    C    :              ... dimensionless number density
    """
    def f(eta):
        return evaluate_n(+1, eta, beta) - C
    return newton(f, 0.0)



def eval_pairs(D, kT, Ye, sgn):
    """
    Returns the number density, pressure, and internal (n,p,u) energy for
    either electrons or positrons depending on sgn (+/-).

    Parameters:
    --------------------------------------------------------

    D   : density (g/cm^3)
    kT  : temperature (MeV)
    Ye  : proton/electron fraction
    sgn : (+) for electrons, (-) for positrons

    Returns:
    --------------------------------------------------------

    n   : number density (1/fm^3)
    p   : pressure (MeV/fm^3)
    u   : internal energy (MeV/fm^3)
    """
    c2 = LIGHT_SPEED*LIGHT_SPEED
    Volume = np.power(np.pi, 2) * np.power(HBAR_C/ELECTRON_MASS, 3)
    Energy = ELECTRON_MASS

    Erest = D * c2 * FM3_TO_CM3 / MEV_TO_ERG
    C = Volume * Ye * Erest / ATOMIC_MASS_UNIT

    beta = ELECTRON_MASS / kT
    eta = solve_eta_pairs(beta, C)

    nk = "number_density"
    pk = "pressure"
    uk = "internal_energy"

    n = (1.0    / Volume) * evaluate_term(nk, sgn, eta, beta)
    p = (Energy / Volume) * evaluate_term(pk, sgn, eta, beta)
    u = (Energy / Volume) * evaluate_term(uk, sgn, eta, beta)

    return n,p,u



def eos_eval(D, kT, Ye, component):
    """
    Returns the number density, pressure, and internal energy given the density,
    temperature, and internal energy.

    Parameters:
    --------------------------------------------------------

    component : A term or list of terms to be included. May contain:

    ["electrons", "positrons", "photons", "neutrinos"]

    D   : density (g/cm^3)
    kT  : temperature (MeV)
    Ye  : proton/electron fraction
    """
    if type(component) is str: component = [component]

    n = 0.0
    p = 0.0
    u = 0.0

    if "electrons" in component:
        ev = eval_pairs(D, kT, Ye, +1)
        n += ev[0]; p += ev[1]; u += ev[2]

    if "positrons" in component:
        ev = eval_pairs(D, kT, Ye, -1)
        n += ev[0]; p += ev[1]; u += ev[2]

    if "photons" in component:
        a = pow(np.pi, 2) / (15*np.power(HBAR_C, 3))
        n += 0.0
        p += a * np.power(kT, 4) / 3.0
        u += a * np.power(kT, 4)

    if "cold_electrons" in component:
        # http://scienceworld.wolfram.com/physics/ElectronDegeneracyPressure.html

        Erest = D * LIGHT_SPEED*LIGHT_SPEED * FM3_TO_CM3 / MEV_TO_ERG
        ne = Ye * Erest / ATOMIC_MASS_UNIT

        # adding the factor 1/8pi below brings cold electrons very close to the
        # limiting case of real electrons

        num = np.power(np.pi, 2) * np.power(HBAR_C, 2) / (8*np.pi)
        den = 5.0 * ELECTRON_MASS
        las = np.power(3.0/np.pi, 2./3.) * np.power(ne, 5./3.)

        n += ne
        p += (num/den) * las
        u += 0.0 # not implemented yet

    return n, p, u


def convert_beta_to_T(beta):
    """
    Convenience function, returns the temperature (in MeV) given beta :=
    mc^2/kT.
    """
    return ELECTRON_MASS / beta
