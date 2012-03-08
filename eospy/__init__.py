#!/usr/bin/env python

import numpy as np
from scipy.integrate import quadpack
from scipy.optimize import newton


# don't worry about overflows in exp
np.seterr(over='ignore')


def pdf_fermion(x, sgn, eta, beta):
    """
    Returns the integrand for fermion number density.

    Parameters:
    --------------------------------------------------------

    x     : pc/mc^2   ... independent variable
    eta   : mu/kT     ... degeneracy parameter
    beta  : mc^2/kT   ... unitless inverse temperature
    
    Notes
    --------------------------------------------------------

    Volume             ... pi^2 (h/mc)^3
    Energy             ... m c^2
    """
    return x*x / (np.exp(beta*(np.sqrt(1 + x*x) - sgn*1.0) - sgn*eta) + 1.0)

def ferm_number_density (*args): return 1.0
def ferm_pressure       (*args): return np.power(x,2) / (np.sqrt(1 + x*x))
def ferm_internal_energy(*args): return np.sqrt(1 + x*x) - args[1]*1



def pdf_fermion_massless(x, sgn, eta, beta):
    """
    Returns the integrand for fermion number density, in the massless limit.

    Parameters
    --------------------------------------------------------

    x     : pc/kT     ... independent variable
    eta   : mu/kT     ... degeneracy parameter
    beta  :           ... dummy argument

    Notes
    --------------------------------------------------------

    Volume             ... pi^2 (hc/kT)^3
    Energy             ... kT
    """
    return 0.5 * sgn * (x*x) / (1 + np.cosh(x - sgn*eta))

def fmml_number_density (*args): return 1.0
def fmml_pressure       (*args): return x
def fmml_internal_energy(*args): return x


def evaluate_term(key, sgn, eta, beta, massless=False):
    """
    Evaluates the EOS variable 'key', which is one of 'number_density',
    'pressure', or 'internal energy'.

    Notes
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
    return quadpack.quad(f, 0.0, quadpack.Inf)[0]



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
        return evaluate_n(+1, eta, beta) - evaluate_n(-1, eta, beta) - C
    return newton(f, 0.0)


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
