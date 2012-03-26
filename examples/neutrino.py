#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np


def pressure(D, T, Y, verbose=False):

    electron = FermiDiracElectrons(D, T, Y)
    positron = FermiDiracPositrons(D, T, Y)
    nucleons = NucleonsShenEos3(D, T, Y)

    mu_n = nucleons.chemical_potential('neutrons')
    mu_p = nucleons.chemical_potential('protons')
    mu_e = electron.chemical_potential()
    mu_t = positron.chemical_potential()

    if verbose:
        print "chemical potentials:"
        print "e-:", mu_e
        print "e+:", mu_t
        print "p :", mu_p
        print "n :", mu_n

    mu_nu_p = mu_e + mu_p - mu_n
    mu_nu_n = mu_t + mu_n - mu_p

    if verbose:
        print "nu+:", mu_nu_p
        print "nu-:", mu_nu_n

    neutrino = NeutrinoComponent(+1, mu_nu_p, T)
    aeutrino = NeutrinoComponent(-1, mu_nu_n, T)

    if verbose:
        print "pressures:"
        print "e- :", electron.pressure()
        print "e+ :", positron.pressure()
        print "p+n:", nucleons.pressure()
        print "nu+:", neutrino.pressure()
        print "nu-:", aeutrino.pressure()

    p = [c.pressure() for c in electron, positron, nucleons, neutrino, aeutrino]
    return p


D = 1e13
T0, T1 = 0.1, 200.0
Ye = 0.08
temp = np.logspace(np.log10(T0), np.log10(T1), 140)


result = [pressure(D, T, Ye, verbose=True) for T in temp]

tex = [r"$e_-$", r"$e_+$", r"$n$ (shen)", r"$\nu_e$", r"$\bar{\nu}_e$"]
ls = ['-.', '-.', '-', ':', '--']

for i in range(5):
    plt.loglog(temp, [r[i] for r in result], ls[i], lw=2.0, label=tex[i])


plt.title(r"Various EOS components for $\rho=10^{13} \ \rm{g/cm^3}$")
plt.ylim(1e-12, 10000.0)
plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
plt.legend(loc='lower right')
plt.show()


