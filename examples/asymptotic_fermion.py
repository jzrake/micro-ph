#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np



"""
Plots the pressure for various components of the EOS as a function of
density.
"""
temp = np.logspace(5.0, 10.0, 3) # in Kelvins
dens = np.logspace(-3, 14, 100)
Ye = 0.08

def pressure(D, T, Y, Term):
    eos = Term(D, T, Y)
    return eos.pressure()

for T in temp:
    kT = T * BOLTZMANN_CONSTANT
    p = [pressure(D, kT, Ye, FermiDiracElectrons) for D in dens]
    plt.loglog(dens, p, '-', lw=0.8, label=r"$T=10^{%d}\rm{K}$" % np.log10(T))

cold  = [pressure(D, kT, Ye, ColdElectrons) for D in dens]
plt.loglog(dens, cold, c='k', ls=':', label=r"cold $e_-$ $\propto n_e^{5/3}$")

dense = [pressure(D, kT, Ye, DenseElectrons) for D in dens]
plt.loglog(dens, dense, c='r', ls='--', lw=3.0,
           label=r"dense $e_-$ $\propto n_e^{4/3}$")

plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
plt.legend(loc='lower right')
plt.show()

