#!/usr/bin/env python

import numpy
import quantities as pq
from matplotlib import pyplot as plt
from eospy.physics import *


"""
Plots the pressure for various components of the EOS as a function of
density.
"""
Yp = 0.08
temp = numpy.logspace( 5.0, 10.0,   3) * pq.K
dens = numpy.logspace(-3.0, 14.0, 100) * pq.g/pq.cm**3

ele_true = lambda T: [FermiDiracElectrons(Yp*n/shen.amu, T) for n in dens]
ele_cold = [ColdElectrons(Yp*n/shen.amu).pressure().rescale('MeV/fm^3') for n in dens]
ele_dens = [DenseElectrons(Yp*n/shen.amu).pressure().rescale('MeV/fm^3') for n in dens]

plt.loglog(dens, ele_cold, c='k', ls=':', label=r"cold $e_-$ $\propto n_e^{5/3}$")
plt.loglog(dens, ele_dens, c='r', ls='--',
           lw=3.0, label=r"dense $e_-$ $\propto n_e^{4/3}$")

for T in temp:
    p = [c.pressure() for c in ele_true(T)]
    plt.loglog(dens, p, '-', lw=0.8,
               label=r"$T=10^{%d}\rm{K}$" % np.log10(T.magnitude))

plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
plt.legend(loc='lower right')
plt.show()

