#!/usr/bin/env python

import numpy
import quantities as pq
from matplotlib import pyplot as plt
from eospy.physics import *
from eospy import shen


"""
Compares the pressure contribution from many different components, to that
of Shen nucleons.
"""
D = 1e13 * pq.g/pq.cm**3
Yp = 0.08
np = Yp * D / shen.amu
T0, T1 = 4.0, 200.0
temp = numpy.logspace(numpy.log10(T0), numpy.log10(T1), 80) * pq.MeV


ele = [FermiDiracElectrons(np, T) for T in temp]
pos = [FermiDiracPositrons(np, T) for T in temp]
pho = [BlackbodyPhotons(T) for T in temp]
nuc = [NucleonsShenEos3(D, T, Yp) for T in temp]

for comp, tex, ls in [(ele, r"$e_-$", '-.'),
                      (pos, r"$e_+$", '-.'),
                      (nuc, r"$n$ (shen)", '-'),
                      (pho, r"$\gamma$", '-x')]:

    print "Working out component '%s'" % tex
    p = [c.pressure().rescale('MeV/fm^3') for c in comp]
    plt.loglog(temp, p, ls, lw=2.0, label=tex)

plt.title(r"Various EOS components for $\rho=10^{13} \ \rm{g/cm^3}$")
plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
plt.legend(loc='lower right')
plt.show()
