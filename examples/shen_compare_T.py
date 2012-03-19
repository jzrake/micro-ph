#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np


if __name__ == "__main__":

    D = 1e13
    T0, T1 = 10.0, 300.0
    Ye = 0.08
    temp = np.logspace(np.log10(T0), np.log10(T1), 60)

    def pressure(D, T, Y, Term):
        eos = Term(D, T, Y)
        return eos.pressure()

    for comp, tex, ls in [(FermiDiracElectrons, r"$e_-$", '-.'),
                          (FermiDiracPositrons, r"$e_+$", '-.'),
                          (NucleonsShenEos3, r"$n$ (shen)", '-'),
                          (BlackbodyPhotons, r"$\gamma$", '-x')]:

        print "Working out component '%s'" % comp
        p = [pressure(D, T, Ye, comp) for T in temp]
        plt.loglog(temp, p, ls, label=tex)

    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='lower right')
    plt.show()
