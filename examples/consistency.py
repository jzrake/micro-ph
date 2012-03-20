#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np


def PlotConsistency(eos):
    """
    Checks the thermodynamic consistency relation:

    p = rho^2 d(u)/d(rho) + T dp/dT
    """

    def consistency(D, T, Y):
        p    = eos.pressure(D, T, Y)
        dpdT = eos.pressure(D, T, Y, derivative='T')
        dedD = eos.specific_internal_energy(D, T, Y, derivative='D')

        return abs(1.0 - (D**2 * dedD + T*dpdT)/p)

    temp = np.logspace(5.0, 10.0, 3) # in Kelvins
    dens = np.logspace(-3, 14, 10)
    Ye = 0.08

    for T in temp:
        kT = T * BOLTZMANN_CONSTANT
        c = np.array([consistency(D, kT, Ye) for D in dens])
        plt.loglog(dens, c, '-o', label=r"$T=10^{%d}\rm{K}$" % np.log10(T))

    plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
    plt.ylabel(r"$1 - \frac{\rho^2}{p} \frac{\partial E}{\partial \rho} - " +
               r"\frac{T}{p} \frac{\partial p}{\partial T}$", fontsize=18)
    plt.legend(loc='lower right')
    plt.show()



if __name__ == "__main__":
    """
    Checks the thermodynamic consistency relation:

    p = rho^2 d(u)/d(rho) + T dp/dT
    """

    eos = EquationOfStateEvaluator([FermiDiracElectrons, FermiDiracPositrons])
    #eos = EquationOfStateEvaluator([NucleonsShenEos3])
    PlotConsistency(eos)
