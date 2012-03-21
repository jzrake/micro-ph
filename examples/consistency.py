#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np


def PlotConsistency(eos, T0, T1, D0, D1, Ye=0.08, Kelvin=True, N=60):
    """
    Checks the thermodynamic consistency relation:

    p = rho^2 d(u)/d(rho) + T dp/dT
    """

    def consistency(D, T, Y):
        p    = eos.pressure(D, T, Y)
        dpdT = eos.pressure(D, T, Y, derivative='T')
        dedD = eos.specific_internal_energy(D, T, Y, derivative='D')

        return abs(1.0 - (D**2 * dedD + T*dpdT)/p)

    temp = np.logspace(T0, T1, 3)
    dens = np.logspace(D0, D1, N)
    Ye = 0.08

    for T in temp:
        kT = T * BOLTZMANN_CONSTANT if Kelvin else T
        c = np.array([consistency(D, kT, Ye) for D in dens])
        lab = r"$T=10^{%d}\rm{K}$" % np.log10(T) if Kelvin \
            else r"$T=%2.1f\rm{MeV}$" % T
        plt.loglog(dens, c, '-o', label=lab)

    plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
    plt.ylabel(r"$1 - \frac{\rho^2}{p} \frac{\partial E}{\partial \rho} - " +
               r"\frac{T}{p} \frac{\partial p}{\partial T}$", fontsize=18)
    plt.legend(loc='best')
    plt.show()


#eos = EquationOfStateEvaluator([FermiDiracElectrons, FermiDiracPositrons])
eos = EquationOfStateEvaluator([BlackbodyPhotons])
PlotConsistency(eos, 5, 11, -3, 14)

#eos = EquationOfStateEvaluator([NucleonsShenEos3])
#eos.set_numerical_derivative_step(1e-4)
#PlotConsistency(eos, np.log10(0.5), np.log10(260.0), 6, 15, Kelvin=False, N=5)

