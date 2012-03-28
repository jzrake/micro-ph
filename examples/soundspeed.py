#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np



def TestGamma():
    D = 4.0
    T0, T1 = 0.2, 250.0
    Ye = 0.08
    temp = np.logspace(np.log10(T0), np.log10(T1), 40)


    eos = EquationOfStateEvaluator([IdealAdiabatic])
    gamma1 = [eos.gamma_effective(D,T,Ye,method=1)-1.4 for T in temp]
    gamma2 = [eos.gamma_effective(D,T,Ye,method=2)-1.4 for T in temp]
    gamma3 = [eos.gamma_effective(D,T,Ye,method=3)-1.4 for T in temp]

    plt.semilogx(temp, gamma1, '-o', lw=1.5, label='method 1')
    plt.semilogx(temp, gamma2, '-s', lw=1.5, label='method 2')
    plt.semilogx(temp, gamma3, '-x', lw=1.5, label='method 3')

    plt.title(r"Error in $\Gamma$ for an adiabatic equation of state")

    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"$\Delta \Gamma_{\rm{eff}}$", fontsize=16)
    plt.legend(loc='upper left')
    plt.show()


TestGamma()
