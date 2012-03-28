#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np



def TestGamma():
    D = 4.0
    T0, T1 = 0.2, 250.0
    Ye = 0.08
    temp = np.logspace(np.log10(T0), np.log10(T1), 4)


    eos = EquationOfStateEvaluator([IdealAdiabatic])
    #gamma1 = [eos.gamma_effective(D,T,Ye,method=1) for T in temp]
    gamma2 = [eos.gamma_effective(D,T,Ye,method=2) for T in temp]
    gamma3 = [eos.gamma_effective(D,T,Ye,method=3) for T in temp]

    #plt.semilogx(temp, gamma1, lw=1.5, label='method 1')
    plt.semilogx(temp, gamma2, lw=1.5, label='method 2')
    plt.semilogx(temp, gamma3, lw=1.5, label='method 3')

    print '2', gamma2
    print '3', gamma3

    plt.title(r"Various EOS components for $\rho=10^{13} \ \rm{g/cm^3}$")
    #plt.ylim(1e-12, 10000.0)
    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='lower right')
    plt.show()


TestGamma()
