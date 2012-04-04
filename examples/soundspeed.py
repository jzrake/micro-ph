#!/usr/bin/env python

from eospy.physics import *
from eospy import eos
from matplotlib import pyplot as plt
from matplotlib import ticker
import numpy as np
import quantities as pq


def TestGamma():
    D = 4.0
    T0, T1 = 0.2, 250.0
    Ye = 0.08
    temp = np.logspace(np.log10(T0), np.log10(T1), 40)


    gas = EquationOfStateEvaluator([IdealAdiabatic])
    gamma1 = [gas.gamma_effective(D,T,Ye,method=1)-1.4 for T in temp]
    gamma2 = [gas.gamma_effective(D,T,Ye,method=2)-1.4 for T in temp]
    gamma3 = [gas.gamma_effective(D,T,Ye,method=3)-1.4 for T in temp]

    plt.semilogx(temp, gamma1, '-o', lw=1.5, label='method 1')
    plt.semilogx(temp, gamma2, '-s', lw=1.5, label='method 2')
    plt.semilogx(temp, gamma3, '-x', lw=1.5, label='method 3')

    plt.title(r"Error in $\Gamma$ for an adiabatic equation of state")

    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"$\Delta \Gamma_{\rm{eff}}$", fontsize=16)
    plt.legend(loc='upper left')
    plt.show()



def TestGamma2():
    D = 1e13
    T0, T1 = 0.2, 250.0
    Ye = 0.08
    temp = np.logspace(np.log10(T0), np.log10(T1), 4)

    gas = eos.ShenNucleons()
    gamma = [[gas.gamma_effective(D*pq.g/pq.cm**3, T*pq.MeV, Ye, method=m) for
              T in temp] for m in [1,2,3]]

    print '1', gamma[0]
    print '2', gamma[1]
    print '3', gamma[2]

    plt.semilogx(temp, gamma[0], '-o', lw=1.5, label='method 1')
    plt.semilogx(temp, gamma[1], '-s', lw=1.5, label='method 2')
    plt.semilogx(temp, gamma[2], '-x', lw=1.5, label='method 3')

    plt.title(r"Error in $\Gamma$ for an adiabatic equation of state")

    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"$\Delta \Gamma_{\rm{eff}}$", fontsize=16)
    plt.legend(loc='upper left')
    plt.show()



def TestGamma3():
    D0, D1 = 1e6, 1e15
    T0, T1 = 0.2, 250.0
    Ye = 0.08
    temp = np.logspace(np.log10(T0), np.log10(T1), 32)
    dens = np.logspace(np.log10(D0), np.log10(D1), 32)

    gas = eos.ShenNucleons()
    gamma = np.array([[gas.gamma_effective(D*pq.g/pq.cm**3, T*pq.MeV, Ye, method=1) for
                       T in temp] for D in dens])

    extent = [np.log10(D0), np.log10(D1), np.log10(T0), np.log10(T1)]
    aspect = (extent[1] - extent[0]) / (extent[3] - extent[2])

    majorFormatter = ticker.FuncFormatter(
        lambda x, pos: r"$10^{%d}$" % int(x) if np.floor(x) == x else "")

    def do_image(C, title=None, Yi=10):
        fig = plt.figure()
        plt.imshow(C.T, origin='lower', extent=extent, aspect=aspect,
                   interpolation='nearest')

        plt.colorbar()
        fig.suptitle(title, fontsize=14)
        fig.axes[0].get_xaxis().set_major_formatter(majorFormatter)
        fig.axes[0].get_yaxis().set_major_formatter(majorFormatter)
        plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
        plt.ylabel(r"$T \ \rm{MeV}$", fontsize=16)

    do_image(gamma, title=r"$\Gamma_{\rm{eff}}$ with method 1, $Y_p=0.08$")
    plt.show()


TestGamma3()
