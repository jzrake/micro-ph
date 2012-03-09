#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np
from eospy import *


def test_large_beta():
    eta = np.linspace(-10.0, 10.0, 5)
    beta = np.logspace(-2.0, 1.5, 40)
    T = [convert_beta_to_T(b) for b in beta]

    for e in eta:
        n = [evaluate_term("number_density", -1, e, B) for B in beta]
        plt.loglog(beta, n, label=r"$\eta = %2.1f$" % e)

    plt.xlabel(r"$\beta = m_e c^2 / k_B T$", fontsize=16)
    plt.ylabel(r"$n_+(\eta,\beta)$", fontsize=16)
    plt.legend(loc='best')
    plt.show()


def test_compare_pressure(D=1e13, Ye=0.08, T0=5.0, T1=80.0):
    """
    Plots the pressure for various components of the EOS.
    """
    temp = np.logspace(np.log10(T0), np.log10(T1), 20)

    for comp, tex, ls in [("positrons", r"$e_+$", '--'),
                          ("electrons", r"$e_-$", '-.'),
                          ("photons", r"$\gamma$", '-x'),
                          ("cold_electrons", r"$e_-$, cold", '-o')]:

        p = [eos_eval(D, T, Ye, comp)[1] for T in temp]
        plt.loglog(temp, p, ls, label=tex)

    plt.xlabel(r"$k_B T$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='lower right')
    plt.show()


def test_load_shen():
    from eospy import shen
    """
    table = shen.load_eos3("data/eos3.tab")
    shen.write_hdf5(table, "data/shen.hdf5")
    """
    table = shen.read_hdf5("data/shen.hdf5", cols=['log10_rhoB', 'logT', 'p', 'Yp'])

    iY = 0
    iD = 80

    D = 10**table['log10_rhoB'][0,iD,iY]
    T = 10**table['logT'      ][:,iD,iY]
    Y =     table['Yp'        ][0,iD,iY]
    p =     table['p'         ][:,iD,iY]
    plt.loglog(T, p, '-.', label="nucleons")

    for comp, tex, ls in [("positrons", r"$e_+$", '--'),
                          ("electrons", r"$e_-$", '-.'),
                          ("photons", r"$\gamma$", '-x'),
                          ("cold_electrons", r"$e_-$, cold", '-o')]:

        p = [eos_eval(D, Ti, Y, comp)[1] for Ti in T]
        plt.loglog(T, p, ls, label=tex)

    plt.xlabel(r"$k_B T$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='best')
    plt.show()


def test_sample():
    from eospy import shen
    table = shen.read_hdf5("data/shen.hdf5", cols=['log10_rhoB', 'logT', 'p', 'Yp'])
    print shen.sample(table, 'p', 0.1, 10**5.1, 0.01)
    #print shen.sample(None, 13, 10.0, 0.08)


if __name__ == "__main__":
    #test_compare_pressure(T0=2.0, T1=90)
    #test_load_shen()
    test_sample()
