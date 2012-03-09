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


def test_compare_pressure(D=1e13):
    temp = np.logspace(0.5, 1.8, 20)

    pos = [eos_eval(D, T, 0.08, ["positrons"])[1] for T in temp]
    plt.loglog(temp, pos, '-x', label=r"$e_+$")

    pho = [eos_eval(D, T, 0.08, ["photons"])[1] for T in temp]
    plt.loglog(temp, pho, '--', label=r"$\gamma$")

    ele = [eos_eval(D, T, 0.08, ["electrons"])[1] for T in temp]
    plt.loglog(temp, ele, '-o', label=r"$e_-$")

    cle = [eos_eval(D, T, 0.08, ["cold_electrons"])[1] for T in temp]
    plt.loglog(temp, cle, '-x', label=r"$e_-$, cold")

    plt.xlabel(r"$k_B T$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='lower right')
    plt.show()


def test_load_shen():
    from eospy import shen
    table = shen.load_eos3("data/eos3_mini.tab")

    shen.write_hdf5(table, "data/shen.hdf5")
    return

    dens = 10**table[1,:,0,var_index['log10_rhoB']]
    pres =     table[1,:,0,var_index['p']]

    plt.xlabel(r"$k_B T$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.loglog(dens, pres, '-.')
    plt.show()


if __name__ == "__main__":
    #test_compare_pressure()
    test_load_shen()
