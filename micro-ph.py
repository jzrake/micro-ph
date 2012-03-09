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
    temp = np.logspace(1.0, 1.8, 4)
    #pos = [eos_eval(D, T, 0.08, ["positrons"])[1] for T in temp]
    ele = [eos_eval(D, T, 0.08, ["electrons"])[1] for T in temp]
    pho = [eos_eval(D, T, 0.08, ["hot_pairs"])[1] for T in temp]
    #plt.loglog(temp, pos, '-x', label=r"$e_+$")
    plt.loglog(temp, ele, '-o', label=r"$e_-$")
    plt.loglog(temp, pho, '--', label=r"$\gamma$")

    plt.xlabel(r"$k_B T$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='lower right')
    plt.show()


if __name__ == "__main__":
    test_compare_pressure()
