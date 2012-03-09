#!/usr/bin/env python

from eospy import *


def test_large_beta():
    from matplotlib import pyplot as plt
    import numpy as np

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



if __name__ == "__main__":
    #print evaluate_term("number_density", +1, 1, 1)
    #print evaluate_term("number_density", -1, 100.0, 1000.0, massless=False)

    #print eos_eval(1e13, 0.1, 0.08, ["photons"])
    print eos_eval(1e13, 0.5, 0.08, ["electrons"])
    #print eos_eval(1e13, 0.1, 0.08, ["positrons"])
