

from matplotlib import pyplot as plt
import numpy as np

import shen
import physics


__all__ = ["test_compare_pressure",
           "test_load_shen",
           "test_sample"]


def test_pressure_vs_density(Ye=0.08):
    """
    Plots the pressure for various components of the EOS.
    """
    temp = np.logspace(5.0, 10.0, 6) # in Kelvins
    dens = np.logspace(-3, 14, 100)

    for T in temp:
        kT = T * physics.BOLTZMANN_CONSTANT
        npumu = [physics.eos(D, kT, Ye, "electrons") for D in dens]
        p  = [n[1] for n in npumu]
        plt.loglog(dens, mu, '-', lw=1.5)
        plt.loglog(dens, p, '-', lw=1.5, label=r"$T=10^{%d}\rm{K}$" % np.log10(T))

    cold = [physics.eos(D, kT, Ye, "cold_electrons")[1] for D in dens]

    plt.loglog(dens, cold, c='k', ls=':', label=r"cold $e_-$")
    plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='lower right')
    plt.show()



def test_compare_pressure(D=1e13, T0=5.0, T1=80.0, Ye=0.08, terms="all"):
    """
    Plots the pressure for various components of the EOS.
    """
    temp = np.logspace(np.log10(T0), np.log10(T1), 60)

    for comp, tex, ls in [("nucleons", r"$n$ (shen)", '-'),
                          ("positrons", r"$e_+$", '--'),
                          ("electrons", r"$e_-$", '-.'),
                          ("photons", r"$\gamma$", '-x'),
                          ("cold_electrons", r"$e_-$, cold", ':')]:

        if comp not in terms and terms is not "all": continue

        print "Working out component '%s'" % comp
        p = [physics.eos(D, T, Ye, comp)[1] for T in temp]
        plt.loglog(temp, p, ls, label=tex)

    plt.xlabel(r"$k_B T$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='lower right')
    plt.show()


def test_sample():
    table = shen.read_hdf5("data/eos3.h5", cols=['log10_rhoB', 'logT', 'p', 'Yp'])
    print shen.sample(table, 'p', 10**5.1, 0.1, 0.01)
    print shen.sample(table, 'p', 1e13, 40.0, 0.08)
