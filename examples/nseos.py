#!/usr/bin/env python

import numpy as np
import h5py
import quantities as pq
from matplotlib import pyplot as plt
from eospy.physics import *
from eospy.cache import memoized
from eospy import eos


def do_pressure(D, N=12):
    Yp = 0.08
    T0, T1 = 1.0, 200.0
    temp = np.logspace(np.log10(T0), np.log10(T1), N) * pq.MeV
    gas = eos.NeutronStarEos()

    result = [[t.pressure() for t in gas.build_terms([D, T, Yp])] for T in temp]
    tex = gas.tex_description()
    ls = ['--', '-.', '-.', '-', ':', '--']
    lw = [1, 2, 2, 2, 2, 2]

    plt.figure()
    for i in range(6):
        plt.loglog(temp, [r[i] for r in result], ls[i], lw=lw[i], label=tex[i])

    plt.title(r"Various EOS components for $Y_p=0.08$ $\rho=10^{%d} \ \rm{g/cm^3}$" %
              int(np.log10(D.magnitude)))
    plt.ylim(1e-10, 1e2)
    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"$p(\rho,T) \ \rm{MeV/fm^3}$", fontsize=16)
    plt.legend(loc='lower right')



def do_soundspeed(gas, D, N=12, method=2):
    Yp = 0.08
    T0, T1 = 1.0, 200.0
    temp = np.logspace(np.log10(T0), np.log10(T1), N) * pq.MeV

    result = [(gas.sound_speed(D, T, Yp, method=method) / pq.c).rescale('')
              for T in temp]
    print result

    plt.semilogx(temp, result, '-o', mfc='none', lw=1.5,
                 label=r"$\rho=10^{%d} \ \rm{g/cm^3}$"%
                 int(np.log10(D.magnitude)))
    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"$c_s/c$", fontsize=16)
    plt.legend(loc='lower right')



def do_gamma(gas, D, N=12, method=2):
    Yp = 0.08
    T0, T1 = 1.0, 200.0
    temp = np.logspace(np.log10(T0), np.log10(T1), N) * pq.MeV

    result = [gas.gamma_effective(D, T, Yp, method=method) for T in temp]
    print result

    plt.semilogx(temp, result, '-o', mfc='none', lw=1.5,
                 label=r"$\rho=10^{%d} \ \rm{g/cm^3}$"%
                 int(np.log10(D.magnitude)))
    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"$\Gamma_{\rm{eff}}$", fontsize=16)
    plt.legend(loc='upper right')



def Pressure():
    D = 1e13 * pq.g/pq.cm**3
    N = 96
    do_pressure(1e-1*D, N)
    do_pressure(1e+0*D, N)
    do_pressure(1e+1*D, N)
    plt.show()



def SoundSpeed():
    gas = eos.NeutronStarEos()
    memoized.load_cache('eos.cache')

    D = 1e13 * pq.g/pq.cm**3
    N = 96

    for m in [2,3]:
        do_soundspeed(gas, 1e-1*D, N, method=m)
        do_soundspeed(gas, 1e+0*D, N, method=m)
        do_soundspeed(gas, 1e+1*D, N, method=m)
        plt.savefig("nseos_cs_method%d.png" % m)
        plt.clf()

        do_gamma(gas, 1e-1*D, N, method=m)
        do_gamma(gas, 1e+0*D, N, method=m)
        do_gamma(gas, 1e+1*D, N, method=m)
        plt.savefig("nseos_gamma_method%d.png" % m)
        plt.clf()

    memoized.save_cache('eos.cache')



def ChemicalPotential():
    Yp = 0.08
    D = 1e13 * pq.g/pq.cm**3
    T0, T1 = 12.0, 200.0
    temp = np.logspace(np.log10(T0), np.log10(T1), 64) * pq.MeV

    gas = eos.NeutronStarEos()

    nuc = [gas.get_term('nucleons', [D, T, Yp]) for T in temp]
    ele = [gas.get_term('electrons', [D, T, Yp]) for T in temp]
    pos = [gas.get_term('positrons', [D, T, Yp]) for T in temp]
    nut = [gas.get_term('neutrinos', [D, T, Yp]) for T in temp]
    aut = [gas.get_term('antineutrinos', [D, T, Yp]) for T in temp]

    mu_protons = [x.chemical_potential('protons') for x in nuc]
    mu_neutrons = [x.chemical_potential('neutrons') for x in nuc]
    mu_electrons = [x.chemical_potential() for x in ele]
    mu_positrons = [x.chemical_potential() for x in pos]
    mu_nu = [x.chemical_potential() for x in nut]
    mu_nu_bar = [x.chemical_potential() for x in aut]

    plot = plt.plot

    plot(temp, mu_protons, lw=2, ls=':', label=r"protons")
    plot(temp, mu_neutrons, lw=2, ls='--', label=r"neutrons")
    plot(temp, mu_electrons, lw=2, ls=':', label=r"electrons")
    plot(temp, mu_positrons, lw=2, ls='--', label=r"positrons")
    plot(temp, mu_nu, lw=2, ls=':', label=r"electron neutrinos")
    plot(temp, mu_nu_bar, lw=2, ls='--', label=r"anti-electron neutrinos")

    plt.xlabel(r"$k_B T \ \rm{MeV}$", fontsize=16)
    plt.ylabel(r"chemical potential $\rm{MeV}$", fontsize=16)
    plt.legend(loc='upper right')
    plt.show()


def WriteTable(ND=16, NT=16):
    Yp = 0.08
    D0, D1 = 0.5e12, 0.5e15
    T0, T1 = 0.1, 100.0
    dens = np.linspace(D0, D1, ND) * pq.g/pq.cm**3
    temp = np.linspace(T0, T1, NT) * pq.MeV

    gas = eos.AdiabaticGasWithDensity()
    gas._mp = pq.constants.proton_mass
    #gas = eos.NeutronStarEos()
    memoized.load_cache('eos.cache')

    def eval_func(func, units):
        return np.array([[getattr(gas, func)(D, T, Yp).rescale(units).
                          magnitude for T in temp] for D in dens])

    table = h5py.File("adeos.h5", 'w')

    table["density"    ] = np.array([[D.magnitude for T in temp] for D in dens])
    table["temperature"] = np.array([[T.magnitude for T in temp] for D in dens])

    table["pressure"] = eval_func("pressure", "MeV/fm^3")
    table["internal_energy"] = eval_func("internal_energy", "MeV/fm^3")
    table["sound_speed"] = eval_func("sound_speed", "c")

    print np.array([[D.magnitude for T in temp] for D in dens])
    print np.array([[T.magnitude for T in temp] for D in dens])

    memoized.save_cache('eos.cache')


#ChemicalPotential()
#Pressure()
#SoundSpeed()
WriteTable(ND=32, NT=32)

