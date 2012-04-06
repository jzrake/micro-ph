#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
import quantities as pq
from eospy.physics import *
from eospy.eos import *


def PlotConsistency(eos, T0, T1, D0, D1, Ye=0.08, kelvin=True, N=60, title=None):
    """
    Checks the thermodynamic consistency relation:

    p = rho^2 d(u)/d(rho) + T dp/dT
    """

    def consistency(D, T, Y):
        print "evaluating with", D, T, Y
        p    = eos.pressure(D, T, Y)
        dpdT = eos.derivative('pressure', 'T', D, T, Y)
        dedD = eos.derivative('specific_internal_energy', 'D', D, T, Y)
        return abs(1.0 - (D**2 * dedD + T*dpdT)/p)

    temp = np.logspace(T0, T1, 3)
    dens = np.logspace(D0, D1, N)

    for T in temp:
        c = np.array([consistency(D*pq.g/pq.cm**3, T*pq.MeV, Ye) for D in dens])
        lab = pq.MeV.dimensionality.latex
        lab = r"$10^{%d}$" % np.log10(T) + " " + lab
        plt.loglog(dens, c, '-o', label=lab)

    if title: plt.title(title)
    plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
    plt.ylabel(r"$1 - \frac{\rho^2}{p} \frac{\partial E}{\partial \rho} - " +
               r"\frac{T}{p} \frac{\partial p}{\partial T}$", fontsize=18)
    plt.legend(loc='best')
    plt.show()



def ShenConsistency():
    """
    Uses different derivative estimates to check the thermodynamic consistency
    relation

    p = rho^2 d(u)/d(rho) + T dp/dT.
    
    """
    cols = ['log10_rhoB', 'p', 'Eint', 'logT', 'nB', 'Yp', 'S', 'F']
    table = shen.read_hdf5("data/eos3.h5", cols=cols)

    logT = table['logT']
    D = 10**table['log10_rhoB']
    T = 10**table['logT']
    Y = table['Yp']
    F = table['F'] + 938.0
    E = table['Eint'] + 931.494
    S = table['S']
    Z = abs(1.0 - (E - T*S)/F)
    print "average error in free energy definition:", Z.mean()

    p = table['p']
    nB = table['nB']

    print "real pressure:                    ", p[100,80,10]

    # 5-point stencil
    # --------------------------------------------------------------------------
    dEdnB   = shen.derivative5(E, nB, 0)
    dpdT    = shen.derivative5(p, T, 1)
    dpdlogT = shen.derivative5(p, logT, 1)

    Plin5 = nB**2 * dEdnB + T*dpdT
    Plog5 = nB**2 * dEdnB + dpdlogT / np.log(10)

    print "5-point stencil, lin spacing in T:", Plin5[100,80,10]
    print "5-point stencil, log spacing in T:", Plog5[100,80,10]

    # 2-point stencil
    # --------------------------------------------------------------------------
    dEdnB   = shen.derivative2(E, nB, 0)
    dpdT    = shen.derivative2(p, T, 1)
    dpdlogT = shen.derivative2(p, logT, 1)

    Plin2 = nB**2 * dEdnB + T*dpdT
    Plog2 = nB**2 * dEdnB + dpdlogT / np.log(10)

    print "2-point stencil, lin spacing in T:", Plin2[100,80,10]
    print "2-point stencil, log spacing in T:", Plog2[100,80,10]

    def do_image(P, title=None, Yi=10):
        plt.figure()
        plt.imshow(np.log10(abs(1.0 - P/p))[2:-2,2:-2,Yi].T, origin='lower',
                   interpolation='nearest')
        plt.title(title + r" $Y_p=%3.2f$" % Y[0,0,Yi])
        plt.colorbar()
        plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
        plt.ylabel(r"$T \ \rm{MeV}$", fontsize=16)

    do_image(Plin2, title='Linear 2-point stencil', Yi=7)
    do_image(Plog2, title='Log 2-point stencil', Yi=7)

    do_image(Plin5, title='Linear 5-point stencil', Yi=7)
    do_image(Plog5, title='Log 5-point stencil', Yi=7)
    plt.show()



def ShenPressure():
    """
    Uses the 5-point derivative to check the definition of the pressure in terms
    of the Helmholtz free energy.
    """
    cols = cols=['log10_rhoB', 'p', 'F', 'nB', 'Eint', 'logT', 'S']
    table = shen.read_hdf5("data/eos3.h5", cols=cols)

    extent = [table['log10_rhoB'][2,0,0], table['log10_rhoB'][-2,0,0],
              table['logT'][0,2,0], table['logT'][0,-2,0]]
    aspect = (extent[1] - extent[0]) / (extent[3] - extent[2])


    majorFormatter = ticker.FuncFormatter(
        lambda x, pos: r"$10^{%d}$" % int(x) if np.floor(x) == x else "")

    def do_image(C, title=None, Yi=10):
        fig = plt.figure()

        plt.imshow(C[2:-2,2:-2,Yi].T, origin='lower', interpolation='nearest',
                   extent=extent, aspect=aspect)

        plt.colorbar()
        fig.suptitle(title, fontsize=14)
        fig.axes[0].get_xaxis().set_major_formatter(majorFormatter)
        fig.axes[0].get_yaxis().set_major_formatter(majorFormatter)
        plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
        plt.ylabel(r"$T \ \rm{MeV}$", fontsize=16)

    tex = r"$(p - n_B^2 \frac{\partial F}{\partial n_B})/p$"

    dFdnB = shen.derivative2(table['F'], table['nB'], 0)
    p = dFdnB * table['nB']**2
    do_image(np.log10(abs(1.0 - p/table['p'])), title=r"2-point stencil " + tex)

    dFdnB = shen.derivative5(table['F'], table['nB'], 0)
    p = dFdnB * table['nB']**2
    do_image(np.log10(abs(1.0 - p/table['p'])), title=r"5-point stencil " + tex)

    plt.show()


def ShenSoundSpeed(image=True):
    """
    Checks consistency of sound speeds using different formulas and difference
    stencils.
    """
    cols = cols=['log10_rhoB', 'p', 'F', 'nB', 'Eint', 'logT', 'S']
    table = shen.read_hdf5("data/eos3.h5", cols=cols)

    extent = [table['log10_rhoB'][2,0,0], table['log10_rhoB'][-2,0,0],
              table['logT'][0,2,0], table['logT'][0,-2,0]]
    aspect = (extent[1] - extent[0]) / (extent[3] - extent[2])


    majorFormatter = ticker.FuncFormatter(
        lambda x, pos: r"$10^{%d}$" % int(x) if np.floor(x) == x else "")

    Di = 79 # 10^13 g/cm^3
    Yi = 7  # 0.08

    def do_image(C, title=None, vmin=None, vmax=None):
        fig = plt.figure()

        plt.imshow(C[2:-2,2:-2,Yi].T, origin='lower', interpolation='nearest',
                   extent=extent, aspect=aspect, vmin=vmin, vmax=vmax)

        plt.colorbar()
        fig.suptitle(title, fontsize=14)
        fig.axes[0].get_xaxis().set_major_formatter(majorFormatter)
        fig.axes[0].get_yaxis().set_major_formatter(majorFormatter)
        plt.xlabel(r"$\rho \ \rm{g/cm^3}$", fontsize=16)
        plt.ylabel(r"$T \ \rm{MeV}$", fontsize=16)

    def do_csplot(x, y, *args, **kwargs):
        plt.loglog(x[Di,2:-2,Yi], y[Di,2:-2,Yi], lw=1.5, *args, **kwargs)
        plt.ylabel(r"$c_s/c$", fontsize=16)
        plt.xlabel(r"$T \ \rm{MeV}$", fontsize=16)

    def do_gmplot(x, y, *args, **kwargs):
        plt.semilogx(x[Di,2:-2,Yi], y[Di,2:-2,Yi], lw=1.5, *args, **kwargs)
        plt.ylabel(r"$\Gamma_{\rm{eff}}$", fontsize=16)
        plt.xlabel(r"$T \ \rm{MeV}$", fontsize=16)

    tex = r"$(p - n_B^2 \frac{\partial F}{\partial n_B})/p$"
    shen.append_sound_speeds(table)

    if image:
        do_image(table['gamma1'], title=r"$\Gamma_{\rm{eff}} \ \rm{method 1}$")
        do_image(table['gamma2'], title=r"$\Gamma_{\rm{eff}} \ \rm{method 2}$")
        do_image(table['gamma3'], title=r"$\Gamma_{\rm{eff}} \ \rm{method 3}$")

        do_image(table['cs1'], title=r"$c_s \ \rm{method 1}$")
        do_image(table['cs2'], title=r"$c_s \ \rm{method 2}$")
        do_image(table['cs3'], title=r"$c_s \ \rm{method 3}$")
    else:
        plt.figure()
        do_csplot(10**table['logT'], table['cs1'], label='method 1', ls='--')
        do_csplot(10**table['logT'], table['cs2'], label='method 2', ls=':')
        do_csplot(10**table['logT'], table['cs3'], label='method 3', ls='-.')
        plt.title(r"Sound speed from different methods at $\rho=10^{%d}\ \rm{g/cm^3}$"
                  % table['log10_rhoB'][Di,0,0], fontsize=14)
        plt.legend(loc='upper left')

        plt.figure()
        do_gmplot(10**table['logT'], table['gamma1'], label='method 1', ls='--')
        do_gmplot(10**table['logT'], table['gamma2'], label='method 2', ls=':')
        do_gmplot(10**table['logT'], table['gamma3'], label='method 3', ls='-.')
        plt.title(r"Effective $\Gamma$ from different methods at $\rho=10^{%d}\ \rm{g/cm^3}$"
                  % table['log10_rhoB'][Di,0,0], fontsize=14)
        plt.legend(loc='upper left')
    plt.show()


#ShenPressure()
#ShenConsistency()
ShenSoundSpeed(image=False)

#gas = ElectronPositronGas()
#gas = AdiabaticGasWithDensity()

#D = 1e13 * pq.g / pq.cm**3
#T = 40.0 * pq.MeV
#print gas.pressure(D, T, 0.08).rescale('MeV/fm^3')

#PlotConsistency(gas, -2, 2, -3, 14, N=4)

#eos = EquationOfStateEvaluator([BlackbodyPhotons])
#PlotConsistency(eos, -2, 2, -3, 14, kelvin=False, N=20, title=r"$e_+/e_-$ pair consistency")

#eos = EquationOfStateEvaluator([NucleonsShenEos3])
#eos.set_numerical_derivative_step(1e-4)
#PlotConsistency(eos, np.log10(0.5), np.log10(260.0), 6, 15, Kelvin=False, N=5)

