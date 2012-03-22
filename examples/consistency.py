#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np


def PlotConsistency(eos, T0, T1, D0, D1, Ye=0.08, Kelvin=True, N=60):
    """
    Checks the thermodynamic consistency relation:

    p = rho^2 d(u)/d(rho) + T dp/dT
    """

    def consistency(D, T, Y):
        p    = eos.pressure(D, T, Y)
        dpdT = eos.pressure(D, T, Y, derivative='T')
        dedD = eos.specific_internal_energy(D, T, Y, derivative='D')

        return abs(1.0 - (D**2 * dedD + T*dpdT)/p)

    temp = np.logspace(T0, T1, 3)
    dens = np.logspace(D0, D1, N)

    for T in temp:
        kT = T * BOLTZMANN_CONSTANT if Kelvin else T
        c = np.array([consistency(D, kT, Ye) for D in dens])
        lab = r"$T=10^{%d}\rm{K}$" % np.log10(T) if Kelvin \
            else r"$T=%2.1f\rm{MeV}$" % T
        plt.loglog(dens, c, '-o', label=lab)

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
        plt.imshow(np.log10(abs(1.0 - P/p)[2:-2,2:-2,Yi]).T, origin='lower',
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


ShenConsistency()


#eos = EquationOfStateEvaluator([FermiDiracElectrons, FermiDiracPositrons])
#eos = EquationOfStateEvaluator([BlackbodyPhotons])
#PlotConsistency(eos, 5, 11, -3, 14)

#eos = EquationOfStateEvaluator([NucleonsShenEos3])
#eos.set_numerical_derivative_step(1e-4)
#PlotConsistency(eos, np.log10(0.5), np.log10(260.0), 6, 15, Kelvin=False, N=5)

