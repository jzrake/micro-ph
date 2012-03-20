#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np


if __name__ == "__main__":
    """
    Checks the thermodynamic consistency relation:

    p = rho^2 d(u)/d(rho) + T dp/dT
    """

    eos = EquationOfStateEvaluator([FermiDiracElectrons, FermiDiracPositrons])

    D, T, Y = 1e13, 40.0, 0.08

    p = eos.pressure(D, T, Y)
    dpdT = eos.pressure(D, T, Y, derivative='T')
    dedD = eos.specific_internal_energy(D, T, Y, derivative='D')

    print p - (D**2 * dedD + T*dpdT)
