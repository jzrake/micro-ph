#!/usr/bin/env python

import numpy as np
from eospy.physics import *
from eospy import quantities as pq
import math




class EquationOfStateBuilder(object):

    def build_terms(self, args):
        return [IdealAdiabatic(*args)]

    def get_vars(self):
        return { 'n': 0, 'T': 1 }



D = 1.2 * pq.kg / pq.meter**3.0
n = D / (28*pq.constants.proton_mass)
T = 293 * pq.Kelvin

builder = EquationOfStateBuilder()
gas = EquationOfStateEvaluator(builder)

print gas.number_density(n, T).rescale('1/cm^3')
print gas.pressure(n, T).rescale('atm')
print gas.entropy(n, T).rescale('MeV/K')

print gas.derivative('entropy', 'n', n, T).rescale('(MeV/K)/(1/cm^3)')
print gas.gamma_effective1(n, T)
print gas.gamma_effective2(n, T)
print gas.gamma_effective3(n, T)


photon_gas = BlackbodyPhotons(40.0 * pq.MeV)
print photon_gas.pressure()

