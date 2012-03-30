#!/usr/bin/env python

import numpy as np
import quantities as pq
from eospy.physics import *
from eospy.eos import *



D = 1.2 * pq.kg / pq.meter**3.0
n = D / (28*pq.constants.proton_mass)
T = 293 * pq.Kelvin

gas = AdiabaticGas()

print gas.number_density(n, T).rescale('1/cm^3')
print gas.pressure(n, T).rescale('atm')
print gas.entropy(n, T).rescale('MeV/K')

print gas.derivative('entropy', 'n', n, T).rescale('(MeV/K)/(1/cm^3)')
print gas.gamma_effective(n, T, method=1)
print gas.gamma_effective(n, T, method=2)
print gas.gamma_effective(n, T, method=3)

print "sound speed:", gas.sound_speed(n, T)

photon_gas = BlackbodyPhotons(40.0 * pq.MeV)
print photon_gas.pressure()

neutrino_gas = NeutrinoComponent(1.0 * pq.MeV, 40.0 * pq.MeV)
print neutrino_gas.pressure()

electron_gas = FermionComponent(1.0 * pq.MeV, 40.0 * pq.MeV)
print electron_gas.pressure()

cold_electrons = ColdElectrons(1.0 / pq.cm**3)
print cold_electrons.pressure().rescale('MeV/cm^3')

dense_electrons = DenseElectrons(1.0 / pq.cm**3)
print dense_electrons.pressure().rescale('MeV/cm^3')

nucleons = NucleonsShenEos3(1e13*pq.g/pq.cm**3, 40.0*pq.MeV, 0.08)
print nucleons.pressure()
print nucleons.number_density()
