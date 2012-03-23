#!/usr/bin/env python

from eospy.physics import *
from matplotlib import pyplot as plt
import numpy as np


D = 1e13
T = 40.0

electron = FermiDiracElectrons(D, T, 0.08)
positron = FermiDiracPositrons(D, T, 0.08)
nucleons = NucleonsShenEos3(D, T, 0.08)

mu_n = nucleons.chemical_potential('neutrons')
mu_p = nucleons.chemical_potential('protons')
mu_e = electron.chemical_potential()
mu_t = positron.chemical_potential()

print "chemical potentials:"
print "e-:", mu_e
print "e+:", mu_t
print "p :", mu_p
print "n :", mu_n

mu_nu_p = mu_e + mu_p - mu_n
mu_nu_n = mu_t + mu_n - mu_p

print "nu+:", mu_nu_p
print "nu-:", mu_nu_n

neutrino = NeutrinoComponent(+1, mu_nu_p, T)
aeutrino = NeutrinoComponent(-1, mu_nu_n, T)

print "pressures:"
print "e- :", electron.pressure()
print "e+ :", positron.pressure()
print "p+n:", nucleons.pressure()
print "nu+:", neutrino.pressure()
print "nu-:", aeutrino.pressure()



