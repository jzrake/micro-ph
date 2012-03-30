#!/usr/bin/env python

import numpy as np
from eospy.physics import *
from eospy import quantities as pq
import math

class MyEos(object):

    _units = {
        'number_density': 1/pq.m**3,
        'pressure': pq.pascal,
        'internal_energy': pq.J/pq.m**3,
        'specific_internal_energy': pq.J,
        'entropy': pq.J/pq.K
        }
    _num_deriv_dx = 1e-8
    _vars = { 'n':0, 'T': 1 }

    def _build_terms(self, args):
        return [IdealAdiabatic(*args)]

    def _call_attr(self, args, attr):
        return sum([getattr(term, attr)() for term in
                    self._build_terms(args)], 0.0 * self._units[attr])

    def number_density(self, *args):
        return self._call_attr(args, 'number_density')

    def pressure(self, *args):
        return self._call_attr(args, 'pressure')

    def internal_energy(self, *args):
        return self._call_attr(args, 'internal_energy')

    def specific_internal_energy(self, *args):
        return self._call_attr(args, 'specific_internal_energy')

    def entropy(self, *args):
        return self._call_attr(args, 'entropy')

    def derivative(self, attr, var, *args):
        dx = self._num_deriv_dx
        X0, X1 = list(args), list(args)

        n = self._vars[var]

        X1[n] = X1[n]*(1 + dx)
        X0[n] = X0[n]*(1 - dx)

        f = lambda x: self._call_attr(x, attr)
        return (f(X1) - f(X0)) / (X1[n] - X0[n])






D = 1.2 * pq.kg / pq.meter**3.0
n = D / (28*pq.constants.proton_mass)
T = 293 * pq.Kelvin

gas = MyEos()

print gas.number_density(n, T).rescale('1/cm^3')
print gas.pressure(n, T).rescale('atm')
print gas.entropy(n, T).rescale('MeV/K')

print gas.derivative('entropy', 'D', n, T).rescale('(MeV/K)/(1/cm^3)')
