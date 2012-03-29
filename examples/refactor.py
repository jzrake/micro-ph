#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np
from eospy.physics import *
from eospy.units import *


class MyEos(object):

    _output_units = { EnergyDensity: 'atm',
                      NumberDensity: '1/cm^3',
                      MassDensity: 'kg/m^3',
                      Temperature: 'K' }

    _variable_dict = {'n': 0, 'T': 1 }
    _input_units = [NumberDensity, Temperature]
    _num_deriv_dx = 1e-8

    def _build_terms(self, A):
        n = A[0]
        T = A[1]
        return [IdealAdiabatic(n, T, gamma=1.4)]

    def _sample(self, A, key, derivative):
        """
        In pseudo code, for example if key is 'pressure':

        return sum([term.pressure(A) for term in terms])

        A is an array of arguments for the EOS, each of which must be a
        dimensional quantity. Return value is in physics base (SI) units.
        """
        if not derivative:
            return sum([getattr(term, key)(asobj=False)
                        for term in self._build_terms(A)])

        elif len(derivative) == 1:
            n = self._variable_dict[derivative]
            dx = self._num_deriv_dx
            X0, X1 = list(A), list(A)

            X0[n] = X0[n].scale(1.0 - dx)
            X1[n] = X1[n].scale(1.0 + dx)

            f = lambda a: self._sample(a, key, None)
            return (f(X1) - f(X0)) / (X1[n].val - X0[n].val)

        elif len(derivative) == 2:
            raise NotImplementedError(
                "Second derivative calculation not written.")
        else:
            raise ValueError(
                "Derivative string must be 0, 1, or 2 characters")


    def _call_sample(self, A, key, unitClass, derivative):
        B = [U(a) for U,a in zip(self._input_units, A)]
        return unitClass(self._sample(B, key, derivative),
                         default_unit=self._output_units[unitClass])

    def pressure(self, n, T, derivative=None):
        return self._call_sample([n, T], 'pressure',
                                 EnergyDensity, derivative)





rho = MassDensity([1.2, 'kg/m^3'])
n = rho.to_number_density(ProtonMass.scale(28.0))
T = RoomTemperature


gas = IdealAdiabatic(n, T, gamma=1.4)

print gas.pressure()
print gas.entropy()

print n, n.measured_in('1/cm^3')
print gas.specific_internal_energy()
print 5./2. * T.convert_to('J')

eos = MyEos()


n = MassDensity([1.2, 'kg/m^3']).to_number_density(ProtonMass.scale(28.0))
T = [293.0, 'K']

print eos.pressure(n, T)
print eos.pressure(n, T, derivative='T')

#print eos.pressure(1e13, 40.0, 0.08)
#print eos.
