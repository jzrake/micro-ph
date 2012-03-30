

import physics
import quantities as pq


class AdiabaticGas(physics.EquationOfStateEvaluator):

    _vars = { 'n': 0, 'T': 1 }

    def build_terms(self, args):
        gas = physics.IdealAdiabatic(*args)
        gas.particle_mass = 28*pq.constants.proton_mass
        return [gas]

    def get_vars(self, key):
        return { 'n': 0, 'T': 1 }[key]


class ElectronPositronGas(physics.EquationOfStateEvaluator):

    _vars = { 'n': 0, 'T': 1 }

    def build_terms(self, args):
        ele = physics.FermiDiracElectrons(*args)
        pos = physics.FermiDiracPositrons(*args)
        return [ele, pos]
