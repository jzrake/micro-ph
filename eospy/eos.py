
import quantities as pq
import physics
import cache

class AdiabaticGas(physics.EquationOfStateEvaluator):

    _density_var = 'n'
    _vars = { 'n': 0, 'T': 1 }

    @cache.memoized()
    def build_terms(self, args):
        gas = physics.IdealAdiabatic(*args)
        gas.particle_mass = 28*pq.constants.proton_mass
        return [gas]


class AdiabaticGasWithDensity(physics.EquationOfStateEvaluator):

    _density_var = 'D'
    _vars = { 'D': 0, 'T': 1, 'Y': 2 }
    _mp = 28*pq.constants.proton_mass

    @cache.memoized()
    def build_terms(self, args):
        """
        Builds the EOS of an ideal gas with adiabatic EOS.
        """
        D, T, Y = args
        n = D / self._mp
        gas = physics.IdealAdiabatic(n, T)
        gas.particle_mass = self._mp
        return [gas]


class ElectronPositronGas(physics.EquationOfStateEvaluator):

    _density_var = 'D'
    _vars = { 'D': 0, 'T': 1, 'Y': 2 }

    @cache.memoized()
    def build_terms(self, args):
        D, T, Yp = args
        np = Yp * D / pq.constants.electron_mass
        ele = physics.FermiDiracElectrons(np, T)
        pos = physics.FermiDiracPositrons(np, T)
        return [ele, pos]
