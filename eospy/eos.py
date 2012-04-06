
import quantities as pq
import physics
import shen
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
        np = Yp * D / shen.amu
        ele = physics.FermiDiracElectrons(np, T)
        pos = physics.FermiDiracPositrons(np, T)
        return [ele, pos]


class ShenNucleons(physics.EquationOfStateEvaluator):

    _density_var = 'D'
    _vars = { 'D': 0, 'T': 1, 'Y': 2 }

    @cache.memoized()
    def build_terms(self, args):
        D, T, Yp = args
        nuc = physics.NucleonsShenEos3(D, T, Yp)
        return [nuc]


class NeutronStarEos(physics.EquationOfStateEvaluator):

    _density_var = 'D'
    _vars = { 'D': 0, 'T': 1, 'Y': 2 }
    _components = {'photons': 0,
                   'electrons': 1,
                   'positrons': 2,
                   'nucleons': 3,
                   'neutrinos': 4,
                   'antineutrinos': 5}

    @cache.memoized()
    def build_terms(self, args):
        D, T, Yp = args
        np = Yp * D / shen.amu

        photons = physics.BlackbodyPhotons(T)
        electron = physics.FermiDiracElectrons(np, T)
        positron = physics.FermiDiracPositrons(np, T)
        nucleons = physics.NucleonsShenEos3(D, T, Yp)

        mu_n = nucleons.chemical_potential('neutrons')
        mu_p = nucleons.chemical_potential('protons')
        mu_e = electron.chemical_potential()
        mu_t = positron.chemical_potential()

        # these will be negatives of one another
        mu_nu     = mu_e + mu_p - mu_n # electrons neutrinos
        mu_nu_bar = mu_t + mu_n - mu_p # anti-electron neutrinos

        neutrino = physics.NeutrinoComponent(mu_nu, T)
        aeutrino = physics.NeutrinoComponent(mu_nu_bar, T)

        return [photons, electron, positron, nucleons, neutrino, aeutrino]

    def get_term(self, key, args, tex=False):
        terms = self.build_terms(args)
        n = self._components[key]
        if not tex:
            return terms[n]
        else:
            return terms[n], self.tex_description()[n]

    def tex_description(self):
        return [r"$\gamma$", r"$e_-$", r"$e_+$",
                r"$n$ (shen)", r"$\nu_e$", r"$\bar{\nu}_e$"]
