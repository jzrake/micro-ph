

LIGHT_SPEED        = 2.997924580e+08 # m/s
ATOMIC_MASS_UNIT   = 1.660538860e-27 # kg
PROTON_MASS        = 1.672621580e-27 # kg
KB                 = 1.380648813e-23 # J/K

GRAM = 1e-03 # kg
CM3  = 1e-06 # m^3
FM3  = 1e-45 # m^3
EV   = 1.602176487e-19 # J
MEV  = 1.602176487e-13 # J
ERG  = 1e-7 # J



class DimensionalQuantity(object):
    
    def __init__(self, U, default_unit=None):
        if type(U) in [tuple, list]:
            self.val = U[0] * self._units[U[1]]
        elif isinstance(U, DimensionalQuantity):
            self.val = U.val # already in SI
        else: # assume float
            self.val = U # assume SI units

        if default_unit is not None:
            assert(default_unit in self._units)
            self._default_unit = default_unit

    def convert_to(self, unit=None):
        if unit is None: unit = self._default_unit
        return self.val / self._units[unit]

    def scale(self, v):
        return type(self)(self.val*v)

    def measured_in(self, unit):
        return type(self)(self.val, default_unit=unit)

    def __str__(self):
        return "%4.3e [%s]" % (self.convert_to(), self._default_unit)



class Mass(DimensionalQuantity):
    _default_unit = 'kg'
    _units = {
        'kg'  : 1.0, # SI are base units
        'g'   : GRAM,
        'MeV' : MEV / (LIGHT_SPEED * LIGHT_SPEED)
        }


class NumberDensity(DimensionalQuantity):

    _default_unit = '1/m^3'
    _units = {
        '1/m^3'   : 1.0, # SI are base units
        '1/cm^3'  : 1.0 / CM3,
        '1/fm^3'  : 1.0 / FM3
        }


class MassDensity(DimensionalQuantity):

    _default_unit = 'kg/m^3'
    _units = {
        'kg/m^3'   : 1.0, # SI are base units
        'g/cm^3'   : GRAM / CM3,
        'amu/cm^3' : ATOMIC_MASS_UNIT / CM3,
        'MeV/fm^3' : (MEV / FM3) / (LIGHT_SPEED * LIGHT_SPEED),
        'g/fm^3'   : GRAM / FM3 }

    def to_number_density(self, m):
        """
        Converts the mass density to a number density, given the particle mass.
        """
        n = self.val / Mass(m).convert_to('kg')
        return NumberDensity(n)


class EnergyDensity(DimensionalQuantity):

    _default_unit = 'J/m^3'
    _units = {
        'J/m^3'     : 1.0, # SI are base units
        'Pa'        : 1.0, # Pascal, alias for 'J/m^3
        'kPa'       : 1e3, # kilo-Pascal
        'atm'       : 1.0132e5, # atmospheres
        'erg/cm^3'  : ERG / CM3,
        'MeV/fm^3'  : MEV / FM3 }


class Temperature(DimensionalQuantity):

    _default_unit = 'K'
    _units = {
        'K'        : 1.0, # Kelvins are base units
        'J'        : 1.0 / KB, 
        'eV'       : EV / KB,
        'MeV'      : MEV / KB }


class Entropy(DimensionalQuantity):

    _default_unit = 'J/K'
    _units = {
        'J/K'      : 1.0, # Joule/Kelvin is base units
        'erg/K'    : ERG,
        'eV/K'     : EV,
        'MeV/K'    : MEV }


BoltzmannConstant = Entropy(KB)
RoomTemperature = Temperature(293.0)
AtmosphericPressure = EnergyDensity([1.0, 'atm'])
AtmosphericNumberDensity = NumberDensity([0.02504e27, '1/m^3'])
ProtonMass = Mass(PROTON_MASS)
AtomicMassUnit = Mass(ATOMIC_MASS_UNIT)
