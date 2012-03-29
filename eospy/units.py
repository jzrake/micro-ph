

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
    
    def __init__(self, U):
        if type(U) is float:
            self.val = U # assume SI units
        elif type(U) in [tuple, list]:
            self.val = U[0] * self._units[U[1]]
        elif isinstance(U, DimensionalQuantity):
            self.val = U.val
        else:
            raise ValueError(
                "expected float, list, tuple, or DimensionalQuantity")

    def convert_to(self, unit=None):
        if unit is None: unit = self.default_unit
        return self.val / self._units[unit]

    def __call__(self, unit=None):
        return self.convert_to(unit)

    def __str__(self):
        return str(self.convert_to(self.default_unit)) + " " + self.default_unit



class Mass(DimensionalQuantity):
    default_unit = 'kg'
    _units = {
        'kg'  : 1.0, # SI are base units
        'g'   : GRAM,
        'MeV' : MEV / (LIGHT_SPEED * LIGHT_SPEED)
        }


class NumberDensity(DimensionalQuantity):

    default_unit = '1/m^3'
    _units = {
        '1/m^3'   : 1.0, # SI are base units
        '1/cm^3'  : 1.0 / CM3,
        '1/fm^3'  : 1.0 / FM3
        }


class MassDensity(DimensionalQuantity):

    default_unit = 'kg/m^3'
    _units = {
        'kg/m^3'   : 1.0, # SI are base units
        'g/cm^3'   : GRAM / CM3,
        'amu/cm^3' : ATOMIC_MASS_UNIT / CM3,
        'MeV/fm^3' : (MEV / FM3) / (LIGHT_SPEED * LIGHT_SPEED),
        'g/fm^3'   : GRAM / FM3 }


class EnergyDensity(DimensionalQuantity):

    default_unit = 'J/m^3'
    _units = {
        'J/m^3'     : 1.0, # SI are base units
        'Pa'        : 1.0, # Pascal, alias for 'J/m^3
        'kPa'       : 1e3, # kilo-Pascal
        'atm'       : 1.0132e5, # atmospheres
        'erg/cm^3'  : ERG / CM3,
        'MeV/fm^3'  : MEV / FM3 }


class Temperature(DimensionalQuantity):

    default_unit = 'K'
    _units = {
        'K'        : 1.0, # Kelvins are base units
        'J'        : 1.0 / KB, 
        'eV'       : EV / KB,
        'MeV'      : MEV / KB }


class Entropy(DimensionalQuantity):

    default_unit = 'J/K'
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
