

LIGHT_SPEED        = 2.997924580e+8 # m/s
ATOMIC_MASS_UNIT   = 1.660538860e-27 # kg
KB                 = 1.380648813e-23 # J/K

GRAM = 1e-03 # kg
CM3  = 1e-06 # m^3
FM3  = 1e-45 # m^3
EV   = 1.602176487e-19 # J
MEV  = 1.602176487e-13 # J
ERG  = 1e-7 # J


class DimensionalQuantity(object):
    
    def __init__(self, val, unit):
        self.val = val * self._units[unit]

    def convert_to(self, unit):
        return self.val / self._units[unit]



class NumberDensity(DimensionalQuantity):

    _units = {
        '1/m^3'   : 1.0, # SI are base units
        '1/cm^3'  : 1.0 / CM3,
        '1/fm^3'  : 1.0 / FM3
        }


class MassDensity(DimensionalQuantity):

    _units = {
        'kg/m^3'   : 1.0, # SI are base units
        'g/cm^3'   : GRAM / CM3,
        'amu/cm^3' : ATOMIC_MASS_UNIT / CM3,
        'MeV/fm^3' : (MEV / FM3) / (LIGHT_SPEED * LIGHT_SPEED),
        'g/fm^3'   : GRAM / FM3 }


class EnergyDensity(DimensionalQuantity):

    _units = {
        'J/m^3'     : 1.0, # SI are base units
        'Pa'        : 1.0, # Pascal, alias for 'J/m^3'
        'erg/cm^3'  : ERG / CM3,
        'MeV/fm^3'  : MEV / FM3 }


class Temperature(DimensionalQuantity):

    _units = {
        'K'        : 1.0, # Kelvins are base units
        'J'        : 1.0 / KB, 
        'eV'       : EV / KB,
        'MeV'      : MEV / KB }


class Entropy(DimensionalQuantity):

    _units = {
        'J/K'      : 1.0, # Joule/Kelvin is base units
        'erg/K'    : ERG,
        'eV/K'     : EV,
        'MeV/K'    : MEV }

