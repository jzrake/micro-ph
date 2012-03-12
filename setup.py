#!/usr/bin/env python

from numpy.distutils.core import setup, Extension
setup(
    ext_modules = [Extension('fdfunc',
                             ['src/fdfunc/cloutman.f90',
                              'src/fdfunc/fermi_dirac_quadrature.f90'])],
)
