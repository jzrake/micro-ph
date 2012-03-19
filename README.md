
# Introduction

Micro-ph is a collection of microphysics code, accessible through a well-written
Python module, for developing astrophysical equations of state. Special
attention is paid to relativistic physics needed at high densities and
temperatures relevant to neutron star material.

The Python code in this project is intended to be used *off-line* to generate a
tabulated EOS which can then be used from inside a C/C++ application. Therefore,
the EOS generating functions don't have to be brutally fast as they're only
executed once.

The project intends to provide a minimal C-library to read, interpolate, invert,
and take derivatives from within a C-application. This functionality will be
available very soon.


# Features

## Physics

* Relativistic and degenerate electron/positron pairs

    * Fermi-Dirac integrals evaluated either with Scipy integration, or through
      a [Fortran module](http://cococubed.asu.edu/code_pages/fermi_dirac.shtml)
      due to Frank Timmes.

    * Asymptotic forms for cold electrons in the 5/3 and 4/3 (high density)
      limit.

* Photons (radiation pressure)
* Nucleons, via the tabulated *Shen* equation of state
* Thermodynamic consistency relations


## Data formats

* Reads and writes ASCII and HDF5
* Contains a parser for the Shen *eos3* table


# Planned features

## Physics

* Neutrinos in beta-equilibrium with pair plasma
* Photo-ionization
* Full EOS of [Timmes](http://cococubed.asu.edu/code_pages/eos.shtml)
* Full EOS of [Lattimer & Swesty](http://www.astro.sunysb.edu/dswesty/lseos.html)

