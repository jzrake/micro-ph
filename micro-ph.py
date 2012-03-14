#!/usr/bin/env python

if __name__ == "__main__":
    from eospy import *

    #shen.rebuild_hdf5("data/eos3.tab", "data/eos3.h5")
    #tests.test_compare_pressure(T0=2.0, T1=340.0, Ye=0.08, terms='all')
    #fdfunc.test_fdintegrand(use_timmes=True)
    #fdfunc.test_fdintegral(use_timmes=False)
    #tests.test_pressure_vs_density()
    tests.test_thermodynamic_consistency(1e13, 40.0, 0.08)
    #tests.test_derivatives()
