#!/usr/bin/env python

if __name__ == "__main__":
    from eospy import *

    #shen.rebuild_hdf5("data/eos3.tab", "data/eos3.h5")
    tests.test_compare_pressure(T0=2.0, T1=340.0)
    #fdfunc.test_fdfunc1_sr()
    #tests.test_pressure_vs_density()
