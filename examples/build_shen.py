#!/usr/bin/env python

from eospy import shen

shen.rebuild_hdf5("data/eos3.tab", "data/eos3.h5", cols="all")

