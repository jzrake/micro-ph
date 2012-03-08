#!/usr/bin/env python


if __name__ == "__main__":
    from _eospy import set_verbose, pressure

    set_verbose(0)
    print [pressure(1e13, T, 0.08, photons=True) for T in [10,20,40,80]]


