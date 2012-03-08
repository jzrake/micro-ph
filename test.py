#!/usr/bin/env python


if __name__ == "__main__":
    from _eospy import set_verbose, pressure
    import numpy as np
    from matplotlib import pyplot as plt

    set_verbose(2)
    pressure(1e10, 0.01, 0.08, electrons=True)
    #pressure(1e10, 0.01, 0.08, positrons=True)

    Ts = np.logspace(-2.0, 2.0, 1000)

    phot = [pressure(1e13, T, 0.08, photons=True) for T in Ts]
    elec = [pressure(1e13, T, 0.08, electrons=True) for T in Ts]

    plt.loglog(Ts, phot)
    plt.show()
