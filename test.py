#!/usr/bin/env python


if __name__ == "__main__":
    from _eospy import set_verbose, pressure
    import numpy as np
    from matplotlib import pyplot as plt

    set_verbose(2)
    for T in np.logspace(0.5, 2.0, 10):
        print T, pressure(1e13, T, 0.08, electrons=True), pressure(1e10, T, 0.08, positrons=True)

    exit()


    Ts = np.logspace(-2.0, 2.0, 1000)

    phot = [pressure(1e13, T, 0.08, photons=True) for T in Ts]    

    plt.loglog(Ts, phot)
    plt.show()
