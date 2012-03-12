


from matplotlib import pyplot as plt
import numpy as np


__all__ = ["test_fdfunc1"]


def fdfunc1(x, n, eta, beta):
    return np.power(x,0.5*n) * np.sqrt(1 + x/(2*beta)) / (np.exp(x+beta-eta) + 1)

def fdfunc2(x, n, eta):
    return x**n / (np.exp(x-eta) + 1)


def fdfunc(n, eta, beta):
    """
    http://adsabs.harvard.edu/abs/1971A%26A....13..209B
    """
    pass



def test_fdfunc1():

    for e in [1.0, 2.0, 10.0, 100.0]:
        x = np.linspace(0,5*e,1000)
        y = [fdfunc2(xi, 0.0, e) for xi in x]
        plt.plot(x/e, y, lw=2.0, label=r"$\eta=%2.1f$"%e)

    plt.title("Non-relativistic Fermi-Dirac functions")
    plt.xlabel(r"$\epsilon/\mu$", fontsize=18)
    plt.ylabel(r"$\frac{1}{e^{(\epsilon-\mu)/kT}+1}$", fontsize=24)
    plt.ylim(-0.1, 1.1)
    plt.legend(loc='best')
    plt.show()
