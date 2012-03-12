


from matplotlib import pyplot as plt
import numpy as np


__all__ = [ ]


def fdfunc1_sr(x, n, eta, beta):
    """
    http://cococubed.asu.edu/code_pages/fermi_dirac.shtml
    """
    return np.power(x,n) * np.sqrt(1 + 0.5*x*beta) / (np.exp(x-eta) + 1)


def fdfunc1_nr(x, n, eta):
    return np.power(x,n) / (np.exp(x-eta) + 1)


def fdfunc(n, eta, beta):
    """
    http://adsabs.harvard.edu/abs/1971A%26A....13..209B
    """
    pass



def test_fdfunc1_nr():

    for e in [1.0, 2.0, 10.0, 100.0]:
        x = np.linspace(0,5*e,1000)
        y = [fdfunc1_nr(xi, 0.0, e) for xi in x]
        plt.plot(x/e, y, lw=2.0, label=r"$\eta=%2.1f$"%e)

    plt.title("Non-relativistic Fermi-Dirac functions")
    plt.xlabel(r"$\epsilon/\mu$", fontsize=18)
    plt.ylabel(r"$\frac{1}{e^{(\epsilon-\mu)/kT}+1}$", fontsize=24)
    plt.ylim(-0.1, 1.1)
    plt.legend(loc='best')
    plt.show()


def test_fdfunc1_sr():

    for c,e in zip("rgby", [1.0, 2.0, 10.0, 100.0]):
        x = np.linspace(0,5*e,1000)
        y0 = [fdfunc1_nr(xi, 0.0, e) for xi in x]
        y1 = [fdfunc1_sr(xi, 0.0, e, 10.0) for xi in x]
        plt.plot(x/e, y0, lw=1.2, ls='-.', c=c)
        plt.plot(x/e, y1, lw=1.6, c=c, label=r"$\eta=%2.1f$"%e)

    plt.title("Relativistic Fermi-Dirac functions")
    plt.xlabel(r"$\epsilon/\mu$", fontsize=18)
    plt.ylabel(r"$\frac{\sqrt{1+x\beta/2}}{e^{(\epsilon-\mu)/kT}+1}$", fontsize=24)
    plt.ylim(-0.2, 2.0)
    plt.legend(loc='best')
    plt.show()
