


from matplotlib import pyplot as plt
import numpy as np
import timmes.fdfunc


__all__ = [ ]


def fdfunc1(x, n, eta, beta, use_timmes=False):
    """
    http://cococubed.asu.edu/code_pages/fermi_dirac.shtml
    http://iopscience.iop.org/0067-0049/125/1/277
    """
    if use_timmes:
        return timmes.fdfunc.fdfunc1(x, [n, eta, beta])[0]
    else:
        return np.power(x,n) * np.sqrt(1 + 0.5*x*beta) / (np.exp(x-eta) + 1)


def test_fdintegrand(use_timmes=False):

    for c,e in zip("rgby", [1.0, 2.0, 10.0, 100.0]):
        x = np.linspace(0,5*e,1000)
        y0 = [fdfunc1(xi, 0.0, e, 0.0, use_timmes) for xi in x]
        y1 = [fdfunc1(xi, 0.0, e, 0.1, use_timmes) for xi in x]
        lab = r"$\beta \rightarrow 0$" if c=='r' else None
        plt.plot(x/e, y0, lw=1.2, ls='-.', c=c, label=lab)
        plt.plot(x/e, y1, lw=1.6, c=c, label=r"$\eta=%2.1f$"%e)

    plt.text(2.0, 1.25, r"$\beta=\frac{kT}{mc^2}=0.1$", fontsize=22)
    plt.title("Relativistic Fermi-Dirac functions")
    plt.xlabel(r"$\epsilon/\mu$", fontsize=18)
    plt.ylabel(r"$\frac{\sqrt{1+x\beta/2}}{e^{(\epsilon-\mu)/kT}+1}$", fontsize=24)
    plt.ylim(-0.2, 2.0)
    plt.legend(loc='best')
    plt.show()

