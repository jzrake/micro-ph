
"""
 * -----------------------------------------------------------------------------
 * FILE: shen.py
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * DESCRIPTION:
 *
 * This code provides a class which reads the 3d tabulated nuclear equation of
 * state developed by Shen (2011). It selects a single, global value of the
 * proton fraction Yp, and reads a 2d slice of the table in logD := log10(rho)
 * and logT := log10(T).
 *
 * REFERENCES:
 *
 * http://user.numazu-ct.ac.jp/~sumi/eos/table2/guide_EOS3.pdf
 *
 * CAUTION:
 *
 * (1) For densities greater than 10^14 gm/cm^3 the pressure at constant density
 *     can become multi-valued in T, making the inverse lookup impossible.
 *
 * (2) The sound speed can become imaginary for high densities and low
 *     temperatures.
 *
 * CONVENTIONS:
 *
 * The unit conversion facility works by providing member functions such as
 * units.GramsPerCubicCentimeter(). If the variable D_phys contains the density
 * in gm/cm^3 units, then D_code := D_phys * units.GramsPerCubicCentimeter().
 *
 * All public member functions receive the density (D) in code units defined
 * through the PhysicalUnits class instance, and the log10 of temperature (T) in
 * MeV. They also return the pressure (p), internal energy (u) and sound speed
 * (squared, cs2) in code units, but the log10 of temperature in MeV. Private
 * member functions all receive log10 of the quantities in the physical units
 * used by Shen's table, which are MeV for temperature, gm/cm^3 for density, and
 * MeV/fm^3 for pressure and internal energy density.
 *
 * -----------------------------------------------------------------------------
"""


from scipy import ndimage
import numpy as np
import h5py
import physics



_col_names = ['logT', 'log10_rhoB', 'nB', 'Yp', 'F', 'Eint', 'S', 'A', 'Z',
             'MN', 'Xn', 'Xp', 'Xa', 'XA', 'p', 'un', 'up', 'ML', 'XL']
_var_index = { n:i for i,n in enumerate(_col_names) }


def load_eos3(fname):
    """
    Parses the ascii file 'fname' for the *eos3* tabulated nuclear equation of
    state due to Shen et. al. (1998). May take a couple minutes to run.

    Returns:
    --------------------------------------------------------

    A (110 x 91 x 65 x 19) numpy array. The axes are (D,T,Yp,q) where q is the
    quantity given by the _col_names list or _var_index dictionary.

    Notes:
    --------------------------------------------------------

    See the following URL for the eos3 user guide:

    http://user.numazu-ct.ac.jp/~sumi/eos/table2/guide_EOS3.pdf
    """
    iD, iT, iY = -1, -1, -1

    shen_file = open(fname, 'r')
    table_dims = (110, 91, 65, 19)
    table = np.zeros(table_dims)

    for line in shen_file:

        if len(line) <= 3:
            # the empty line after a Yp block
            iY += 1
            iD = -1
            continue

        elif line.startswith(' cccccccccccc'):
            action = 'log10(Temp)_header'

        elif action == 'log10(Temp)_header':
            action = 'read_log10(Temp)'

        elif action == 'read_log10(Temp)':
            action = 'read_data_line'
            logT, Tstr = line.split()
            T = float(Tstr)
            iT += 1
            iY = -1
            iD = -1
            print "parsing block with T=%f MeV" % T

        elif action == 'read_data_line':
            data = [float(logT)] + [float(x) for x in line.split()]
            iD += 1
            table[iD, iT, iY] = data

    return table



def write_hdf5(table, fname, cols="all"):
    """
    Writes the Shen table, represented as the numpy array 'table' to the hdf5
    file 'fname'. One dataset is created for each column requested in 'cols',
    which defaults to all available data.
    """
    print "writing Shen table to", fname

    h5f = h5py.File(fname, "w")

    cols = [c for c in _col_names if c in cols or cols == "all"]

    for col in _col_names:
        h5f[col] = table[:,:,:,_var_index[col]]

    h5f.close()
    return None



def read_hdf5(fname, cols="all"):
    # print "reading Shen table from", fname

    h5f = h5py.File(fname, "r")
    table = { }

    cols = [c for c in _col_names if c in cols or cols == "all"]

    for col in cols:
        table[col] = h5f[col].value

    h5f.close()
    return table



def rebuild_hdf5(ascii, h5, cols="all"):
    """
    Convenience function, eads the ASCII table from 'ascii' and writes 'cols' to
    hdf5 formatted 'h5'.
    """
    table = load_eos3(ascii)
    write_hdf5(table, h5, cols)



def sample(table, col, D, T, Y, order=3):
    """
    Samples the EOS table using interpolation.
    """

    nD, nT, nY = table['Yp'].shape

    logD0, logD1 = table['log10_rhoB'][0,0,0], table['log10_rhoB'][-1,0,0]
    logT0, logT1 = table['logT'      ][0,0,0], table['logT'      ][0,-1,0]
    Y0   ,    Y1 = table['Yp'        ][0,0,0], table['Yp'        ][0,0,-1]

    x = (np.log10(D) - logD0) / (logD1 - logD0) * (nD-1)
    y = (np.log10(T) - logT0) / (logT1 - logT0) * (nT-1)
    z = (Y - Y0) / (Y1 - Y0) * (nY-1)

    return ndimage.map_coordinates(table[col], [[x],[y],[z]], order=order)[0]



def derivative2(data, xs, ax):
    f0 = np.roll(data, +1, axis=ax)
    f1 = np.roll(data, -1, axis=ax)
    x0 = np.roll(xs, +1, axis=ax)
    x1 = np.roll(xs, -1, axis=ax)
    return (f1 - f0) / (x1 - x0)



def derivative5(data, xs, ax):
    f0 = np.roll(data, -2, axis=ax)
    f1 = np.roll(data, -1, axis=ax)
    f3 = np.roll(data, +1, axis=ax)
    f4 = np.roll(data, +2, axis=ax)

    x0 = np.roll(xs, -2, axis=ax)
    x1 = np.roll(xs, -1, axis=ax)
    x3 = np.roll(xs, +1, axis=ax)
    x4 = np.roll(xs, +2, axis=ax)

    return (-f0 + 8*f1 - 8*f3 + f4) / (-x0 + 8*x1 - 8*x3 + x4)



if __name__ == "__main__":
    """
    Runs some tests.
    """
    import unittest
    print "testing", __file__

    class TestShenSampling(unittest.TestCase):

        def test_splines(self):
            table = read_hdf5("data/eos3.h5", cols=['logT', 'log10_rhoB', 'Yp', 'nB'])

            for i in range(0,table['Yp'].shape[0],10):
                D = 10**table['log10_rhoB'][i,0,0]
                T = 10**table['logT'      ][i,0,0]
                Y =     table['Yp'        ][i,0,0]
                n =     table['nB'        ][i,0,0]
                self.assertAlmostEqual(n, sample(table, 'nB', D, T, Y))

            for i in range(0,table['Yp'].shape[1],10):
                D = 10**table['log10_rhoB'][0,i,0]
                T = 10**table['logT'      ][0,i,0]
                Y =     table['Yp'        ][0,i,0]
                n =     table['nB'        ][0,i,0]
                self.assertAlmostEqual(n, sample(table, 'nB', D, T, Y))

            for i in range(0,table['Yp'].shape[2],10):
                D = 10**table['log10_rhoB'][0,0,i]
                T = 10**table['logT'      ][0,0,i]
                Y =     table['Yp'        ][0,0,i]
                n =     table['nB'        ][0,0,i]
                self.assertAlmostEqual(n, sample(table, 'nB', D, T, Y))

        def test_derivatives(self):
            """
            Examines the difference in the derivative estimate when taken by
            finite differencing the table directly, versus spline-sampling
            nearby points.
            """
            table = read_hdf5("data/eos3.h5", cols=['logT', 'log10_rhoB', 'Yp', 'p'])

            dpdT2 = derivative2(table['p'], 10**table['logT'], 1)
            dpdT5 = derivative5(table['p'], 10**table['logT'], 1)

            D = 10**table['log10_rhoB'][15,15,15]
            T = 10**table['logT'      ][15,15,15]
            Y =     table['Yp'        ][15,15,15]
            p =     table['p'         ][15,15,15]

            T0 = T*(1-1e-8)
            T1 = T*(1+1e-8)
            p0 = sample(table, 'p', D, T0, Y)
            p1 = sample(table, 'p', D, T1, Y)

            # multiplying by T/p is just to normalize the result near 1.0
            # self.assertAlmostEqual((T/p) * (p1-p0) / (T1-T0), (T/p)*dpdT[15,15,15], places=2)
            print (T/p) * (p1-p0) / (T1-T0), (T/p)*dpdT2[15,15,15], (T/p)*dpdT5[15,15,15]

        def test_free_energy(self):
            """
            Verifies equation (82) in the Shen user guide.
            """
            t = read_hdf5("data/eos3.h5", cols=['p', 'F', 'nB', 'un', 'up', 'Yp'])
            C = (t['un']*(1-t['Yp']) + t['up']*t['Yp'] - t['p']/t['nB']) - t['F']
            print "number of times F-consistency violated by more than 1%:",\
                len(t['F'][np.where(abs(C/t['F']) > 1e-2)])

        def test_pressure(self):
            """
            Uses the 5-point derivative to check the definition of the pressure
            in terms of the Helmholtz free energy.
            """
            table = read_hdf5("data/eos3.h5", cols=['p', 'F', 'nB', 'Eint', 'logT', 'S'])
            dFdnB = derivative5(table['F'], table['nB'], 0)
            print (dFdnB * table['nB']**2)[10,10,10], table['p'][10,10,10]

        def test_consistency(self):
            """
            Uses the 5-point derivative to check the thermodynamic consistency
            relation

            p = rho^2 d(u)/d(rho) + T dp/dT.

            """
            table = read_hdf5("data/eos3.h5", cols=['p', 'Eint', 'logT', 'nB'])
            dEdnB = derivative5(table['Eint'], table['nB'], 0)
            dpdlogT = derivative5(table['p'], table['logT'], 1)
            print "p consistency:", (table['p'] - dEdnB * table['nB']**2 - dpdlogT)[10,10,10]

    unittest.main()

