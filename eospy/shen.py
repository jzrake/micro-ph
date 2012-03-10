
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


col_names = ['logT', 'log10_rhoB', 'nB', 'Yp', 'F', 'Eint', 'S', 'A', 'Z',
             'MN', 'Xn', 'Xp', 'Xa', 'XA', 'p', 'un', 'up', 'ML', 'XL']
var_index = { n:i for i,n in enumerate(col_names) }


def load_eos3(fname):
    """
    Loads the 'eos3' tabulated nuclear equation of state due to Shen from the
    ascii 'fname'.

    Returns:
    --------------------------------------------------------

    A (110 x 91 x 65 x 19) numpy array. The axes are (D,T,Yp,q) where q is the
    quantity given by the col_names list or var_index dictionary.

    Notes:
    --------------------------------------------------------

    See the following URL for the eos3 user guide:

    http://user.numazu-ct.ac.jp/~sumi/eos/table2/guide_EOS3.pdf
    """
    import numpy as np

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
            D = float(data[var_index['log10_rhoB']])
            iD += 1
            table[iD, iT, iY] = data

    return table


def write_hdf5(table, fname, cols="all"):
    """
    Writes the Shen table, represented as the numpy array 'table' to the hdf5
    file 'fname'. One dataset is created for each column requested in 'cols',
    which defaults to all available data.
    """
    import h5py
    print "writing Shen table to", fname

    h5f = h5py.File(fname, "w")

    cols = [c for c in col_names if c in cols or cols == "all"]

    for col in col_names:
        h5f[col] = table[:,:,:,var_index[col]]


def read_hdf5(fname, cols="all"):
    import h5py
    print "reading Shen table from", fname

    h5f = h5py.File(fname, "r")
    table = { }

    cols = [c for c in col_names if c in cols or cols == "all"]

    for col in cols:
        table[col] = h5f[col].value

    return table


def rebuild_hdf5_table(ascii, h5, cols="all"):
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
    logD0, logD1 = table['log10_rhoB'][0,0,0], table['log10_rhoB'][-1,0,0]
    logT0, logT1 = table['logT'      ][0,0,0], table['logT'      ][0,-1,0]
    Y0   ,    Y1 = table['Yp'        ][0,0,0], table['Yp'        ][0,0,-1]

    x = (np.log10(D) - logD0) * (logD1 - logD0)
    y = (np.log10(T) - logT0) * (logT1 - logT0)
    z = (Y - Y0) * (Y1 - Y0)

    return ndimage.map_coordinates(table[col], [[x],[y],[z]], order=order)
