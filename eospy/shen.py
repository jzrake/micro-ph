

col_names = ['log10_rhoB', 'nB', 'Yp', 'F', 'Eint', 'S', 'A', 'Z',
                 'MN', 'Xn', 'Xp', 'Xa', 'XA', 'p', 'un', 'up', 'ML', 'XL']
var_index = { n:i for i,n in enumerate(col_names) }

def load_eos3(fname):
    """
    Loads the 'eos3' tabulated nuclear equation of state due to Shen from the
    ascii 'fname'.

    Returns:
    --------------------------------------------------------

    A (90 x 110 x 65 x 18) numpy array. The axes are (T,D,Yp,q) where q is the
    quantity given by the col_names list or var_index dictionary.

    Notes:
    --------------------------------------------------------

    user guide at http://user.numazu-ct.ac.jp/~sumi/eos/table2/guide_EOS3.pdf
    """
    import numpy as np

    iT, iD, iY = -1, -1, -1

    shen_file = open(fname, 'r')
    table_dims = (2, 110, 65, 18)
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
            data = [float(x) for x in line.split()]
            D = float(data[var_index['log10_rhoB']])
            iD += 1
            table[iT, iD, iY] = data

    return table


def write_hdf5(table, fname):
    import h5py
    h5f = h5py.File(fname, "w")

    for col in col_names:
        h5f[col] = table[:,:,:,var_index[col]]
