import copy
import struct
import sys

import numpy as np
from scipy.linalg import block_diag

from pysrc.SHC import SHC


def _read_BIN(file, mode='packed'):
    '''
    Read the binary file containing symmetric/full or block diagonal matrices and associated vectors and parameters.

    Usage:
    dat = read_BIN(file)
    dat = read_BIN(file, mode = 'packed')
    dat = read_BIN(file, mode = 'full')

    Inputs:
    file -> [str] Input file

    Parameters:
    mode -> [optional, str, default = 'packed'] Available options are 'packed' or 'full' form for the filter matrices.
    If 'packed', the matrix remains in packed form (dat['pack1'] field). If 'full', the matrix expands to its full form (dat['mat1'] field).
    Warning: the 'full' option may cause excessive RAM memory use with large matrices.

    Outputs:
    dat -> [dic]: Dictionary with the file content

    Notice:
    This program is translated from the matlab/octave source code read_BIN.m written by Roelof Rietbroek 2016.
    For more information, please refer to https://github.com/strawpants/GRACE-filter
    '''
    if mode == 'packed':
        unpack = False
    elif mode == 'full':
        unpack = True  # unpack matrix in full size
    else:
        raise Exception("Only 'packed' or 'full' are avaliable.")

    # ckeck endian
    endian = sys.byteorder

    if endian == 'little':
        # open the binary file in little endian
        f = open(file, 'rb')
    else:
        raise Exception('The endian of the binary file is little, but the endian of OS is big.')

    dat = {}
    # read the data version and type from the binary file
    dat['version'] = f.read(8).decode().strip()
    dat['type'] = f.read(8).decode()
    dat['descr'] = f.read(80).decode().strip()

    for key in ['nints', 'ndbls', 'nval1', 'nval2']:
        dat[key] = struct.unpack('<I', f.read(4))[0]

    for key in ['pval1', 'pval2']:
        dat[key] = struct.unpack('<I', f.read(4))[0]

    dat['nvec'], dat['pval2'] = 0, 1
    dat['nread'], dat['nval2'] = 0, dat['nval1']

    # read additional nblocks parameter
    nblocks = struct.unpack('<i', f.read(4))[0]

    lists = f.read(dat['nints'] * 24).decode().split()
    for element in lists:
        dat[element] = struct.unpack('<i', f.read(4))[0]

    lists = f.read(dat['ndbls'] * 24).decode().replace(':', '').split()
    for element in lists:
        dat[element] = struct.unpack('<d', f.read(8))[0]

    # side description meta data
    lists = f.read(dat['nval1'] * 24).decode()
    dat['side1_d'] = [(lists[i:i + 24]).replace('         ', '') for i in range(0, len(lists), 24)]

    # type specific meta data
    dat['blockind'] = np.array(struct.unpack('<' + str(nblocks) + 'i', f.read(4 * nblocks)))

    dat['side2_d'] = dat['side1_d']

    # read matrix data
    npack1 = dat['pval1'] * dat['pval2']
    dat['pack1'] = np.array(struct.unpack('<' + str(npack1) + 'd', f.read(8 * npack1)))

    f.close()  # close file

    if not unpack: return dat

    sz = dat['blockind'][0]
    dat['mat1'] = dat['pack1'][:sz ** 2].reshape(sz, sz).T

    shift1 = shift2 = sz ** 2

    for i in range(1, nblocks):
        sz = dat['blockind'][i] - dat['blockind'][i - 1]
        shift2 = shift1 + sz ** 2
        dat['mat1'] = block_diag(dat['mat1'], dat['pack1'][shift1:shift2].reshape(sz, sz).T)
        shift1 = shift2
    del dat['pack1']

    return dat


def _getblock(filter_type: str):
    dir = '../data/ddk-data/'
    if filter_type == 'DDK1':
        Wbd = _read_BIN(dir + 'Wbd_2-120.a_1d14p_4')
    elif filter_type == 'DDK2':
        Wbd = _read_BIN(dir + 'Wbd_2-120.a_1d13p_4')
    elif filter_type == 'DDK3':
        Wbd = _read_BIN(dir + 'Wbd_2-120.a_1d12p_4')
    elif filter_type == 'DDK4':
        Wbd = _read_BIN(dir + 'Wbd_2-120.a_5d11p_4')
    elif filter_type == 'DDK5':
        Wbd = _read_BIN(dir + 'Wbd_2-120.a_1d11p_4')
    elif filter_type == 'DDK6':
        Wbd = _read_BIN(dir + 'Wbd_2-120.a_5d10p_4')
    elif filter_type == 'DDK7':
        Wbd = _read_BIN(dir + 'Wbd_2-120.a_1d10p_4')
    elif filter_type == 'DDK8':
        Wbd = _read_BIN(dir + 'Wbd_2-120.a_5d9p_4')
    else:
        raise Exception('Currently, only DDK1~DDK8 are feasible.')

    return Wbd


class DDKFilter:
    def __init__(self, filter_type: str):
        self.Wbd = _getblock(filter_type)
        pass

    def ApplyTo(self, shc: SHC):
        W = self.Wbd
        shc_filtered = copy.deepcopy(shc)

        cilm = np.array([shc.Cnm2d, shc.Snm2d])

        lmax = cilm.shape[1] - 1

        # Extract the minimum and maximum degree supported by the filter matrix
        lmaxfilt, lminfilt = W['Lmax'], W['Lmin']

        # Determine the output maximum degree (limited by either the filter or input data)
        lmaxout = min(lmax, lmaxfilt)

        # Reserve space for output (will have same size as input) and set to zero
        cilm_filtered = np.zeros_like(cilm)

        # Loop parameter indicating the previous block number and the end position in the packed matrix of the previous block
        lastblckind, lastindex = 0, 0

        # loop over the available blocks
        for iblk in range(W['Nblocks']):
            # Get the degree of the block from the block index
            degree = (iblk + 1) // 2

            # Break loop if the degrees of the block are larger than the degrees of the input
            if degree > lmaxout: break
            trig = (iblk + int(iblk > 0) + 1) % 2

            # Compute the size of the side of the stored block
            sz = W['blockind'][iblk] - lastblckind

            # Initialize the filter order block to a unit diagonal matrix
            blockn = np.identity(lmaxfilt + 1 - degree)

            # Minimum (stored) degree for this particular block (may be limited by the mininum degree supported by the filter)
            lminblk = max(lminfilt, degree)

            shift = lminblk - degree
            # unpack the stored filterblock (vector) in a fully occupied order block matrix
            blockn[shift:, shift:] = W['pack1'][lastindex:lastindex + sz ** 2].reshape(sz, sz).T

            # Filter the input coefficients (this is in fact just a matrix vector multiplication)
            if trig:
                cilm_filtered[0, degree:lmaxout + 1, degree] = np.dot(
                    blockn[:lmaxout + 1 - degree, :lmaxout + 1 - degree],
                    cilm[0, degree:lmaxout + 1, degree])
            else:
                cilm_filtered[1, degree:lmaxout + 1, degree] = np.dot(
                    blockn[:lmaxout + 1 - degree, :lmaxout + 1 - degree],
                    cilm[1, degree:lmaxout + 1, degree])

            # Prepare the loop variables for next block
            lastblckind = W['blockind'][iblk]
            lastindex = lastindex + sz ** 2

        shc_filtered.Cnm2d, shc_filtered.Snm2d = cilm_filtered[0], cilm_filtered[1]
        return shc_filtered


def demo():
    from pysrc.LoadL2Product import demo_load_by_year
    shcs = demo_load_by_year(2005, 2005, institute='GFZ')[0]
    shc = shcs[10]

    print(shc.begin_date)
    ddkfilter = DDKFilter('DDK8')
    shc_filtered = ddkfilter.ApplyTo(shc)
    pass


if __name__ == '__main__':
    demo()
