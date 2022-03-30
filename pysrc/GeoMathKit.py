import gzip
import json

import numpy as np

from pysrc.Grid import Grid
from pysrc.Setting import EarthModel


class GeoMathKit:
    @staticmethod
    def un_gz(file_name):

        # acquire the filename and remove the postfix
        f_name = file_name.replace(".gz", "")
        # start uncompress
        g_file = gzip.GzipFile(file_name)
        # read uncompressed files and write down a copy without postfix
        open(f_name, "wb+").write(g_file.read())
        g_file.close()

    @staticmethod
    def CS_2dTo1d(CS: np.ndarray):
        """
        Transform the CS in 2-dimensional matrix to 1-dimemsion vectors
        example:
        00
        10 11
        20 21 22
        30 31 32 33           =>           00 10 11 20 21 22 30 31 32 33 ....

        :param CS:
        :return:
        """
        shape = np.shape(CS)
        assert len(shape) == 2
        index = np.nonzero(np.tril(np.ones(shape)))

        return CS[index]

    @staticmethod
    def CS_1dTo2d(CS: np.ndarray):
        """

        C00 C10 C20 =>     C00
                           C10 C20

        :param CS: one-dimension array
        :return: two dimension array
        """

        def index(N):
            n = (np.round(np.sqrt(2 * N))).astype(np.int) - 1
            m = N - (n * (n + 1) / 2).astype(np.int) - 1
            return n, m

        CS_index = np.arange(len(CS)) + 1
        n, m = index(CS_index)

        dim = index(len(CS))[0] + 1
        CS2d = np.zeros((dim, dim))
        CS2d[n, m] = CS

        return CS2d

    @staticmethod
    def getIndex(n: int, m: int):
        """
        index of Cnm at one-dimension array
        :param n: degree
        :param m: order
        :return:
        """
        assert m <= n

        return int(n * (n + 1) / 2 + m)

    @staticmethod
    def getCoLatLoninRad(lat, lon):
        '''
        :param lat: geophysical coordinate in degree
        :param lon: geophysical coordinate in degree
        :return: Co-latitude and longitude in rad
        '''

        theta = (90. - lat) / 180. * np.pi
        phi = lon / 180. * np.pi

        return theta, phi

    @staticmethod
    def getPnm(lat, Nmax: int, option=0):
        """
        get legendre function up to degree/order Nmax in Lat.
        :param lat: Co-latitude if option=0, unit[rad]; geophysical latitude if option = others, unit[degree]
        :param Nmax:
        :param option:
        :return:
        """

        if option != 0:
            lat = (90. - lat) / 180. * np.pi

        NMmax = int((Nmax + 1) * (Nmax + 2) / 2)

        if type(lat) is np.ndarray:
            Nsize = np.size(lat)
        else:
            Nsize = 1

        Pnm = np.zeros((NMmax, Nsize))

        Pnm[GeoMathKit.getIndex(0, 0)] = 1

        Pnm[GeoMathKit.getIndex(1, 1)] = np.sqrt(3) * np.sin(lat)

        '''For the diagonal element'''
        for n in range(2, Nmax + 1):
            Pnm[GeoMathKit.getIndex(n, n)] = np.sqrt((2 * n + 1) / (2 * n)) * np.sin(lat) * Pnm[
                GeoMathKit.getIndex(n - 1, n - 1)]

        for n in range(1, Nmax + 1):
            Pnm[GeoMathKit.getIndex(n, n - 1)] = np.sqrt(2 * n + 1) * np.cos(lat) * Pnm[
                GeoMathKit.getIndex(n - 1, n - 1)]

        for n in range(2, Nmax + 1):
            for m in range(n - 2, -1, -1):
                Pnm[GeoMathKit.getIndex(n, m)] = \
                    np.sqrt((2 * n + 1) / ((n - m) * (n + m)) * (2 * n - 1)) \
                    * np.cos(lat) * Pnm[GeoMathKit.getIndex(n - 1, m)] \
                    - np.sqrt((2 * n + 1) / ((n - m) * (n + m)) * (n - m - 1) * (n + m - 1) / (2 * n - 3)) \
                    * Pnm[GeoMathKit.getIndex(n - 2, m)]

        return Pnm

    @staticmethod
    def shrink(data, rows, cols):
        return data.reshape(rows, int(data.shape[0] / rows), cols, int(data.shape[1] / cols)).sum(axis=1).sum(axis=2)

    @staticmethod
    def gridSum(signal, area):
        if type(signal) is Grid:
            signal = signal.map

        if np.shape(signal) != np.shape(area):
            if np.shape(signal) == (360, 720) and np.shape(area) == (180, 360):
                signal = GeoMathKit.shrink(signal, 180, 360) / 4
            elif np.shape(signal) == (180, 360) and np.shape(area) == (360, 720):
                area = GeoMathKit.shrink(area, 180, 360) / 4
            else:
                raise Exception('shapes of signal and area grid are must be equal.')

        ac = GeoMathKit.acreage_like(signal)
        return np.sum(signal * area * ac) / np.sum(area * ac)

    @staticmethod
    def acreage_like(grid):
        grid_space = 180 / np.shape(grid)[0]
        thetas = np.array([[np.radians(i * grid_space)] * len(grid[0]) for i in range(len(grid))])

        return np.sin(thetas)

    @staticmethod
    def keepland_c(grid):
        """
        keep the land signal while satisfying (quality) conservation.
        :param grid: spacing MUST be 0.5 * 0.5(360 * 720) or 1 * 1(180 * 360) [deg]
        :return:
        """
        shp = np.shape(grid)
        assert shp == (360, 720) or shp == (
            180, 360), 'grid spacing MUST be 0.5 * 0.5(360 * 720) or 1 * 1(180 * 360) [deg].'
        if shp == (360, 720):
            land_mask = np.load('../data/grids/ocean_maskGrid.dat(360,720).npy')
            ocean_mask = np.load('../data/grids/land_maskGrid.dat(360,720).npy')
        else:
            land_mask = np.load('../data/grids/land_mask(180,360).npy')
            ocean_mask = np.load('../data/grids/ocean_mask(180,360).npy')

        land = grid * ocean_mask
        total_ocean = sum(sum(grid * land_mask))
        ocean_fix = total_ocean * land_mask / sum(sum(land_mask))
        return land + ocean_fix

    @staticmethod
    def getAcreage(area, earth_model=EarthModel.general):

        delta_theta = np.radians(180 / np.shape(area)[0])
        if earth_model is EarthModel.general:
            with open('../data/json/EarthModels.json', 'r+') as f:
                dic = json.load(f)['general']
        radius_e = dic['radius_e']

        ac = GeoMathKit.acreage_like(area) * radius_e ** 2 * delta_theta ** 2
        return np.sum(area * ac)

    @staticmethod
    def getGridIndex(lon, lat, gs):
        return int((lat + 90) / gs), int((lon + 180) / gs)

    @staticmethod
    def xyz2grd(xyz):
        """
        xyz to grid
        :param xyz: np.ndarray or list, [ [lon0, lat0, z0], [lon1, lat1, z1], ... ], lons and lats' spacing are equal.
        :return: np.ndarray, grid
        """
        gs = max(np.abs(xyz[0][0] - xyz[1][0]), np.abs(xyz[0][1] - xyz[1][1]))
        grid = np.zeros((int(180 / gs), int(360 / gs)))
        for i in xyz:
            lon, lat = i[0], i[1]
            l, m = GeoMathKit.getGridIndex(lon, lat, gs)
            grid[l][m] = i[2]
        return grid
