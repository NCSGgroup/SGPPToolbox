import numpy as np

from pysrc.Grid import Grid
from pysrc.GeoMathKit import GeoMathKit
from pysrc.RefEllipsoid import EllipsoidType, RefEllipsoid
from pysrc.SHC import SHC


class Harmonic:
    """
    Harmonic analysis and synthesis: Ordinary 2D integration for computing Spherical Harmonic coefficients
    """

    def __init__(self, Parallel: int = -1):
        self._parallel = Parallel
        self._ellipsoid = RefEllipsoid(EllipsoidType.gif48)
        pass

    def setEllipsoid(self, ell: EllipsoidType):
        self._ellipsoid = RefEllipsoid(ell)
        return self

    def analysis(self, Nmax: int, Inner, lat, lon, Pnm):
        if type(Inner) is Grid:
            Inner = Inner.map
        nlat = len(lat)
        nlon = len(lon)
        theta, phi = GeoMathKit.getCoLatLoninRad(lat, lon)

        term1 = np.zeros(Nmax + 1)

        for l in range(0, Nmax + 1):
            term1[l] = 1 + 2 * l

        NMmax = int((Nmax + 1) * (Nmax + 2) / 2)

        factor1 = 2 * np.pi / nlon
        factor2 = 0.25 / nlat
        Cnm, Snm = np.zeros(NMmax), np.zeros(NMmax)

        Am, Bm = np.zeros((Nmax + 1, nlat)), np.zeros((Nmax + 1, nlat))
        I_new = Inner.reshape(-1, nlon)

        for m in range(Nmax + 1):
            Am[m] = factor1 * np.array(I_new * np.mat(np.cos(m * phi)).T).flatten()
            Bm[m] = factor1 * np.array(I_new * np.mat(np.sin(m * phi)).T).flatten()

        thetaS = np.tile(np.sin(theta), (GeoMathKit.getIndex(Nmax, Nmax) + 1, 1))
        Qnm = Pnm * thetaS

        for n in range(Nmax + 1):
            indexM = np.arange(n + 1)
            Cnm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1] = factor2 * np.sum(
                Qnm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1] * Am[indexM], 1)
            Snm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1] = factor2 * np.sum(
                Qnm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1] * Bm[indexM], 1)
            pass
        return GeoMathKit.CS_1dTo2d(Cnm), GeoMathKit.CS_1dTo2d(Snm)

    def analysis_for_Grid(self, Nmax: int, Inner: Grid, Pnm):
        grid = Inner.map
        lat = Inner.lat
        lon = Inner.lon

        shc = SHC(self.analysis(Nmax, grid, lat, lon, Pnm))
        shc.begin_date, shc.end_date = Inner.begin_date, Inner.end_date
        shc.type = Inner.type

        return shc

    def job(self, n):
        """
        Do not call this function as it is only used for parallel computing.
        :param n:
        :return:
        """
        Nmax, nlat, nlon, I, factor1, phi, Cnm, Snm, factor2, Pnm, theta = self.__input

        Am, Bm = np.zeros((Nmax + 1, nlat)), np.zeros((Nmax + 1, nlat))
        I_new = I[n].reshape(-1, nlon)

        thetaS = np.tile(np.sin(theta), (GeoMathKit.getIndex(Nmax, Nmax) + 1, 1))
        Qnm = Pnm * thetaS

        for m in range(n + 1):
            Am[m] = factor1 * np.array(I_new * np.mat(np.cos(m * phi)).T).flatten()
            Bm[m] = factor1 * np.array(I_new * np.mat(np.sin(m * phi)).T).flatten()

        indexM = np.arange(n + 1)
        Cnm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1] = factor2 * np.sum(
            Qnm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1] * Am[indexM], 1)
        Snm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1] = factor2 * np.sum(
            Qnm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1] * Bm[indexM], 1)

        c = list(Cnm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1])
        s = list(Snm[GeoMathKit.getIndex(n, 0):GeoMathKit.getIndex(n, n) + 1])
        return n, c, s

    def synthesis(self, Cnm, Snm, Nmax, lat, lon):
        """
        A two step synthesis method, see the paper GJI (Nico Sneew)
        :param Cnm: in general, it should be the geo-potential coefficients sorted in 2 dimension.
        :param Snm:
        :param Nmax: Max degree of harmonic expansion.
        :param lat: geophysical latitude in unit "degree" [dimension : N]
        :param lon: geophysical latitude in unit "degree"[dimension : M]
        :return: grid (nlat*nlon) [dimension N*M]
        """

        lat, lon = GeoMathKit.getCoLatLoninRad(lat, lon)
        nlat = np.size(lat)
        nlon = np.size(lon)
        Pnm = GeoMathKit.getPnm(lat, Nmax)
        Am, Bm = np.zeros((Nmax + 1, nlat)), np.zeros((Nmax + 1, nlat))
        for m in range(Nmax + 1):
            for l in range(m, Nmax + 1):
                index = GeoMathKit.getIndex(l, m)
                Am[m] = Am[m] + Pnm[index] * Cnm[l][m]
                Bm[m] = Bm[m] + Pnm[index] * Snm[l][m]

        Fout = 0

        for m in range(Nmax + 1):
            co = np.cos(m * lon)
            so = np.sin(m * lon)
            Fout = Fout + np.mat(Am[m]).T * co + np.mat(Bm[m]).T * so

        return np.array(Fout)

    def synthesis_for_SHC(self, shc: SHC, lat, lon):
        Cnm, Snm = shc.Cnm2d, shc.Snm2d
        Nmax = shc.nmax

        grid = Grid(self.synthesis(Cnm, Snm, Nmax, lat, lon))
        grid.type = shc.type
        grid.begin_date, grid.end_date = shc.begin_date, shc.end_date

        return grid



