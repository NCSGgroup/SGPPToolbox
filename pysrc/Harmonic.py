import numpy as np

from pysrc.GeoMathKit import GeoMathKit


class Harmonic:
    """
    Harmonic analysis and synthesis: Ordinary 2D integration for computing Spherical Harmonic coefficients
    """

    def __init__(self, lat, lon, PnmMatrix, Nmax: int):
        self.lat, self.lon = GeoMathKit.getCoLatLoninRad(lat, lon)
        self.nlat, self.nlon = len(lat), len(lon)
        self.PnmMatrix = PnmMatrix

        m = np.arange(Nmax + 1)
        self.g = m[:, None] @ self.lon[None, :]

        self.factor1 = np.ones((self.nlat, Nmax + 1))
        self.factor1[:, 0] += 1
        self.factor1 = 1 / (self.factor1 * self.nlon)

        self.factor2 = np.ones((Nmax + 1, Nmax + 1))
        self.factor2[:, 0] += 1
        self.factor2 *= np.pi / (2 * self.nlat)
        pass

    def analysis(self, Inners: list):
        Inners = np.array(Inners)

        g = self.g.T
        co = np.cos(g)
        so = np.sin(g)

        Am = np.einsum('pij,jm->pim', Inners, co, optimize='greedy') * self.factor1
        Bm = np.einsum('pij,jm->pim', Inners, so, optimize='greedy') * self.factor1

        Cnms = np.einsum('pim,ilm,i->plm', Am, self.PnmMatrix, np.sin(self.lat), optimize='greedy') * self.factor2
        Snms = np.einsum('pim,ilm,i->plm', Bm, self.PnmMatrix, np.sin(self.lat), optimize='greedy') * self.factor2
        return Cnms, Snms

    def synthesis(self, Cnms: list, Snms: list):
        assert len(Cnms) == len(Snms)
        Cnms = np.array(Cnms)
        Snms = np.array(Snms)

        Am = np.einsum('ijk,ljk->ilk', Cnms, self.PnmMatrix)
        Bm = np.einsum('ijk,ljk->ilk', Snms, self.PnmMatrix)

        co = np.cos(self.g)
        so = np.sin(self.g)

        Fout = Am @ co + Bm @ so

        return Fout
