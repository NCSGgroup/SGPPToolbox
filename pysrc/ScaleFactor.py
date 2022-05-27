from pysrc.Decorrelation import *
from pysrc.GaussianFilter import *
from pysrc.GeoMathKit import GeoMathKit
from pysrc.Harmonic import Harmonic
from pysrc.SHC import SHC
from pysrc.Setting import *


class ScaleFactor:
    def __init__(self, area: np.ndarray):
        self.area = area
        self.grid_space = 180 / np.shape(area)[0]

        self.dec = None
        self.Gs = None
        pass

    def setFilter(self, de_correlation_filter=None, Gaussian_filter=None):
        self.dec = de_correlation_filter
        self.Gs = Gaussian_filter
        return self

    def getFactor(self, nmax=60):
        lat = np.arange(-90, 90, self.grid_space)
        lon = np.arange(-180, 180, self.grid_space)
        colat_rad, lon_rad = GeoMathKit.getCoLatLoninRad(lat, lon)
        PnmMatrix = GeoMathKit.getPnmMatrix(colat_rad, nmax)
        Har = Harmonic(lat, lon, PnmMatrix, nmax)

        CS = Har.analysis([self.area])
        Cnm = CS[0][0]
        Snm = CS[1][0]

        if self.dec is not None:
            shc_basin = SHC((Cnm, Snm))
            if self.dec[0] == DecorrelationFilterType.PnMm:
                dec = PnMm(*self.dec[1])
                shc_basin = dec.ApplyTo(shc_basin)

            elif self.dec[0] == DecorrelationFilterType.StableWindow:
                dec = StableWindow(*self.dec[1])
                shc_basin = dec.ApplyTo(shc_basin)

            elif self.dec[0] == DecorrelationFilterType.VariableWindow:
                dec = VariableWindow(*self.dec[1])
                shc_basin = dec.ApplyTo(shc_basin)

            Cnm, Snm = shc_basin.Cnm2d, shc_basin.Snm2d

        if self.Gs is not None:
            if self.Gs[0] == GaussianFilterType.isotropic:
                Gs_filter = IsoGaussian(*self.Gs[1])
                Cnm = Gs_filter.ApplyTo([Cnm])[0]
                Snm = Gs_filter.ApplyTo([Snm])[0]

            elif self.Gs[0] == GaussianFilterType.anisoropic:
                Gs_filter = AniGaussian(*self.Gs[1])
                Cnm = Gs_filter.ApplyTo([Cnm])[0]
                Snm = Gs_filter.ApplyTo([Snm])[0]

            elif self.Gs[0] == GaussianFilterType.fan:
                Gs_filter = Fan(*self.Gs[1])
                Cnm = Gs_filter.ApplyTo([Cnm])[0]
                Snm = Gs_filter.ApplyTo([Snm])[0]

        basin_filtered = Har.synthesis([Cnm], [Snm])[0]

        k = GeoMathKit.gridSum(self.area, self.area) / GeoMathKit.gridSum(basin_filtered, self.area)
        return k


def demo():
    basin = np.load('../data/grids/Caspian_maskGrid.dat(360,720).npy')
    sf = ScaleFactor(basin)
    sf.setFilter(de_correlation_filter=(DecorrelationFilterType.PnMm, (3, 10)),
                 Gaussian_filter=(GaussianFilterType.anisoropic, (300, 500, 15)))
    k = sf.getFactor()
    print(k)


if __name__ == '__main__':
    demo()
