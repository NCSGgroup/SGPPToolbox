from pysrc.Decorrelation import *
from pysrc.GaussianFilter import *
from pysrc.GeoMathKit import GeoMathKit
from pysrc.Harmonic import Harmonic
from pysrc.Setting import *


class Leakage:
    def __init__(self, area: np.ndarray):
        self.area = area
        self.grid_space = 180 / np.shape(area)[0]

        self.origin_dec = None
        self.origin_Gs = None

        self.Gaussian_filter = None
        pass

    def setOriginalFilter(self, de_correlation_filter=None, Gaussian_filter=None):
        self.origin_dec = de_correlation_filter
        self.origin_Gs = Gaussian_filter
        return self

    def setGaussianFilter(self, Gaussian_filter):
        self.Gaussian_filter = Gaussian_filter
        return self

    def ApplyTo(self, *shcs):
        """
        Currently Cnms and Snms should be DataType.EWH
        """
        assert self.Gaussian_filter is not None

        shcs_filtered = shcs
        nmax = shcs[0].nmax

        if self.origin_dec is not None:
            if self.origin_dec[0] == DecorrelationFilterType.PnMm:
                dec = PnMm(*self.origin_dec[1])
                shcs_filtered = dec.ApplyTo(*shcs)

            elif self.origin_dec[0] == DecorrelationFilterType.StableWindow:
                dec = StableWindow(*self.origin_dec[1])
                shcs_filtered = dec.ApplyTo(*shcs)

            elif self.origin_dec[0] == DecorrelationFilterType.VariableWindow:
                dec = VariableWindow(*self.origin_dec[1])
                shcs_filtered = dec.ApplyTo(*shcs)

        Cnms_filtered = [shcs_filtered[i].Cnm2d for i in range(len(shcs_filtered))]
        Snms_filtered = [shcs_filtered[i].Snm2d for i in range(len(shcs_filtered))]

        if self.origin_Gs is not None:
            if self.origin_Gs[0] == GaussianFilterType.isotropic:
                Gs_filter = IsoGaussian(*self.origin_Gs[1])
                Cnms_filtered = Gs_filter.ApplyTo(Cnms_filtered)
                Snms_filtered = Gs_filter.ApplyTo(Snms_filtered)

            elif self.origin_Gs[0] == GaussianFilterType.anisoropic:
                Gs_filter = AniGaussian(*self.origin_Gs[1])
                Cnms_filtered = Gs_filter.ApplyTo(Cnms_filtered)
                Snms_filtered = Gs_filter.ApplyTo(Snms_filtered)

            elif self.origin_Gs[0] == GaussianFilterType.fan:
                Gs_filter = Fan(*self.origin_Gs[1])
                Cnms_filtered = Gs_filter.ApplyTo(Cnms_filtered)
                Snms_filtered = Gs_filter.ApplyTo(Snms_filtered)

        lat = np.arange(-90, 90, self.grid_space)
        lon = np.arange(-180, 180, self.grid_space)
        colat_rad, lon_rad = GeoMathKit.getCoLatLoninRad(lat, lon)
        PnmMatrix = GeoMathKit.getPnmMatrix(colat_rad, nmax)
        Har = Harmonic(lat, lon, PnmMatrix, nmax)

        grids = Har.synthesis(Cnms_filtered, Snms_filtered)

        grids_out_basin = grids * (1 - self.area)

        Cnms_out_basin, Snms_out_basin = Har.analysis(grids_out_basin)
        if self.Gaussian_filter[0] == GaussianFilterType.isotropic:
            Gs_filter_second = IsoGaussian(*self.origin_Gs[1])
            Cnms_out_basin = Gs_filter_second.ApplyTo(Cnms_out_basin)
            Snms_out_basin = Gs_filter_second.ApplyTo(Snms_out_basin)

        elif self.Gaussian_filter[0] == GaussianFilterType.isotropic:
            Gs_filter_second = IsoGaussian(*self.origin_Gs[1])
            Cnms_out_basin = Gs_filter_second.ApplyTo(Cnms_out_basin)
            Snms_out_basin = Gs_filter_second.ApplyTo(Snms_out_basin)

        elif self.Gaussian_filter[0] == GaussianFilterType.isotropic:
            Gs_filter_second = IsoGaussian(*self.origin_Gs[1])
            Cnms_out_basin = Gs_filter_second.ApplyTo(Cnms_out_basin)
            Snms_out_basin = Gs_filter_second.ApplyTo(Snms_out_basin)

        grids_out_basin_filtered = Har.synthesis(Cnms_out_basin, Snms_out_basin)

        return np.array([GeoMathKit.gridSum(grids_out_basin_filtered[i], self.area) for i in
                         range(len(grids_out_basin_filtered))])


def demo():
    from main.getTimeSeries import load_GSM_by_year

    shcs = load_GSM_by_year(2005, 2015, L2ProductRelease.RL06, L2instituteType.CSR, GSMMaxDegree.degree60, toEWH=True)[
        1]

    lk = Leakage(np.load('../data/grids/Ocean_maskGrid.dat(360,720).npy'))
    lk.setOriginalFilter((PnMm, (3, 10)), (GaussianFilterType.isotropic, (200,)))
    lk.setGaussianFilter((GaussianFilterType.isotropic, (300,)))
    lk.ApplyTo(*shcs)


if __name__ == '__main__':
    demo()
