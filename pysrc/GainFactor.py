from pysrc.Decorrelation import *
from pysrc.GaussianFilter import IsoGaussian, Fan, AniGaussian
from pysrc.GeoMathKit import GeoMathKit
from pysrc.Harmonic import Harmonic
from pysrc.LeastSquare import curve_fit
from pysrc.SHC import SHC
from pysrc.Setting import *


def _func(x, k):
    return k * x


class GainFactor:
    def __init__(self):
        self.grid_space = None
        self.dec = None
        self.Gs = None

        self.original_models = None

    def setFilter(self, de_correlation_filter=None, Gaussian_filter=None):
        self.dec = de_correlation_filter
        self.Gs = Gaussian_filter
        return self

    def setModel(self, grids):
        self.original_models = grids - np.mean(grids, axis=0)
        self.grid_space = 180 / np.shape(grids[0])[0]

        return self

    def __filter(self, nmax):

        lat = np.arange(-90, 90, self.grid_space)
        lon = np.arange(-180, 180, self.grid_space)
        colat_rad, lon_rad = GeoMathKit.getCoLatLoninRad(lat, lon)
        PnmMatrix = GeoMathKit.getPnmMatrix(colat_rad, nmax)
        Har = Harmonic(lat, lon, PnmMatrix, nmax)

        Cnms, Snms = Har.analysis(self.original_models)

        if self.dec is not None:
            shcs = [SHC((Cnms[i], Snms[i])) for i in range(len(self.original_models))]

            if self.dec[0] == DecorrelatedFilterType.PnMm:
                dec = PnMm(*self.dec[1])
                shcs = dec.ApplyTo(*shcs)

            elif self.dec[0] == DecorrelatedFilterType.StableWindow:
                dec = StableWindow(*self.dec[1])
                shcs = dec.ApplyTo(*shcs)

            elif self.dec[0] == DecorrelatedFilterType.VariableWindow:
                dec = VariableWindow(*self.dec[1])
                shcs = dec.ApplyTo(*shcs)

            Cnms = [shcs[i].Cnm2d for i in range(len(shcs))]
            Snms = [shcs[i].Snm2d for i in range(len(shcs))]

        if self.Gs is not None:
            if self.Gs[0] == GaussianFilterType.isotropic:
                Gs_filter = IsoGaussian(*self.Gs[1])
                Cnms = Gs_filter.ApplyTo(Cnms)
                Snms = Gs_filter.ApplyTo(Snms)

            elif self.Gs[0] == GaussianFilterType.anisoropic:
                Gs_filter = AniGaussian(*self.Gs[1])
                Cnms = Gs_filter.ApplyTo(Cnms)
                Snms = Gs_filter.ApplyTo(Snms)

            elif self.Gs[0] == GaussianFilterType.fan:
                Gs_filter = Fan(*self.Gs[1])
                Cnms = Gs_filter.ApplyTo(Cnms)
                Snms = Gs_filter.ApplyTo(Snms)

        models_filtered = Har.synthesis(Cnms, Snms)

        return models_filtered

    def getBasinFactor(self, basin, nmax=60):

        if np.shape(basin) != np.shape(self.original_models[0]):
            basin = GeoMathKit.shrink(basin, *np.shape(self.original_models[0])) // 4

        models_filtered = self.__filter(nmax)

        basin_signals = [GeoMathKit.gridSum(self.original_models[i], basin) for i in
                         range(len(self.original_models))]
        basin_signals_filtered = [GeoMathKit.gridSum(models_filtered[i], basin) for i in
                                  range(len(models_filtered))]

        z = curve_fit(_func, np.array(basin_signals_filtered), np.array(basin_signals))
        k = z[0][0, 0]
        return k

    def getGridFactor(self, nmax=60):
        models_filtered = self.__filter(nmax)
        model_shape = np.shape(models_filtered[0])

        original_models_1d = np.array([self.original_models[i].flatten() for i in range(len(self.original_models))])
        models_filtered_1d = np.array([models_filtered[i].flatten() for i in range(len(models_filtered))])

        factors_grids = []
        for i in range(len(original_models_1d[0])):
            print('{}/{}'.format(i + 1, len(original_models_1d[0])))
            factors_grids.append(
                curve_fit(_func, np.array(models_filtered_1d[:, i]), np.array(original_models_1d[:, i]))[0][0, 0])

        factors_grids = np.array(factors_grids).reshape(model_shape)

        return factors_grids


def demo():
    import os
    from pysrc.LoadNoah import getTWS
    from pysrc.Plot import plot_grid
    from pysrc.Grid import Grid

    basin = np.load('../data/grids/Amazon_maskGrid.dat(360,720).npy')

    gldas_path = '../data/Noah2.1'
    files = os.listdir(gldas_path)
    files.sort()

    x = []
    grids = []
    for i in range(len(files)):
        yyyymm = files[i].split('.')[1]
        year = int(yyyymm[1:5])
        if not (2005 <= year <= 2019):
            continue
        month = int(yyyymm[5:])
        print(year, month)

        year_frac = year + (month - 1) / 12

        tws = getTWS(os.path.join(gldas_path, files[i]))

        x.append(year_frac)
        grids.append(tws)

    grids = np.array(grids)

    bf = GainFactor().setFilter(de_correlation_filter=(DecorrelatedFilterType.PnMm, (3, 10)),
                                Gaussian_filter=(GaussianFilterType.isotropic, (300,))).setModel(grids)
    k = bf.getGridFactor()
    plot_grid([Grid(k), '', 0, 2])
    print(k)


if __name__ == '__main__':
    demo()
