import numpy as np
from scipy import optimize

from pysrc.Decorrelation import PnMm, VariableWindow, StableWindow
from pysrc.GaussianFilter import IsoGaussian, Fan, AniGaussian
from pysrc.GeoMathKit import GeoMathKit
from pysrc.Harmonic import Harmonic
from pysrc.SHC import SHC
from pysrc.Setting import *


class SingleFactor:
    def __init__(self, nmax, area):
        self.ts_true = None
        self.ts_filtered = None
        self.nmax = nmax
        self.Gaussian_type = None
        self.Gaussian_paras = None
        self.dec_type = None
        self.dec_paras = None
        self.area = area
        self.vf = None
        self.vt = None
        self.k = None
        self.gs = None

    @staticmethod
    def __fc(x, a):
        return a * x

    def grids(self, grids):
        assert np.shape(grids[0]) in [(360, 720), (180, 360)]

        self.ts_true = grids
        if np.shape(grids[0]) == (360, 720):
            self.gs = 0.5
        else:
            self.gs = 1
        return self

    def setDecFilter(self, dec_filter_type: DecorrelatedFilterType, *paras):
        self.dec_type = dec_filter_type
        self.dec_paras = paras

        # self.poly_n = paras[0]
        # self.start_order = paras[1]
        #
        # if dec_filter_type == DecorrelatedFilterType.StableWindow:
        #     self.window_len = paras[3]
        #
        # elif dec_filter_type == DecorrelatedFilterType.SlideWindow:
        #     self.window_min = paras[2]
        #     self.window_a = paras[3]
        #     self.window_k = paras[4]
        # elif dec_filter_type == DecorrelatedFilterType.PnMm
        #
        # else:
        #     pass

        return self

    def setGaussianFilter(self, gaussian_filter_type:GaussianFilterType, *paras):
        self.Gaussian_type = gaussian_filter_type
        self.Gaussian_paras = paras

        if Gaussian_type == GaussianFilterType.isotropic:
            self.r = paras[0]

        elif Gaussian_type == GaussianFilterType.fan:
            self.r1 = paras[0]
            self.r2 = paras[1]

        elif Gaussian_type == GaussianFilterType.anisoropic:
            self.r0 = paras[0]
            self.r1 = paras[1]
            self.trunc_m = paras[2]

        return self

    def run(self):
        ts_true = self.ts_true
        average_true = sum(ts_true) / len(ts_true)
        delta_ts_true = np.array([ts_true[i] - average_true for i in range(len(ts_true))])
        lat = np.arange(-90, 90, self.gs)
        lon = np.arange(-180, 180, self.gs)
        colat_rad, lon_rad = GeoMathKit.getCoLatLoninRad(lat, lon)
        PnmMatrix = GeoMathKit.getPnmMatrix(colat_rad, 60, 0)
        Har = Harmonic(lat, lon, PnmMatrix, self.nmax)

        if self.Gaussian_type == GaussianFilterType.isotropic:
            filter_mat = IsoGaussian(self.r).getFilter(self.nmax)

        elif self.Gaussian_type == GaussianFilterType.fan:
            filter_mat = Fan(self.r1, self.r2).getFilter(self.nmax)

        elif self.Gaussian_type == GaussianFilterType.anisoropic:
            filter_mat = AniGaussian(self.r0, self.r1, self.trunc_m).getFilter(self.nmax)

        else:
            filter_mat = np.ones((self.nmax + 1, self.nmax + 1))

        delta_ts_filtered = []

        for i in range(len(delta_ts_true)):
            print('\rtime series:', i, end='')

            Cnm, Snm = Har.analysis(Nmax=60, Inner=delta_ts_true[i], lat=lat,
                                    lon=lon, Pnm=Pnm)
            shc = SHC(CS=(Cnm, Snm))

            if self.dec_type == DecorrelatedFilterType.PnMm:
                poly_n = self.dec_paras[0]
                start_order = self.dec_paras[1]

                pnmm = PnMm(poly_n, start_order).ApplyTo(shc)
                Cnm, Snm = pnmm.Cnm2d, pnmm.Snm2d

            if self.dec_type == DecorrelatedFilterType.StableWindow:
                poly_n = self.poly_n
                start_order = self.start_order
                window_len = self.window_len

                window_dec = StableWindow(poly_n, start_order, window_len).ApplyTo(shc)
                Cnm, Snm = window_dec.Cnm2d, window_dec.Snm2d

            if self.dec_type == DecorrelatedFilterType.SlideWindow:
                poly_n = self.poly_n
                start_order = self.start_order
                window_a = self.window_a
                window_k = self.window_k
                window_min = self.window_min

                window_dec = SlideWindow(poly_n, start_order, window_min, window_a, window_k).ApplyTo(shc)
                Cnm, Snm = window_dec.Cnm2d, window_dec.Snm2d

            signal_filtered = Har.synthesis(Cnm * filter_mat, Snm * filter_mat, self.nmax, lat, lon)
            delta_ts_filtered.append(signal_filtered)

        delta_ts_filtered = np.array(delta_ts_filtered)

        self.vf = np.array(
            [GeoMathKit.gridSum(delta_ts_filtered[i], self.area) for i in range(len(delta_ts_filtered))])
        self.vt = np.array(
            [GeoMathKit.gridSum(delta_ts_true[i], self.area) for i in range(len(delta_ts_filtered))])
        self.k = optimize.curve_fit(self.__fc, self.vf, self.vt)[0][0]

        return self


class GridFactor:
    def __init__(self):
        self.ts_true = None
        self.ts_filtered = None
        self.nmax = None
        self.DecFilterType = None
        self.GaussianFilterType = None
        self.window_a = None
        self.window_k = None
        self.window_min = None
        self.window_len = None
        self.poly_n = None
        self.start_order = None
        self.r = None
        self.r0 = None
        self.r1 = None
        self.r2 = None
        self.trunc_m = None
        self.area = None
        self.gs = None

    @staticmethod
    def __fc(x, a):
        return a * x

    @staticmethod
    def __curve_fit_maps(filtered_maps, true_maps):
        factor_map = np.zeros_like(filtered_maps[0])
        for i in range(len(filtered_maps[0])):
            print('calculating factor map... {}/360'.format(i))
            for j in range(len(filtered_maps[0][i])):
                v_f = []
                v_t = []
                for k in range(len(filtered_maps)):
                    v_f.append(filtered_maps[k][i][j])
                    v_t.append(true_maps[k][i][j])
                factor_map[i][j] = optimize.curve_fit(GridFactor.__fc, v_f, v_t)[0][0]
        return factor_map

    def setTimeSeries(self, timeseries):
        self.ts_true = timeseries
        assert np.shape(timeseries[0]) in [(360, 720), (180, 360)]
        if np.shape(timeseries[0]) == (360, 720):
            self.gs = 0.5
        else:
            self.gs = 1
        return self

    def setNmax(self, n):
        self.nmax = n
        return self

    def setDecFilter(self, dec_type: DecorrelatedFilterType, paras: dict):
        self.dec_type = dec_type
        self.poly_n = paras['n']
        self.start_order = paras['m']

        if dec_type == DecorrelatedFilterType.SlideWindow:
            self.window_a = paras['a']
            self.window_k = paras['k']
            self.window_min = paras['min']

        if dec_type == DecorrelatedFilterType.StableWindow:
            self.window_len = paras['len']

        if dec_type == DecorrelatedFilterType.PnMm:
            pass

        return self

    def setGaussianFilter(self, Gaussian_type: GaussianFilterType, paras: dict):
        self.Gaussian_type = Gaussian_type

        if Gaussian_type == GaussianFilterType.Iso:
            self.r = paras['r']

        if Gaussian_type == GaussianFilterType.Fan:
            self.r1 = paras['r1']
            self.r2 = paras['r2']

        if Gaussian_type == GaussianFilterType.Ani:
            self.r0 = paras['r0']
            self.r1 = paras['r1']
            self.trunc_m = paras['trunc_m']

        return self

    def getKMap(self):
        ts_true = self.ts_true
        average_true = sum(ts_true) / len(ts_true)
        delta_ts_true = np.array([ts_true[i] - average_true for i in range(len(ts_true))])
        lat = np.arange(-90, 90, self.gs)
        lon = np.arange(-180, 180, self.gs)
        Pnm = GeoMathKit.getPnm(lat, 60, 1)
        Har = Harmonic(Parallel=-1)
        Har.setEllipsoid(EllipsoidType.gif48)

        if self.Gaussian_type == GaussianFilterType.Iso:
            filter_mat = IsoGaussian(self.r).getFilter(self.nmax)

        elif self.Gaussian_type == GaussianFilterType.Fan:
            filter_mat = Fan(self.r1, self.r2).getFilter(self.nmax)

        elif self.Gaussian_type == GaussianFilterType.Ani:
            filter_mat = AniGaussian(self.r0, self.r1).getFilter(self.nmax)

        else:
            filter_mat = np.ones((self.nmax + 1, self.nmax + 1))

        delta_ts_filtered = []

        for i in range(len(delta_ts_true)):
            print('time series:', i)

            Cnm, Snm = Har.analysis(Nmax=60, Inner=delta_ts_true[i], lat=lat,
                                    lon=lon, Pnm=Pnm)
            shc = SHC(CS=(Cnm, Snm))

            if self.dec_type == DecorrelatedFilterType.PnMm:
                poly_n = self.poly_n
                start_order = self.start_order

                pnmm = PnMm(poly_n, start_order).ApplyTo(shc)
                Cnm, Snm = pnmm.Cnm2d, pnmm.Snm2d

            if self.dec_type == DecorrelatedFilterType.StableWindow:
                poly_n = self.poly_n
                start_order = self.start_order
                window_len = self.window_len

                window_dec = StableWindow(poly_n, start_order, window_len).ApplyTo(shc)
                Cnm, Snm = window_dec.Cnm2d, window_dec.Snm2d

            if self.dec_type == DecorrelatedFilterType.SlideWindow:
                poly_n = self.poly_n
                start_order = self.start_order
                window_a = self.window_a
                window_k = self.window_k
                window_min = self.window_min

                window_dec = SlideWindow(poly_n, start_order, window_min, window_a, window_k).ApplyTo(shc)
                Cnm, Snm = window_dec.Cnm2d, window_dec.Snm2d

            signal_filtered = Har.synthesis(Cnm * filter_mat, Snm * filter_mat, self.nmax, lat, lon)
            delta_ts_filtered.append(signal_filtered)

        delta_ts_filtered = np.array(delta_ts_filtered)
        kmap = self.__curve_fit_maps(delta_ts_filtered, delta_ts_true)
        return kmap
