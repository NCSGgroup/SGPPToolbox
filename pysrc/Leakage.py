import numpy as np

from pysrc.GeoMathKit import GeoMathKit
from pysrc.Grid import Grid
from pysrc.Harmonic import Harmonic, EllipsoidType
from pysrc.SHC import SHC


class Leakage:
    def __init__(self):
        self.area = None
        self.grid_space = None
        self.signals = None
        self.signal_type = None  # Grid, SHC
        self.filter_matrix = None
        self.nmax = None
        self.scale_factor = None

    def setArea(self, area):
        if type(area) is Grid:
            self.grid_space = area.grid_space
            self.area = area.map
        else:
            self.grid_space = 180 / np.shape(area)[0]
            self.area = area

        return self

    def setSignal(self, signals: list):
        """
        :param signals:list of SHCs or list of Grids
        """
        assert type(signals[0]) in (SHC, Grid)

        self.signals = signals
        self.signal_type = type(signals[0])
        return self

    def setGaussianFilter(self, filter_matrix):
        self.filter_matrix = filter_matrix
        self.nmax = np.shape(filter_matrix)
        return self

    def getTS(self):

        area = self.area
        matrix = self.filter_matrix

        year_fracs = []
        leakage_signals = []
        signals = []
        area_bar = 1 - area

        lat = np.arange(-90, 90, self.grid_space)
        lon = np.arange(-180, 180, self.grid_space)
        Har = Harmonic(Parallel=-1).setEllipsoid(ell=EllipsoidType.gif48)
        Pnm = GeoMathKit.getPnm(lat, self.nmax, 1)

        grids_global = []

        Cnm_area, Snm_area = Har.analysis(60, area, lat, lon, Pnm)
        area_filtered = Har.synthesis(Cnm_area * matrix, Snm_area * matrix, 60, lat, lon)
        self.scale_factor = GeoMathKit.gridSum(area, area) / GeoMathKit.gridSum(area_filtered, area)

        if self.signal_type is SHC:
            SHCs = self.signals
            for i in range(len(SHCs)):
                shc = SHCs[i]
                grid_global = Har.synthesis_for_SHC(shc, lat, lon)
                grids_global.append(grid_global)

        elif self.signal_type is Grid:
            grids = self.signals
            for i in range(len(grids)):
                grids_global.append(grids[i])

        for i in range(len(grids_global)):
            year_fracs.append(grids_global[i].year_frac_average())

            this_map_bar = grids_global[i].map * area_bar
            Cnm_bar, Snm_bar = Har.analysis(60, this_map_bar, lat, lon, Pnm)
            Cnm_bar *= matrix
            Snm_bar *= matrix
            leakage_map = Har.synthesis(Cnm_bar, Snm_bar, 60, lat, lon)

            this_map = grids_global[i].map * area
            Cnm, Snm = Har.analysis(60, this_map, lat, lon, Pnm)
            Cnm *= matrix
            Snm *= matrix
            signal_map = Har.synthesis(Cnm, Snm, 60, lat, lon)

            leakage_signals.append(GeoMathKit.gridSum(leakage_map, area))
            signals.append(GeoMathKit.gridSum(signal_map, area))

        return np.array(year_fracs), self.scale_factor * (np.array(signals) - np.array(leakage_signals))
