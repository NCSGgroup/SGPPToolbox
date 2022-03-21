import copy

import numpy as np

from pysrc.GaussianFilter import IsoGaussian as isoGs, AniGaussian as aniGs, Fan as Fan
from pysrc.GeoMathKit import GeoMathKit
from pysrc.Harmonic import Harmonic, EllipsoidType
from pysrc.Setting import GaussianFilterType


class ForwardModeling:
    def __init__(self):
        self.model_obs = None
        self.model_tru = None
        self.initial = None
        self.filter = None
        self.gaussian_radius = None
        self.ani_radius_0, self.ani_radius_1 = None, None
        self.fan_radius_0, self.fan_radius_1 = None, None
        self.area = None
        self.thr_time = 100
        self.nmax = 60

    def setModel(self, observed, initial=None):
        if initial is None:
            initial = observed

        self.model_obs = observed
        self.initial = initial
        return self

    def setArea(self, area):
        self.area = area
        return self

    def setFilter(self, filter_method: GaussianFilterType, *rs):
        self.filter = filter_method

        if filter_method == GaussianFilterType.iso:
            self.gaussian_radius = rs[0]

        if filter_method == GaussianFilterType.ani:
            self.ani_radius_0 = rs[0]
            self.ani_radius_1 = rs[1]

        if filter_method == GaussianFilterType.fan:
            self.fan_radius_0 = rs[0]
            self.fan_radius_1 = rs[1]
        return self

    def setNmax(self, n):
        self.nmax = n
        return self

    def setLoopTime(self, time=None):
        if time is not None:
            self.thr_time = time
        return self

    def run(self):
        sp = self.model_obs.grid_space

        lat = np.arange(-90, 90, sp)
        lon = np.arange(-180, 180, sp)
        Pnm = GeoMathKit.getPnm(lat, self.nmax, 1)

        model_tru = copy.deepcopy(self.initial).map * self.area
        Har = Harmonic(Parallel=-1)
        Har.setEllipsoid(EllipsoidType.gif48)

        filter_mat = None
        if self.filter == GaussianFilterType.isotropic:
            filter_mat = isoGs(self.gaussian_radius).getFilter(self.nmax)

        if self.filter == GaussianFilterType.anisoropic:
            filter_mat = aniGs(self.ani_radius_0, self.ani_radius_1).getFilter(self.nmax)

        if self.filter == GaussianFilterType.fan:
            filter_mat = Fan(self.fan_radius_0, self.fan_radius_1).getFilter(self.nmax)

        it = 0
        while True:
            it += 1
            Cnm, Snm = Har.analysis(Nmax=60, Inner=model_tru, lat=lat, lon=lon, Pnm=Pnm)
            model_pre = Har.synthesis(Cnm * filter_mat, Snm * filter_mat, self.nmax, lat, lon)
            model_dif = (self.model_obs.map - model_pre) * self.area

            if it > self.thr_time:
                break

            model_dif = GeoMathKit.keepland_c(model_dif)
            model_tru += model_dif

            print('\riter={}'.format(it), end='')
        print()
        self.model_tru = model_tru
        return self
