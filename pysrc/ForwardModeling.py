import copy

import numpy as np

from pysrc.GaussianFilter import IsoGaussian as isoGs, AniGaussian as aniGs, Fan as Fan
from pysrc.GeoMathKit import GeoMathKit
from pysrc.Grid import Grid
from pysrc.Harmonic import Harmonic
from pysrc.Setting import *


class ForwardModeling:
    def __init__(self):
        self.model_obs = None
        self.model_tru = None
        self.initials = None
        self.filter = None
        self.gaussian_radius = None
        self.ani_radius_0, self.ani_radius_1, self.ani_trunc_m = None, None, None
        self.fan_radius_0, self.fan_radius_1 = None, None
        self.area = None
        self.thr_time = 50
        self.nmax = None

    def setModel(self, observeds: list, initials=None):
        if initials is None:
            initials = observeds

        assert len(observeds) == len(initials)

        self.model_obs = observeds
        self.initials = initials
        return self

    def setArea(self, area):
        self.area = area
        return self

    def setFilter(self, filter_method: GaussianFilterType, *rs):
        self.filter = filter_method

        if filter_method == GaussianFilterType.isotropic:
            self.gaussian_radius = rs[0]

        if filter_method == GaussianFilterType.anisoropic:
            self.ani_radius_0 = rs[0]
            self.ani_radius_1 = rs[1]
            self.ani_trunc_m = rs[2]

        if filter_method == GaussianFilterType.fan:
            self.fan_radius_0 = rs[0]
            self.fan_radius_1 = rs[1]
        return self

    def setNmax(self, n):
        self.nmax = n
        return self

    def setLoopTime(self, time: int):
        self.thr_time = time
        return self

    def run(self):
        if type(self.initials[0]) is Grid:
            initials = np.array([self.initials[i].map for i in range(len(self.initials))])
            sp = self.model_obs.grid_space
        else:
            initials = self.initials
            sp = 180 / np.shape(initials[0])[0]

        lat = np.arange(-90, 90, sp)
        lon = np.arange(-180, 180, sp)
        colat_rad, lon_rad = GeoMathKit.getCoLatLoninRad(lat, lon)

        PnmMartix = GeoMathKit.getPnmMatrix(colat_rad, self.nmax, 0)

        if type(self.model_obs[0]) is Grid:
            model_obs = np.array([self.model_obs[i].map for i in range(len(self.model_obs))])
        else:
            model_obs = self.model_obs

        model_tru = copy.deepcopy(initials) * self.area
        Har = Harmonic(lat, lon, PnmMartix, self.nmax)

        filter = None
        if self.filter == GaussianFilterType.isotropic:
            filter = isoGs(self.gaussian_radius)

        if self.filter == GaussianFilterType.anisoropic:
            filter = aniGs(self.ani_radius_0, self.ani_radius_1, self.ani_trunc_m)

        if self.filter == GaussianFilterType.fan:
            filter = Fan(self.fan_radius_0, self.fan_radius_1)

        it = 0

        while True:
            it += 1

            print('\rforward modeling: iter={}'.format(it), end='')

            model_tru = GeoMathKit.keepland_c(model_tru)

            Cnms, Snms = Har.analysis(Inners=model_tru)
            Cnms_filtered, Snms_filtered = filter.ApplyTo(Cnms), filter.ApplyTo(Snms)

            model_pre = Har.synthesis(Cnms_filtered, Snms_filtered)
            model_dif = (model_obs - model_pre) * self.area

            model_tru += model_dif

            if it >= self.thr_time:
                break

        print()
        self.model_tru = model_tru
        return self
