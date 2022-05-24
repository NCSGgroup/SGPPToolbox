import copy
import json

import numpy as np

from pysrc.SHC import SHC
from pysrc.Setting import EarthModel


def for_SHCs(shcs, matrix):
    shcs_filtered = []
    for i in range(len(shcs)):
        this_shc = copy.deepcopy(shcs[i])
        this_shc *= matrix
        shcs_filtered.append(this_shc)

    if len(shcs) == 1:
        shcs_filtered = shcs_filtered[0]

    return shcs_filtered


def for_vecCS(shcs, matrix):
    if len(shcs) == 1:
        shc = shcs[0]
        return [matrix * shc]
    else:
        shcs = np.array(shcs)
        return np.einsum('lm,plm->plm', matrix, shcs)


class IsoGaussian:
    def __init__(self, r, earth_model: EarthModel = EarthModel.general):
        """
        :param r: Gaussian filter radius in unit [km]
        """
        if earth_model is EarthModel.general:
            with open('../data/json/EarthModels.json', 'r+') as f:
                dic = json.load(f)['general']
        self.radius_e = dic['radius_e']

        self.radius = r * 1000

    def getFilter(self, nmax, output_format='2d'):

        w = np.zeros(nmax + 1)
        b = np.log(2) / (1 - np.cos(self.radius / self.radius_e))
        w[0] = 1
        w[1] = (1 + np.exp(-2 * b)) / (1 - np.exp(-2 * b)) - 1 / b
        for i in range(1, nmax):
            w[i + 1] = -(w[i] * (2 * i + 1)) / b + w[i - 1]

        assert output_format in ['2d', '1d'], 'arg \'format\' must be \'1d\' or \'2d\''
        if output_format == '2d':
            return np.array([[w[n] for i in range(nmax + 1)] for n in range(len(w))])

        elif output_format == '1d':
            return w

    def ApplyTo(self, shcs: list):
        if type(shcs[0]) is SHC:
            nmax = shcs[0].nmax
            matrix = self.getFilter(nmax)
            return for_SHCs(shcs, matrix)

        else:
            nmax = np.shape(shcs[0])[0] - 1
            matrix = self.getFilter(nmax)
            return for_vecCS(shcs, matrix)


class AniGaussian:
    def __init__(self, r0, r1, trunc_m=None, earth_model: EarthModel = EarthModel.general):
        """
        r0, r1 are set in unit [km]
        m is the truncated order
        """
        if earth_model is EarthModel.general:
            with open('../data/json/EarthModels.json', 'r+') as f:
                dic = json.load(f)['general']
        self.radius_e = dic['radius_e']

        self.r0, self.r1 = r0 * 1000, r1 * 1000
        self.trunc_m = trunc_m

    def getFilter(self, nmax):
        """

        :return:
        """
        if self.trunc_m is None:
            self.trunc_m = nmax
        assert self.trunc_m <= nmax

        rr = lambda mm: (self.r1 - self.r0) / self.trunc_m * mm + self.r0

        matrix = np.zeros((nmax + 1, nmax + 1))

        w = np.zeros(nmax + 1)

        for n in range(nmax + 1):
            for m in range(n + 1):
                if m <= self.trunc_m:
                    r = rr(m)
                    b = np.log(2) / (1 - np.cos(r / self.radius_e))
                    w[0] = 1
                    w[1] = (1 + np.exp(-2 * b)) / (1 - np.exp(-2 * b)) - 1 / b
                    for i in range(1, n):
                        w[i + 1] = -(w[i] * (2 * i + 1)) / b + w[i - 1]
                    matrix[n][m] = w[n]

        return matrix

    def ApplyTo(self, shcs: list):
        if type(shcs[0]) is SHC:
            nmax = shcs[0].nmax
            matrix = self.getFilter(nmax)
            return for_SHCs(shcs, matrix)

        else:
            nmax = np.shape(shcs[0])[0] - 1
            matrix = self.getFilter(nmax)
            return for_vecCS(shcs, matrix)


class Fan:
    def __init__(self, r1, r2, earth_model: EarthModel = EarthModel.general):
        """
        r1, r2 are set in unit [km]
        """
        if earth_model is EarthModel.general:
            with open('../data/json/EarthModels.json', 'r+') as f:
                dic = json.load(f)['general']
        self.radius_e = dic['radius_e']

        self.r1 = r1 * 1000
        self.r2 = r2 * 1000

    def getFilter(self, nmax):
        matrix = np.zeros((nmax + 1, nmax + 1))

        w1 = np.zeros(nmax + 1)
        b1 = np.log(2) / (1 - np.cos(self.r1 / self.radius_e))
        w1[0] = 1
        w1[1] = (1 + np.exp(-2 * b1)) / (1 - np.exp(-2 * b1)) - 1 / b1
        for i in range(1, nmax):
            w1[i + 1] = -(w1[i] * (2 * i + 1)) / b1 + w1[i - 1]

        w2 = np.zeros(nmax + 1)
        b2 = np.log(2) / (1 - np.cos(self.r2 / self.radius_e))
        w2[0] = 1
        w2[1] = (1 + np.exp(-2 * b2)) / (1 - np.exp(-2 * b2)) - 1 / b2
        for i in range(1, nmax):
            w2[i + 1] = -(w2[i] * (2 * i + 1)) / b2 + w2[i - 1]

        for n in range(nmax + 1):
            for m in range(n + 1):
                matrix[n][m] = w1[n] * w2[m]

        return matrix

    def ApplyTo(self, shcs: list):
        if type(shcs[0]) is SHC:
            nmax = shcs[0].nmax
            matrix = self.getFilter(nmax)
            return for_SHCs(shcs, matrix)

        else:
            nmax = np.shape(shcs[0])[0] - 1
            matrix = self.getFilter(nmax)
            return for_vecCS(shcs, matrix)
