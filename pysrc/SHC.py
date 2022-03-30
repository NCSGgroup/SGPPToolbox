import copy
import json

import numpy as np

from pysrc.GeoMathKit import GeoMathKit
from pysrc.Setting import DataType, EarthModel
from pysrc.TimeTool import TimeTool


class SHC:
    """
    This class define a data type SHC(Spherical Harmonic Coefficients), usually given by LoadL2Product.py.

    A whole variable of class SHC usually includes these attributes(* means not necessary):
        nmax: int, Maximum degree/order of this coefficients Cnm and Snm.
        Cnm2d / Snm2d: np.ndarray, in 2d, e.g. [[1, 0, ...], [0, 0, ...], ...].
        *earth_model: EarthModel(Enum).
        *type: Datatype(Enum), to show which geophysics signal this SHC stands for, e.g. DataType.geoid
        *begin_date / *end_date: date,date, Of which time this SHC signal is.
        *unused_days: list of date.date, Of which dates that GRACE(-FO) missed in this SHC.
    """

    def __init__(self, CS=None, earth_model=None):
        self.earth_model = EarthModel.general
        self.nmax = None
        self.type = None
        self.begin_date, self.end_date = None, None  # date.date
        self.Cnm2d, self.Snm2d = None, None
        self.unused_days = []

        if CS is not None:
            self.Cnm2d, self.Snm2d = CS[0], CS[1]
            self.nmax = np.shape(self.Cnm2d)[0] - 1

        if earth_model is not None:
            self.earth_model = earth_model

        pass

    def __str__(self):
        return str([self.Cnm2d, self.Snm2d])

    def __add__(self, other):
        assert type(other) is SHC
        assert self.type == other.type and self.nmax == other.nmax

        new = copy.deepcopy(self)

        new.begin_date = min(self.begin_date, other.begin_date)
        new.end_date = max(self.end_date, other.end_date)

        new.Cnm2d = self.Cnm2d + other.Cnm2d
        new.Snm2d = self.Snm2d + other.Snm2d

        return new

    def __sub__(self, other):
        assert type(other) is SHC
        assert self.type == other.type and self.nmax == other.nmax

        new = copy.deepcopy(self)
        new.Cnm2d = self.Cnm2d - other.Cnm2d
        new.Snm2d = self.Snm2d - other.Snm2d
        return new

    def __truediv__(self, other):
        new = copy.deepcopy(self)
        new.Cnm2d /= other
        new.Snm2d /= other
        return new

    def __mul__(self, other: np.ndarray):
        new = copy.deepcopy(self)
        assert np.shape(new.Cnm2d) == np.shape(other)

        new.Cnm2d *= other
        new.Snm2d *= other
        return new

    @staticmethod
    def mean(iters):
        """
        Get the average of some SHCs(usually time continuous)
        :param iters: iters of SHC.
        :return: SHC of the average.
        Note:
            1. Every SHC should be in the same type(e.g. DataType.geoid) and the same size(i.e. nmax).
            2. begin_date of the result is the most early one of inputs, and end_date of which is the latest one.
        """
        result = iters[0]
        for i in range(len(iters) - 1):
            result += iters[i + 1]
        return result / len(iters)

    @staticmethod
    def delta(SHC_list, background):
        """
        Get the SHC anomaly relative to a background one.
        :param SHC_list: list of SHC
        :param background: background SHC
        :return: list of each anomaly of int inputs relative to the background.
        Note:
            1. Every SHC should be in the same type and the same size with the background one.
        """
        return [SHC_list[i] - background for i in range(len(SHC_list))]

    def simple_begin(self):
        """

        :return: str, The begin date in form of yyyymm, e.g. '200805'.
        """
        return TimeTool.yyyymm(self.begin_date)

    def year_frac_begin(self):
        """

        :return: float, The begin date in form of year-frac, e.g. 2008.331 stands for date 2008-5-1.
        """
        return TimeTool.year_frac(self.begin_date)

    def simple_end(self):
        """

        :return: str, The end date in form of yyyymm.
        """
        return TimeTool.yyyymm(self.end_date)

    def year_frac_end(self):
        """

        :return: float, The begin date in form of year-frac, e.g. 2008.331 stands for date 2008-5-1.
        """
        return TimeTool.year_frac(self.end_date)

    def average_date(self, simple_mode=False):
        """

        :param simple_mode: bool, True for returning the result considering the unused dates, False for not considering them.
        :return: date.date, the average date from beginning date to ending date.
        """
        if simple_mode or not self.unused_days:
            return self.begin_date + (self.end_date - self.begin_date) / 2
        else:
            used_days = TimeTool.list_of_unsed_days(self.begin_date, self.end_date, self.unused_days)
            return TimeTool.average_date_of_series(used_days)

    def simple_average(self, simple_mode=False):
        """

        :param simple_mode: bool, True for returning the result considering the unused dates, False for not considering them.
        :return: str, the average date in yyyymm format from beginning date to ending date.
        """
        return TimeTool.yyyymm(self.average_date(simple_mode=simple_mode))

    def year_frac_average(self, simple_mode=False):
        """

        :param simple_mode: bool, True for returning the result considering the unused dates, False for not considering them.
        :return: float, the average date in year-fraction format from beginning date to ending date.
        """
        return TimeTool.year_frac(self.average_date(simple_mode=simple_mode))

    def convertTypeTo(self, datatype: DataType, ln):
        """
        Convert the data type to another one. Currently, this function supports only DataType.geoid to another.
        :param datatype: DataType, the conversion target.
        :param ln: load love number, given by LoveNumber.LoveNumber.
        :return:
        """
        assert len(ln) >= self.nmax
        assert self.type is DataType.geoid

        if self.earth_model is EarthModel.general:
            with open('../data/json/EarthModels.json', 'r+') as f:
                dic = json.load(f)['general']
        density_water, density_earth, radius_e = dic['density_water'], dic['density_earth'], dic['radius_e']

        if datatype is self.type:
            return self

        elif datatype is DataType.EWH:
            ln = ln[:self.nmax + 1]
            kn = np.array([(2 * n + 1) / (1 + ln[n]) for n in range(len(ln))]) * radius_e * density_earth / (
                    3 * density_water)
            knm = np.array([[kn[l] for i in range(self.nmax + 1)] for l in range(self.nmax + 1)])
            self.Cnm2d *= knm
            self.Snm2d *= knm
            self.type = DataType.EWH
            return self

        elif datatype is DataType.density:
            ln = ln[:self.nmax + 1]
            kn = np.array([(2 * n + 1) / (1 + ln[n]) for n in range(len(ln))]) * radius_e * density_earth / 3
            knm = np.array([[kn[l] for i in range(self.nmax + 1)] for l in range(self.nmax + 1)])
            self.Cnm2d *= knm
            self.Snm2d *= knm
            self.type = DataType.EWH
            return self

        else:
            print('Failed to translate SHC type, check the input')
            return self

    def replace(self, coefficients: dict, *, c00=False, c10=False, c11=False, s11=False, c20=False, c30=False,
                info=True):
        """
        Replace some certain (low-)degree/order coefficients. Currently, it supports degree-1, c20 and c30.
        :param coefficients: dict, coefficients to replace, given by LoadL2Product.LoadLowDegree.
        :param c00: bool, True for replacing c00, False for not.
        :param c10: same as aforementioned.
        :param c11: same as aforementioned.
        :param s11: same as aforementioned.
        :param c20: same as aforementioned.
        :param c30: same as aforementioned.
        :param info: bool, True for print the replacing information if some are not replaced.
        :return:
        """

        if c00:
            assert 'c00' in coefficients.keys()
            flag = False
            for key in coefficients['c00'].keys():
                if int(self.simple_begin()) <= int(key) <= int(self.simple_end()) and coefficients['c00'][key] == \
                        coefficients['c00'][key]:
                    self.Cnm2d[0][0] = coefficients['c00'][key]
                    flag = True
                    break
            if info and not flag:
                print('c00 of {} is not replaced in this situation.'.format(self.begin_date))
        if c10:
            assert 'c10' in coefficients.keys()
            flag = False
            for key in coefficients['c10'].keys():
                if int(self.simple_begin()) <= int(key) <= int(self.simple_end()) and coefficients['c10'][key] == \
                        coefficients['c10'][key]:
                    self.Cnm2d[1][0] = coefficients['c10'][key]
                    flag = True
                    break
            if info and not flag:
                print('c10 of {} is not replaced in this situation.'.format(self.begin_date))

        if c11:
            assert 'c11' in coefficients.keys()
            flag = False
            for key in coefficients['c11'].keys():
                if int(self.simple_begin()) <= int(key) <= int(self.simple_end()) and coefficients['c11'][key] == \
                        coefficients['c11'][key]:
                    self.Cnm2d[1][1] = coefficients['c11'][key]
                    flag = True
                    break
            if info and not flag:
                print('c11 of {} is not replaced in this situation.'.format(self.begin_date))

        if s11:
            assert 's11' in coefficients.keys()
            flag = False
            for key in coefficients['s11'].keys():
                if int(self.simple_begin()) <= int(key) <= int(self.simple_end()) and coefficients['s11'][key] == \
                        coefficients['s11'][key]:
                    self.Snm2d[1][1] = coefficients['s11'][key]
                    flag = True
                    break
            if info and not flag:
                print('s11 of {} is not replaced in this situation.'.format(self.begin_date))

        if c20:
            assert 'c20' in coefficients.keys()
            flag = False
            for key in coefficients['c20'].keys():
                if int(self.simple_begin()) <= int(key) <= int(self.simple_end()) and coefficients['c20'][key] == \
                        coefficients['c20'][key]:
                    self.Cnm2d[2][0] = coefficients['c20'][key]
                    flag = True
                    break
            if info and not flag:
                print('c20 of {} is not replaced in this situation.'.format(self.begin_date))

        if c30:
            assert 'c30' in coefficients.keys()
            flag = False
            for key in coefficients['c30'].keys():
                if int(self.simple_begin()) <= int(key) <= int(self.simple_end()) and coefficients['c30'][key] == \
                        coefficients['c30'][key]:
                    self.Cnm2d[3][0] = coefficients['c30'][key]
                    flag = True
                    break
            if info and not flag:
                print('c30 of {} is not replaced in this situation.'.format(self.begin_date))

        return self

    def CS1d(self):
        """

        :return: tuple of numpy.ndarray, C,S in 1d format. e,g, ([c00, c10, c11, c20, ...], [...])
        """
        Cnm1d = GeoMathKit.CS_2dTo1d(self.Cnm2d)
        Snm1d = GeoMathKit.CS_2dTo1d(self.Snm2d)
        return Cnm1d, Snm1d

    def truncate(self, degree: int):
        """
        Keep Cnm2d, Snm2d coefficients below a smaller maximum degree/order.
        :param degree: int, the new maximum degree/order.
        :return:
        """
        assert self.nmax >= degree
        self.Cnm2d = self.Cnm2d[:degree + 1, :degree + 1]
        self.Snm2d = self.Snm2d[:degree + 1, :degree + 1]
        self.nmax = degree
        return self
