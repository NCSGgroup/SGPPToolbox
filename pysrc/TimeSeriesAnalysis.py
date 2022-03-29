import numpy as np

from pysrc.GeoMathKit import GeoMathKit
from pysrc.Grid import Grid
from pysrc.LeastSquare import curve_fit


def fit_function_with_semiannual(x, a, b, c, d, e, f):
    """
    Linear Trend + Annual + Semiannual
    """
    return a + b * x + c * np.sin(2 * np.pi * x) + d * np.cos(2 * np.pi * x) + e * np.sin(4 * np.pi * x) + f * np.cos(
        4 * np.pi * x)


def fit_function_without_semiannual(x, a, b, c, d):
    """
    Linear Trend + Annual + Semiannual
    """
    return a + b * x + c * np.sin(2 * np.pi * x) + d * np.cos(2 * np.pi * x)


class AreaAnalysis:
    def __init__(self, area):
        self.time_points = []
        self.maps = []
        self.values = []
        self.base = None
        self.trend = None
        self.annual_amplitude = None
        self.annual_phase = None
        self.semiannual_amplitude = None
        self.semiannual_phase = None
        self.delta_trend = None
        self.delta_annual_amplitude = None
        self.delta_annual_phase = None
        self.delta_semiannual_amplitude = None
        self.delta_semiannual_phase = None
        self.area = area
        pass

    @staticmethod
    def __fc(x, a, b, c, d, e, f):
        return a + b * x + c * np.cos(x * (np.pi * 2)) + d * np.sin(x * (np.pi * 2)) + e * np.cos(
            x * (np.pi * 4)) + f * np.sin(x * (np.pi * 4))

    @staticmethod
    def __curve_fit(x, y):
        results, deltas = curve_fit(AreaAnalysis.__fc, x, y)

        return results, deltas

    def setGrids(self, series: list):
        for i in range(len(series)):
            self.time_points.append(series[i].year_frac_average())
            self.maps.append(series[i].map)
            self.values.append(GeoMathKit.gridSum(series[i].map, self.area))

        results, deltas = self.__curve_fit(self.time_points, self.values)
        self.trend = results[1]
        self.annual_amplitude = np.sqrt(results[2] ** 2 + results[3] ** 2)
        self.annual_phase = np.arctan(results[2] / results[3])
        self.semiannual_amplitude = np.sqrt(results[4] ** 2 + results[5] ** 2)
        self.annual_phase = np.arctan(results[4] / results[5])

        self.delta_trend = deltas[1]
        self.delta_annual_amplitude = np.sqrt(deltas[2] ** 2 + deltas[3] ** 2)
        self.delta_annual_phase = np.arctan(deltas[2] / deltas[3])
        self.delta_semiannual_amplitude = np.sqrt(deltas[4] ** 2 + deltas[5] ** 2)
        self.delta_annual_phase = np.arctan(deltas[4] / deltas[5])

        return self


class GridAnalysis:
    def __init__(self):
        self.time_points = []
        self.maps = []
        self.trend = None
        self.annual_amplitude = None
        self.semiannual_amplitude = None
        self.seasonal_amplitude = None
        pass

    @staticmethod
    def __fc(x, a, b, c, d, e, f, g, h):
        return a + b * x + c * np.cos(x * (np.pi * 2)) + d * np.sin(x * (np.pi * 2)) + e * np.cos(
            x * (np.pi * 4)) + f * np.sin(x * (np.pi * 4)) + g * np.cos(
            x * (np.pi * 8)) + h * np.sin(x * (np.pi * 8))

    @staticmethod
    def __curve_fit(x, y):

        a, b, c, d, e, f, g, h = curve_fit(GridAnalysis.__fc, x, y)[0]
        return a, b, c, d, e, f, g, h

    def setGrids(self, series: list):
        for i in range(len(series)):
            self.time_points.append(series[i].year_frac_average())
            self.maps.append(series[i].map)

        trend = Grid(np.zeros_like(self.maps[0]))
        annual_amplitude = Grid(np.zeros_like(self.maps[0]))
        semiannual_amplitude = Grid(np.zeros_like(self.maps[0]))
        seasonal_amplitude = Grid(np.zeros_like(self.maps[0]))

        for i in range(len(self.maps[0])):
            print('\r{}/{}'.format(i + 1, len(self.maps[0])), end='')
            for j in range(len(self.maps[0][i])):
                y = np.array([self.maps[t][i][j] for t in range(len(self.time_points))])
                a, b, c, d, e, f, g, h = self.__curve_fit(self.time_points, y)
                trend.map[i][j] = b
                annual_amplitude.map[i][j] = np.sqrt(c ** 2 + d ** 2)
                semiannual_amplitude.map[i][j] = np.sqrt(e ** 2 + f ** 2)
                seasonal_amplitude.map[i][j] = np.sqrt(g ** 2 + h ** 2)

        self.trend = trend
        self.annual_amplitude = annual_amplitude
        self.semiannual_amplitude = semiannual_amplitude
        self.seasonal_amplitude = seasonal_amplitude
        return self


class For1d:
    def __init__(self):
        self.semi_analysis = True
        self.fit_function = fit_function_with_semiannual

        self.trend = None

        self.annual_amplitude = None
        self.annual_phase = None

        self.semiannual_amplitude = None
        self.semiannual_phase = None

        self.delta_trend = None

        self.delta_annual_amplitude = None
        self.delta_annual_phase = None

        self.delta_semiannual_amplitude = None
        self.delta_semiannual_phase = None

        pass

    def semiannual_on(self, on: bool = True):
        if not on:
            self.semi_analysis = False
            self.fit_function = fit_function_without_semiannual
        else:
            self.semi_analysis = True
            self.fit_function = fit_function_with_semiannual

        return self

    def setSignals(self, times, values):
        """
        :param times: iter, year fractions, for example, [2002., 2002.083, ...]
        :param values: iter
        """
        fit_result = curve_fit(self.fit_function, times, values)
        z = fit_result[0][0]
        delta_z = np.sqrt(np.diag(fit_result[1][0]))

        self.trend = z[1]
        self.delta_trend = delta_z[1]

        self.annual_amplitude = np.sqrt(z[2] ** 2 + z[3] ** 2)
        self.annual_phase = np.degrees(np.arctan(z[3] / z[2]))  # phi = arctan(C/S)

        self.delta_annual_amplitude = np.sqrt(delta_z[2] ** 2 + delta_z[3] ** 2)
        self.delta_annual_phase = np.degrees(np.arctan(delta_z[3] / delta_z[2]))

        if self.semi_analysis:
            self.semiannual_amplitude = np.sqrt(z[4] ** 2 + z[5] ** 2)
            self.semiannual_phase = np.degrees(np.arctan(z[5] / z[4]))

            self.delta_semiannual_amplitude = np.sqrt(delta_z[4] ** 2 + delta_z[5] ** 2)
            self.delta_semiannual_phase = np.degrees(np.arctan(delta_z[5] / delta_z[4]))

        return self
