import numpy as np

from pysrc.GeoMathKit import GeoMathKit
from pysrc.Harmonic import Harmonic
from pysrc.RefEllipsoid import EllipsoidType


def get_leakage_time_series(SHCs: list, area, matrix):
    """
    Gaussian 300km
    :param grids:
    :param area:
    :return:
    """

    year_fracs = []
    signals = []
    area_bar = 1 - area

    lat = np.arange(-90, 90, 0.5)
    lon = np.arange(-180, 180, 0.5)
    Har = Harmonic(Parallel=-1).setEllipsoid(ell=EllipsoidType.gif48)
    Pnm = GeoMathKit.getPnm(lat, 60, 1)

    grids = []
    for i in range(len(SHCs)):
        shc = SHCs[i]
        print('\r', shc.begin_date, end='', sep='')
        grid = Har.synthesis_for_SHC(shc, lat, lon)
        grids.append(grid)

    for i in range(len(grids)):
        print('\r{}'.format(grids[i].begin_date), end='')
        year_fracs.append(grids[i].year_frac_average())
        this_map_bar = grids[i].map * area_bar

        Cnm, Snm = Har.analysis(60, this_map_bar, lat, lon, Pnm)
        Cnm *= matrix
        Snm *= matrix
        leakage_map = Har.synthesis(Cnm, Snm, 60, lat, lon)

        signals.append(GeoMathKit.gridSum(leakage_map, area))
    print()

    return np.array(year_fracs), np.array(signals)


class Leakage:
    def __init__(self, SHCs: list):
        self.SHCs = SHCs

    def setArea(self, area):
        self.area = area
