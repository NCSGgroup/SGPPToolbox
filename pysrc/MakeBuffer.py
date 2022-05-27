import numpy as np

from pysrc.GaussianFilter import IsoGaussian, AniGaussian, Fan
from pysrc.GeoMathKit import GeoMathKit
from pysrc.Grid import Grid
from pysrc.Harmonic import Harmonic
from pysrc.Plot import plot_grid
from pysrc.Setting import GaussianFilterType


def make_buffer(area, gaussian_filter, weight=0.1):
    sp = 180 / np.shape(area)[0]
    nmax = 60

    lat = np.arange(-90, 90, sp)
    lon = np.arange(-180, 180, sp)
    colat_rad, lon_rad = GeoMathKit.getCoLatLoninRad(lat, lon)
    PnmMatrix = GeoMathKit.getPnmMatrix(colat_rad, nmax, 0)
    Har = Harmonic(lat, lon, PnmMatrix, nmax)
    area_complement = 1 - area

    area_complement_C = Har.analysis([area_complement])[0][0]
    area_complement_S = Har.analysis([area_complement])[1][0]

    gaussian_filter_type = gaussian_filter[0]
    gaussian_filter_paras = gaussian_filter[1]

    if gaussian_filter_type is GaussianFilterType.isotropic:
        gsfilter = IsoGaussian(*gaussian_filter_paras)

    elif gaussian_filter_type is GaussianFilterType.anisoropic:
        gsfilter = AniGaussian(*gaussian_filter_paras)

    else:
        gsfilter = Fan(*gaussian_filter_paras)

    area_complement_bar_shc = gsfilter.ApplyTo([area_complement_C, area_complement_S])
    area_complement_bar_grid = Har.synthesis([area_complement_bar_shc[0]], [area_complement_bar_shc[1]])[0]

    area_buffer = area_complement_bar_grid < weight
    return area_buffer * area


if __name__ == '__main__':
    ocean = np.load('../data/grids/ocean_maskGrid.dat(360,720).npy')

    buffered_ocean_300 = make_buffer(ocean, (GaussianFilterType.isotropic, (300,)))
    plot_grid([Grid(buffered_ocean_300), 'gaussian 300km', None, None])

    buffered_ocean_500 = make_buffer(ocean, (GaussianFilterType.isotropic, (500,)))
    plot_grid([Grid(buffered_ocean_500), 'gaussian 500km', None, None])
