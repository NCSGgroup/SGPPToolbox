import os

import numpy as np

from main.getTimeSeries import load_GSM_by_year
from pysrc.Decorrelated import PnMm, StableWindow, SlideWindow
from pysrc.GaussianFilter import IsoGaussian, AniGaussian, Fan
from pysrc.Harmonic import Harmonic
from pysrc.LoadL2Product import LoadL2Product, LoadLowDegree
from pysrc.LoveNumber import LoveNumber, LoveNumberType
from pysrc.Plot import plot_grid
from pysrc.RefEllipsoid import EllipsoidType
from pysrc.SHC import SHC
from pysrc.Setting import *


def get_mean(begin, end, release, institute, max_degree):
    load = load_GSM_by_year(begin, end, release, institute, max_degree, toEWH=False)[0]
    return SHC.mean(load)


def show_month(release, institute, max_degree, year, month, replace, replace_C20_GSFC: bool, toEWH: bool, dec, gauss,
               pic_min, pic_max, pic_title, save_path, save_name, background: SHC = None, area=None):
    L2id = ['BA01', 'BB01'][[GSMMaxDegree.degree60, GSMMaxDegree.degree96].index(max_degree)]
    nmax = [60, 96][[GSMMaxDegree.degree60, GSMMaxDegree.degree96].index(max_degree)]

    path = '../data/L2_SH_Products/{}/{}/{}/{}/{}/'.format(release.name, institute.name, 'GSM', L2id, str(year))

    files_in_this_year = os.listdir(path)

    shc = None
    for i in range(len(files_in_this_year)):
        shc = LoadL2Product(os.path.join(path, files_in_this_year[i])).getSHC()
        if int(shc.begin_date.month) == int(month):
            break

    assert shc is not None, 'file dose not exist, please check.'

    if replace is None:
        replace = []
    c10, c11, s11 = LowDegree.C10 in replace, LowDegree.C11 in replace, LowDegree.S11 in replace
    c20, c30 = LowDegree.C20 in replace, LowDegree.C30 in replace

    lowdegrees = {}
    if replace:
        if institute not in ['CSR', 'JPL', 'GFZ']:
            deg1 = LoadLowDegree('../data/LowDegreeReplace/TN-13_GEOC_CSR_RL06.txt').coefficients()
        else:
            deg1 = LoadLowDegree('../data/LowDegreeReplace/TN-13_GEOC_{}_RL06.txt'.format(institute)).coefficients()
        c20c30 = LoadLowDegree('../data/LowDegreeReplace/TN-14_C30_C20_SLR_GSFC.txt').coefficients()
        lowdegrees = {}
        lowdegrees.update(deg1)
        lowdegrees.update(c20c30)

    if not replace_C20_GSFC:
        c20_dict = LoadLowDegree('../data/LowDegreeReplace/TN-11_C20_SLR_RL06.txt').coefficients()
        lowdegrees.update(c20_dict)

    shc.replace(lowdegrees, c10=c10, c11=c11, s11=s11, c20=c20, c30=c30)

    if background is not None:
        shc -= background

    LN = LoveNumber('../data/Auxiliary/')
    ln = LN.getNumber(nmax, LoveNumberType.Wang)

    if toEWH:
        shc.convertTypeTo(DataType.EWH, ln)

    # de-correlated filter
    print('filtering...')
    decFilter = None
    if dec is not None:
        decFilterType = dec[0]
        dec_paras = dec[1]
        if decFilterType == DecorrelatedFilterType.PnMm:
            decFilter = PnMm(*dec_paras)

        elif decFilterType == DecorrelatedFilterType.StableWindow:
            decFilter = StableWindow(*dec_paras)

        elif decFilterType == DecorrelatedFilterType.SlideWindow:
            decFilter = SlideWindow(*dec_paras)

        if decFilter is not None:
            shcs = decFilter.ApplyTo(shc)

    # Gaussian(-like) low-pass filter
    if gauss is not None:
        GsFilterType = gauss[0]
        Gs_paras = gauss[1]

        if GsFilterType == GaussianFilterType.isotropic:
            GaussianFilter = IsoGaussian(*Gs_paras)

        elif GsFilterType == GaussianFilterType.anisoropic:
            GaussianFilter = AniGaussian(*Gs_paras)

        elif GsFilterType == GaussianFilterType.fan:
            GaussianFilter = Fan(*Gs_paras)

        else:
            GaussianFilter = None

    else:
        GaussianFilter = None

    if GaussianFilter is not None:
        shc_filtered = GaussianFilter.ApplyTo(shc)
    else:
        shc_filtered = shc

    '''harmonic'''
    print('harmonic synthesis...')
    lat = np.arange(-90, 90, 0.5)
    lon = np.arange(-180, 180, 0.5)
    Har = Harmonic(Parallel=-1).setEllipsoid(ell=EllipsoidType.gif48)

    grid = Har.synthesis_for_SHC(shc_filtered, lat, lon)

    if area is None:
        area = [-180, 180, -90, 90]

    if save_path is None:
        save = None
    else:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        if not save_name.endswith('.png'):
            save_name += '.png'
        save = os.path.join(save_path, save_name)
    plot_grid([grid, pic_title, pic_min, pic_max], area=area, save=save)
