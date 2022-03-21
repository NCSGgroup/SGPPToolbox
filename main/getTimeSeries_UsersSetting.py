from main.getTimeSeries import calculate

from pysrc.Setting import *

if __name__ == '__main__':
    basin = Basin.Ocean
    begin, end = 2002, 2017
    release = L2ProductRelease.RL06
    toEWH = True
    institute = L2instituteType.GFZ
    max_degree = L2MaxDegree.degree60
    replace_list = [LowDegree.C10, LowDegree.C11, LowDegree.S11, LowDegree.C20]
    replace_C20_GSFC = True
    dec = (DecorrelatedFilterType.PnMm, (3, 5))
    gauss = (GaussianFilterType.isotropic, (300,))
    leakage_method = LeakageMethod.BufferZone
    gia_model = GIAModel.ICE6G_D

    save_path = '../results/20220321/test_ocean/'
    save_name = 'ocean'
    background = None

    calculate(basin, begin, end, release, toEWH, institute, max_degree, replace_list,
              replace_C20_GSFC, dec, gauss, leakage_method, gia_model, save_path, save_name, background)
