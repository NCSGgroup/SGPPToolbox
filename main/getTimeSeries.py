import os

import matplotlib.pyplot as plt
import numpy as np

from pysrc.Decorrelation import PnMm, VariableWindow, StableWindow
from pysrc.ForwardModeling import ForwardModeling
from pysrc.GainFactor import SingleFactor
from pysrc.GaussianFilter import IsoGaussian, AniGaussian, Fan
from pysrc.GeoMathKit import GeoMathKit
from pysrc.Grid import Grid
from pysrc.Harmonic import Harmonic
from pysrc.Leakage_Wahr2006 import get_leakage_time_series
from pysrc.LoadL2Product import LoadL2Product, LoadLowDegree
from pysrc.LoadNoah import getTWS
from pysrc.LoveNumber import LoveNumberType, LoveNumber
from pysrc.RefEllipsoid import EllipsoidType
from pysrc.SHC import SHC
from pysrc.Setting import *
from pysrc.TimeSeriesAnalysis import For1d as TS1d


def load_GSM_by_year(begin, end, release: L2ProductRelease, institute: L2instituteType, max_degree: GSMMaxDegree,
                     toEWH=True, replace=None, GSFC=True, replace_info=True, background=None):
    """
    A convenient function to load GRACE(-FO) GSM products into SHC.
    :param begin: int, begin year.
    :param end: int, end year.
    :param release: L2ProductRelease, GSM product release, currently, only 'RL06' can be chosen.
    :param institute: L2instituteType, currently, 'CSR', 'GFZ' and 'JPL' can be chosen.
    :param max_degree: L2MaxDegree
    :param toEWH: bool, True for translating GSM data (from goied) into EWH.
    :param replace: list, which coefficients to be place, e.g. ['s11', 'c10', 'c11', 'c20', 'c30']
    :param GSFC: bool, True for choosing c20 and c30 from GSFC, False for that from CSR (if replace)
    :param replace_info: bool, True for printing the message of replacing information, False for not.
    :param background: None or SHC
    :return:
    """

    if replace is None:
        replace = []
    c10, c11, s11 = LowDegree.C10 in replace, LowDegree.C11 in replace, LowDegree.S11 in replace
    c20, c30 = LowDegree.C20 in replace, LowDegree.C30 in replace

    all_files = []
    rl = release.name
    ins = institute.name
    L2id = ['BA01', 'BB01'][[GSMMaxDegree.degree60, GSMMaxDegree.degree96].index(max_degree)]
    nmax = [60, 96][[GSMMaxDegree.degree60, GSMMaxDegree.degree96].index(max_degree)]
    for year in range(begin, end + 1):
        path = '../data/L2_SH_Products/{}/{}/{}/{}/{}/'.format(rl, ins, 'GSM', L2id, str(year))
        files_in_this_year = os.listdir(path)
        all_files += [path + files_in_this_year[i] for i in range(len(files_in_this_year))]
    all_files.sort()

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

    if not GSFC:
        c20_dict = LoadLowDegree('../data/LowDegreeReplace/TN-11_C20_SLR_RL06.txt').coefficients()
        lowdegrees.update(c20_dict)

    LN = LoveNumber('../data/Auxiliary/')
    ln = LN.getNumber(nmax, LoveNumberType.Wang)

    shc_list = []
    for i in range(len(all_files)):
        load = LoadL2Product(all_files[i])
        shc = load.getSHC()
        shc.truncate(nmax)
        shc.replace(lowdegrees, c10=c10, c11=c11, s11=s11, c20=c20, c30=c30, info=replace_info)

        if toEWH:
            shc.convertTypeTo(DataType.EWH, ln)

        shc_list.append(shc)

    if background is None:
        background = SHC.mean(shc_list)

    return shc_list, SHC.delta(shc_list, background)


def load_GAX_by_year(begin, end, model_type: L2ProductType, release: L2ProductRelease, institute: L2instituteType,
                     max_degree: GSMMaxDegree, toEWH=True, background=None):
    all_files = []
    rl = release.name
    ins = institute.name
    nmax = [60, 96][[GSMMaxDegree.degree60, GSMMaxDegree.degree96].index(max_degree)]
    L2id = 'BC01'
    for year in range(begin, end + 1):
        path = '../data/L2_SH_Products/{}/{}/{}/{}/{}/'.format(rl, ins, model_type.name, L2id, str(year))
        files_in_this_year = os.listdir(path)
        all_files += [path + files_in_this_year[i] for i in range(len(files_in_this_year))]
    all_files.sort()

    LN = LoveNumber('../data/Auxiliary/')
    ln = LN.getNumber(nmax, LoveNumberType.Wang)

    shc_list = []
    for i in range(len(all_files)):
        load = LoadL2Product(all_files[i])
        shc = load.getSHC()
        shc.truncate(nmax)

        if toEWH:
            shc.convertTypeTo(DataType.EWH, ln)

        shc_list.append(shc)

    if background is None:
        background = SHC.mean(shc_list)

    return shc_list, SHC.delta(shc_list, background)


def calculate(basin: Basin, begin, end, release, institute, max_degree, replace_list, replace_C20_GSFC: bool,
              dec, gauss, leakage_method, gia_model, gax_model, save_path, save_name, background=None):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    area = np.load('../data/grids/{}_maskGrid.dat(360,720).npy'.format(basin.name))

    # load & replace low-degree
    print('loading files...')

    shcs = load_GSM_by_year(begin, end, release=release, institute=institute, max_degree=max_degree, toEWH=True,
                            replace=replace_list, GSFC=replace_C20_GSFC, replace_info=False, background=background)[1]
    if gax_model is not None:
        shcs_gax = load_GAX_by_year(begin, end, model_type=gax_model, release=release, institute=institute,
                                    max_degree=max_degree, toEWH=True, background=None)[0]
        for i in range(len(shcs)):
            assert shcs[i].begin_date == shcs_gax[i].begin_date, 'GSM or GAX files may not be complete, please check.'
            shcs[i] -= shcs_gax[i]

    dates = np.array([shcs[i].average_date() for i in range(len(shcs))])
    year_fracs = np.array([shcs[i].year_frac_average() for i in range(len(shcs))])
    ewhs = []

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

        elif decFilterType == DecorrelatedFilterType.VariableWindow:
            decFilter = VariableWindow(*dec_paras)

        if decFilter is not None:
            shcs = decFilter.ApplyTo(*shcs)

    # Gaussian(-like) low-pass filter
    GsFilterType = gauss[0]
    Gs_paras = gauss[1]
    GaussianFilter = None
    if GsFilterType == GaussianFilterType.isotropic:
        GaussianFilter = IsoGaussian(*Gs_paras)

    elif GsFilterType == GaussianFilterType.anisoropic:
        GaussianFilter = AniGaussian(*Gs_paras)

    elif GsFilterType == GaussianFilterType.fan:
        GaussianFilter = Fan(*Gs_paras)

    if GaussianFilter is not None:
        shcs_filtered = GaussianFilter.ApplyTo(*shcs)
    else:
        shcs_filtered = shcs

    '''harmonic'''
    print('harmonic synthesis...')
    lat = np.arange(-90, 90, 0.5)
    lon = np.arange(-180, 180, 0.5)
    Har = Harmonic(Parallel=-1).setEllipsoid(ell=EllipsoidType.gif48)
    Pnm = GeoMathKit.getPnm(lat, 60, 1)

    grids = []
    for i in range(len(shcs_filtered)):
        shc = shcs_filtered[i]
        print('\r', shc.begin_date, end='', sep='')
        grid = Har.synthesis_for_SHC(shc, lat, lon)
        grids.append(grid)
    print('')

    # leakage
    if leakage_method == LeakageMethod.ForwardModeling:
        print('forward modeling...')
        area_complement = 1 - area
        ac_area = GeoMathKit.getAcreage(area)
        ac_area_complement = GeoMathKit.getAcreage(area_complement)

        for i in range(len(grids)):
            print(grids[i].begin_date)
            fm = ForwardModeling().setModel(grids[i]).setLoopTime(50).setFilter(GsFilterType, *Gs_paras).setArea(
                area_complement).setNmax(60)
            fm.run()
            ewhs.append(-GeoMathKit.gridSum(fm.model_tru, area_complement) * ac_area_complement / ac_area)

        ewhs = np.array(ewhs)
        np.save(os.path.join(save_path, save_name).replace('\\', '/'), np.array([year_fracs, ewhs]))

    elif leakage_method == LeakageMethod.BufferZone:
        assert basin == Basin.Ocean, 'buffer zone method is currently available to ocean.'
        ocean_buffered = np.load('../data/grids/ocean_300km-buffer(360,720).npy')
        print('buffer...')
        for i in range(len(grids)):
            ewhs.append(GeoMathKit.gridSum(grids[i], ocean_buffered))

        ewhs = np.array(ewhs)
        np.save(os.path.join(save_path, save_name).replace('\\', '/'), np.array([year_fracs, ewhs]))

    elif leakage_method == LeakageMethod.Wahr2006:
        print('leakage: wahr\'s method')
        for i in range(len(grids)):
            ewhs.append(GeoMathKit.gridSum(grids[i], area))
        ewhs = np.array(ewhs)

        area_filtered_shc = Har.analysis_for_Grid(60, Grid(area), Pnm) * GaussianFilter.getFilter(60)
        area_filtered = Har.synthesis_for_SHC(area_filtered_shc, lat, lon)
        k = GeoMathKit.gridSum(area, area) / GeoMathKit.gridSum(area_filtered, area)

        ewh_leakages = get_leakage_time_series(shcs, area, GaussianFilter.getFilter(60))[1]

        ewhs = k * (ewhs - ewh_leakages)

    elif leakage_method == LeakageMethod.GainFactor:
        print('leakage: gain factor')
        assert basin != Basin.Ocean, 'gain factor method is not currently available to ocean.'

        noah_path = '../data/Noah2.1/'
        noah_files = os.listdir(noah_path)
        noah_files.sort()

        tws_ts = []
        for i in range(len(noah_files)):
            this_year_month = noah_files[i].split('_')[2].split('.')[1][1:7]
            if int(this_year_month[:4]) not in range(begin, end + 1):
                continue
            print('\rget tws...{}'.format(this_year_month), end='')
            tws_ts.append(getTWS(os.path.join(noah_path, noah_files[i])))
        nmax = 60

        sf = SingleFactor().setTimeSeries(tws_ts).setNmax(nmax)
        sf.setDecFilter(dec)
        sf.setGaussianFilter(gauss)
        sf.setArea(area)

        sf.run()

        ewhs = []
        for i in range(len(grids)):
            ewhs.append(GeoMathKit.gridSum(grids[i], area))
        ewhs = np.array(ewhs)

        ewhs *= sf.k

    # GIA
    if gia_model == GIAModel.ICE5G_A:
        gia = np.load('../data/GIAModel/ice5ga_05grids.npy')

    elif gia_model == GIAModel.ICE6G_C:
        gia = np.load('../data/GIAModel/ice6gc_05grids.npy')

    elif gia_model == GIAModel.ICE6G_D:
        gia = np.load('../data/GIAModel/ice6gd_05grids.npy')

    elif gia_model == GIAModel.Caron2018:
        gia = np.load('../data/GIAModel/Caron18_05grids.npy')

    elif gia_model == GIAModel.Caron2019:
        gia = np.load('../data/GIAModel/Caron_Ivins19_05grids.npy')

    else:
        gia = np.zeros_like(area)

    gia_trend = GeoMathKit.gridSum(gia, area)

    ewhs -= gia_trend * (year_fracs - year_fracs[0]) / 1000

    # get trend et al.
    ts_analisis = TS1d().setSignals(year_fracs, ewhs)

    # save files

    save_name_npy = save_name + '.npy'
    save_name_txt = save_name + '.txt'
    save_name_png = save_name + '.png'
    np.save(os.path.join(save_path, save_name_npy).replace('\\', '/'), np.array([year_fracs, ewhs]))

    with open(os.path.join(save_path, save_name_txt).replace('\\', '/'), 'w+') as f:
        f.write('basin: {}\n'.format(basin.name))

        f.write('begin year: {}\nend year: {}\n'.format(begin, end))

        f.write('institute: {}\n'.format(institute))

        f.write('replace: ')
        for i in range(len(replace_list)):
            f.write('{} '.format(replace_list[i]))
        f.write('\n')

        if dec is not None:
            f.write('de-correlated filter: {}\t{}\n'.format(dec[0].name, dec[1]))

        f.write('Gaussian(-like) filter: {}\t{}\n'.format(GsFilterType.name, Gs_paras))

        f.write('leakage: {}\n'.format(leakage_method))

        if gia_model is not None:
            f.write('GIA: {}\n'.format(gia_model.name))
        else:
            f.write('No GIA model applied.')

        f.write(
            '\ntrend = {} mm/year \nsigma(trend) = {} mm/year \nannual amplitude = {} mm \nsigma(annual amplitude) =  {} mm \nsemi-annual amplitude = {} mm \nsigma(semi-annual amplitude) = {} mm\n'.format(
                ts_analisis.trend * 1000, ts_analisis.delta_trend * 1000,
                ts_analisis.annual_amplitude * 1000, ts_analisis.delta_annual_amplitude * 1000,
                ts_analisis.semiannual_amplitude * 1000, ts_analisis.delta_semiannual_amplitude * 1000))
        f.write('=' * 20 + 'END OF HEAD' + '=' * 20 + '\n')
        for i in range(len(year_fracs)):
            f.write('{}\t{}\n'.format(dates[i], ewhs[i]))

    plt.plot(year_fracs, ewhs)
    plt.savefig(os.path.join(save_path, save_name_png), dpi=450)
    plt.close()


if __name__ == '__main__':
    a, b, = load_GAX_by_year(2018, 2019, L2ProductType.GAA, L2ProductRelease.RL06, L2instituteType.GFZ,
                             max_degree=GSMMaxDegree.degree60, toEWH=True, background=None)
    pass
