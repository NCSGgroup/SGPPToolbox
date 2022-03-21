from main.getTimeSeries import calculate

from pysrc.Setting import *

# ==================== user's setting ==================== #
basin = Basin.Yangtze
begin, end = 2002, 2017
release = L2ProductRelease.RL06
institute = L2instituteType.GFZ
max_degree = L2MaxDegree.degree60
replace_list = [LowDegree.C10, LowDegree.C11, LowDegree.S11, LowDegree.C20]
replace_C20_GSFC = True
dec = (DecorrelatedFilterType.PnMm, (3, 5))
gauss = (GaussianFilterType.anisoropic, (300, 500, 15))
leakage_method = LeakageMethod.BufferZone
gia_model = GIAModel.ICE6G_D

save_path = '../results/test_yangtze/'
save_name = 'yangtze'
background = None
# ================ end of user's setting ================= #

if __name__ == '__main__':
    calculate(basin, begin, end, release, institute, max_degree, replace_list, replace_C20_GSFC, dec, gauss,
              leakage_method, gia_model, save_path, save_name, background)
