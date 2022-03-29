from main.getTimeSeries import calculate

from pysrc.Setting import *

# ==================== user's setting ==================== #

basin = Basin.Greenland
begin, end = 2005, 2015
release = L2ProductRelease.RL06
institute = L2instituteType.CSR
max_degree = GSMMaxDegree.degree60
replace_list = [LowDegree.C10, LowDegree.C11, LowDegree.S11, LowDegree.C20, LowDegree.C30]
replace_C20_GSFC = True
background = None

dec = (DecorrelatedFilterType.PnMm, (5, 7))
gauss = (GaussianFilterType.isotropic, (300,))
leakage_method = LeakageMethod.Wahr2006

gia_model = GIAModel.Caron2018
gax_model = None

save_path = '../results/test_greenland/'
save_name = 'greenland'

# ================ end of user's setting ================= #

if __name__ == '__main__':
    calculate(basin, begin, end, release, institute, max_degree, replace_list, replace_C20_GSFC, dec, gauss,
              leakage_method, gia_model, gax_model, save_path, save_name, background)
