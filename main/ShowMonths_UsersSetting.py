from main.ShowMonths import show_month, get_mean
from pysrc.Setting import *

# ==================== user's setting ==================== #

# shc and plot setting
year, month = 2005, 10
release = L2ProductRelease.RL06
institute = L2instituteType.CSR
max_degree = GSMMaxDegree.degree60
replace_list = [LowDegree.C10, LowDegree.C11, LowDegree.S11]
replace_C20_GSFC = True
toEWH = True
dec = (DecorrelationFilterType.PnMm, (4, 6))
gauss = (GaussianFilterType.isotropic, (300,))

area = [-180, 180, -90, 90]  # [lon_start, lon_end, lat_start, lat_end], None for global
pic_maxvalue = 0.3
pic_minvalue = -0.3
pic_title = 'global P4M6 + Gaussian300km'

# background: average GSM of a time range
begin = 2002
end = 2022

# save
save_path = '../results/test_plot/'
save_name = 'global_p4m6_gs300.png'

# ================ end of user's setting ================= #

if __name__ == '__main__':
    bc = get_mean(2002, 2010, L2ProductRelease.RL06, L2instituteType.CSR, GSMMaxDegree.degree60)

    show_month(L2ProductRelease.RL06, L2instituteType.CSR, GSMMaxDegree.degree60, year, month, replace_list,
               replace_C20_GSFC, toEWH, dec, gauss, pic_minvalue, pic_maxvalue, pic_title, save_path, save_name, bc)
