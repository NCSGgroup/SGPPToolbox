from main.DownloadL2Products import downloadL2Products, downloadLowDegrees
from pysrc.Setting import *

# ==================== user's setting ==================== #

begin = 2002
end = 2022
local_path_GRACEL2Products = '../data/L2_SH_products'
institute = L2instituteType.CSR
product = L2ProductType.GSM

local_path_GSMLowDegrees = '../data/LowDegreeReplace'

# ================ end of user's setting ================= #

if __name__ == '__main__':

    downloadL2Products(begin, end, local_path_GRACEL2Products, institute, product)

    downloadLowDegrees(local_path_GSMLowDegrees)
