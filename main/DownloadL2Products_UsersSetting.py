from main.DownloadL2Products import downloadL2Products, downloadLowDegrees
from pysrc.Setting import *

# ==================== user's setting ==================== #

begin = 2002
end = 2022
institute = [L2instituteType.CSR, L2instituteType.JPL]
product = [L2ProductType.GAA, L2ProductType.GAB, L2ProductType.GAC, L2ProductType.GAD]
update_mode = True
local_path_GRACEL2Products = '../data/L2_SH_products'
# The recommended path is '../data/L2_SH_products'; None for not downloading

local_path_GSMLowDegrees = None
# The recommended path is '../data/LowDegreeReplace'; None for not downloading

# ================ end of user's setting ================= #

if __name__ == '__main__':
    if local_path_GRACEL2Products is not None:
        if type(institute) is L2instituteType:
            institute = [institute]

        if type(product) is L2ProductType:
            product = [product]

        for i in range(len(institute)):
            for j in range(len(product)):
                downloadL2Products(begin, end, local_path_GRACEL2Products, institute[i], product[j], update_mode)

    if local_path_GSMLowDegrees is not None:
        downloadLowDegrees(local_path_GSMLowDegrees)
