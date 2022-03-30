import os.path

from pysrc.LoadL2Product import LoadL2Product
from pysrc.Setting import L2instituteType
from pysrc.TimeTool import TimeTool


class LoadGMAM:
    def __init__(self, institute: L2instituteType):

        unfound_flag = False
        if os.path.exists('../data/L2_SH_products/RL06/{}/GAA/BC01/'.format(institute.name)):
            self.gaa_path = '../data/L2_SH_products/RL06/{}/GAA/BC01/'.format(institute.name)
        else:
            unfound_flag = True

        if unfound_flag:
            if os.path.exists('../data/L2_SH_products/RL06/CSR/GAA/BC01/'):
                self.gaa_path = '../data/L2_SH_products/RL06/CSR/GAA/BC01/'
            elif os.path.exists('../data/L2_SH_products/RL06/GFZ/GAA/BC01/'):
                self.gaa_path = '../data/L2_SH_products/RL06/GFZ/GAA/BC01/'
            else:
                self.gaa_path = '../data/L2_SH_products/RL06/JPL/GAA/BC01/'

    def loadby(self, begin, end):
        coefficients = {}
        for year in range(begin, end + 1):
            all_files = os.listdir(os.path.join(self.gaa_path, str(year)))
            for i in range(len(all_files)):
                shc = LoadL2Product(os.path.join(self.gaa_path, str(year), all_files[i])).getSHC()
                date_yyyymm = TimeTool.yyyymm(shc.average_date())
                coefficients[date_yyyymm] = shc.Cnm2d[0][0]

        return {'c00': coefficients}


def demo():
    coefficients = LoadGMAM(L2instituteType.JPL).loadby(2002, 2020)
    print(coefficients)


if __name__ == '__main__':
    demo()
