import numpy as np
from netCDF4 import Dataset

from pysrc.GeoMathKit import GeoMathKit


class LoadNOAH21:
    def __init__(self):
        self.nc = None
        self.keys = None

    def setFile(self, file):
        self.nc = Dataset(file)
        self.keys = self.nc.variables.keys()
        return self

    def get2dData(self, key):
        assert key in self.keys, 'no such key word'
        data = np.array(self.nc.variables[key])[0]
        xyz = []
        for i in range(len(data)):
            for j in range(len((data[i]))):
                lat = int(i - 90) + 30
                lon = j - 180
                xyz.append([lon, lat, data[i][j]])
        xyz = np.array(xyz)
        map = GeoMathKit.xyz2grd(xyz)
        return map


def getTWS(file):
    """
    The GLDAS/Noah soil moisture (SM), snow water equivalent (SWE), and plant canopy water storage (PCSW) are jointly
    used to calculate the TWS variations. doi: 10.1155/2019/3874742
    :param file: path + filename.nc of NOAH
    :return: 1*1 degree TWS map [m]
    """
    nc = LoadNOAH21().setFile(file)
    sm0_10 = nc.get2dData('SoilMoi0_10cm_inst')
    sm10_40 = nc.get2dData('SoilMoi10_40cm_inst')
    sm40_100 = nc.get2dData('SoilMoi40_100cm_inst')
    sm100_200 = nc.get2dData('SoilMoi100_200cm_inst')
    cano = nc.get2dData('CanopInt_inst')
    swe = nc.get2dData('SWE_inst')
    return (sm0_10 + sm10_40 + sm40_100 + sm100_200 + cano + swe) / 1000


if __name__ == '__main__':
    tws = getTWS('../data/Noah2.1/GLDAS_NOAH10_M.A200204.021.nc4.SUB.nc4')
    print(tws)
