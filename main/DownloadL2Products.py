import os

from pysrc.DataCollecting import RetrieveL2SH, RetrieveLowDegree
from pysrc.Setting import *


def downloadL2Products(begin, end, localDir, institute, product, update_mode=False):
    if not os.path.exists(localDir):
        os.makedirs(localDir)
    download = RetrieveL2SH().config(LocalDir=localDir, Institute=institute, Product=product, Sat=SAT.GRACE)
    download.byYear(begin, end, update_mode)

    download = RetrieveL2SH().config(LocalDir=localDir, Institute=institute, Product=product, Sat=SAT.GRACE_FO)
    download.byYear(begin, end, update_mode)


def downloadLowDegrees(localDir):
    if not os.path.exists(localDir):
        os.makedirs(localDir)
    download = RetrieveLowDegree(localDir)
    download.degree_one()
    download.degree_C20_C30()
