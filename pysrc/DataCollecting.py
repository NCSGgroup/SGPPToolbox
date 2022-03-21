import ftplib
import os

import numpy as np

from pysrc.GeoMathKit import GeoMathKit
from pysrc.Setting import L2instituteType, L2ProductType, SAT


class RetrieveL2SH:
    """
    This class is to download GRACE(-FO) Gravity spherical model(GSM) products from ftp://isdcftp.gfz-potsdam.de

    usage:
        download = RetrieveL2SH().config(LocalDir, Institute, Product, Sat)
        download.byYear(begin, end, update)

    parameters:
        config():
            LocalDir: str, local directory where downloaded files are saved.
            Institute: str, currently, 'CSR', 'GFZ' and 'JPL' can be chosen.
            Product: L2ProductType(Enum), currently, L2ProductType.GSM, .GAA, .GAB, .GAC and .GAD can be chosen.
            Sat: SAT(Enum), currently, SAT.GRACE and SAT.GRACE_FO can be chosen.
        byYear():
            begin: begin year of the data to download.
            end: end year of the data to download.
            update: bool, True for ignore the files that already exist, False for downloading all.
    """

    def __init__(self):
        self._LocalDir = None
        self._Ins = None
        self._Product = None
        self._Sat = SAT.GRACE_FO

        host = 'isdcftp.gfz-potsdam.de'
        username = 'anonymous'
        password = 'what'
        self._ftp = ftplib.FTP(host)
        self._ftp.login(username, password)

        pass

    def config(self, LocalDir: str, Institute: L2instituteType, Product: L2ProductType, Sat: SAT):
        self._LocalDir = LocalDir
        assert os.path.exists(LocalDir)

        self._Ins = Institute
        self._Product = Product
        self._Sat = Sat

        if self._Sat == SAT.GRACE_FO:
            self._ftp.cwd('/grace-fo/Level-2/' + Institute.name + '/RL06/')
        elif self._Sat == SAT.GRACE:
            self._ftp.cwd('/grace/Level-2/' + Institute.name + '/RL06/')

        return self

    def byYear(self, begin, end, update: bool = False):
        try:
            years = np.arange(int(begin), int(end) + 1)

            files = self._ftp.nlst()
            files.sort()

            path1 = os.path.join(self._LocalDir, 'RL06', self._Ins.name, self._Product.name)

            for i in range(len(files)):
                str_date = files[i].split('_')[1].split('-')
                start = str_date[0]
                end = str_date[1]
                id = files[i].split('_')[-2]  # degree/order = 60 or 90

                if (self._Product.name in files[i]) and (int(start[:4]) in years):
                    path2 = os.path.join(path1, id, start[:4])
                    if not os.path.exists(path2):
                        os.makedirs(path2)

                    if update and files[i].replace('.gz', '') in os.listdir(path2):
                        continue

                    print('Downloading:  %s' % files[i])

                    final = os.path.join(path2, files[i])

                    with open(final, "wb") as f:
                        file_handle = f.write
                        self._ftp.retrbinary('RETR %s' % files[i], file_handle, blocksize=1024)
                    GeoMathKit.un_gz(final)
                    os.remove(final)

        except Exception as r:
            print('downloading failed, {}, retrying...'.format(r))
            self.byYear(begin, end, update=True)

        pass


class RetrieveLowDegree:
    """
    This class is to download Auxiliary low-degree GSM data from ftp://isdcftp.gfz-potsdam.de

    usage:
        download = RetrieveLowDegree(LocalDir)
        download.degree_one()
        /
        download = RetrieveLowDegree(LocalDir)
        download.degree_one()

    parameters:
        __init__():
            LocalDir: str, local directory where downloaded files are saved.
    """

    def __init__(self, LocalDir: str):

        self._localDir = LocalDir
        assert os.path.exists(LocalDir)

        host = 'isdcftp.gfz-potsdam.de'
        username = 'anonymous'
        password = 'what'
        self._ftp = ftplib.FTP(host)
        self._ftp.login(username, password)

        self._ftp.cwd('/grace-fo/DOCUMENTS/TECHNICAL_NOTES/')

    def degree_one(self):
        files = self._ftp.nlst()
        files.sort()

        for file in files:

            if 'GEOC' in file:
                final = os.path.join(self._localDir, file)
                print('Downloading:  %s' % file)
                with open(final, "wb") as f:
                    file_handle = f.write
                    self._ftp.retrbinary('RETR %s' % file, file_handle, blocksize=1024)

        pass

    def degree_C20_C30(self):

        files = self._ftp.nlst()

        for file in files:

            if 'C20' in file:
                final = os.path.join(self._localDir, file)
                print('Downloading:  %s' % file)
                with open(final, "wb") as f:
                    file_handle = f.write
                    self._ftp.retrbinary('RETR %s' % file, file_handle, blocksize=1024)

        pass


