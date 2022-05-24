import copy
import datetime
import re

import numpy as np

from pysrc.SHC import SHC
from pysrc.Setting import DataType
from pysrc.TimeTool import TimeTool


class LoadL2Product:
    """
    This class is to load GRACE(-FO) L2 product into SHC type.

    usage:
        load = LoadL2Product(filepath)
        load.getSHC()
        /
        load = LoadL2Product(filepath)
        load.getCS2d()

    parameters:
        __init__():
            filepath: filename with directory of file to be loaded.

    return:
        load.getSHC():
            SHC type, for details, see SHC.py.
        load.getCS2d():
            tuple, (Cnm: np.ndarray in 2d, Snm: np.ndarray in 2d)


    """

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.filename = filepath.split('/')[-1]
        pass

    def getSHC(self):
        shc = SHC()
        shc.type = DataType.geoid

        # load date
        filename = self.filename
        pat_key = '(GSM|GAA|GAB|GAC|GAD)-2_([0-9]{7})-([0-9]{7})'
        groups = re.search(pat_key, filename, re.M).groups()  # ('GSM', '2008001', '2008031')

        begin_year_first = datetime.date(int(groups[1][:4]), 1, 1)
        begin_days = int(groups[1][4:])
        begin_date = begin_year_first + datetime.timedelta(begin_days - 1)

        end_year_first = datetime.date(int(groups[2][:4]), 1, 1)
        end_days = int(groups[2][4:])
        end_date = end_year_first + datetime.timedelta(end_days - 1)

        shc.begin_date = begin_date
        shc.end_date = end_date

        # start read file
        with open(self.filepath) as f:
            txt = f.read()
        pat_degree = r'(max_)?degree *(: *)?\d+'
        deg = int(re.search(pat_degree, txt).group().split()[-1])

        # unused days
        pat_unused_days = r'unused_days.*(\[.*\])'
        unues = re.search(pat_unused_days, txt)
        if unues is not None:
            pat_unused_days_each = r'\d{4}-\d{2}-\d{2}'
            unues_each = re.findall(pat_unused_days_each, unues.group())
            unues_days = []
            for i in range(len(unues_each)):
                ymd = unues_each[i].split('-')
                year, month, day = int(ymd[0]), int(ymd[1]), int(ymd[2])
                unues_days.append(datetime.date(year, month, day))
            shc.unused_days = unues_days
        pass

        # Cnm, Snm
        shc.nmax = deg
        shc.Cnm2d = np.zeros((deg + 1, deg + 1))
        shc.Snm2d = np.zeros((deg + 1, deg + 1))

        key1 = 'gfc'
        key2 = 'GRCOF2'
        pat_data = key1 + r'.*\d|' + key2 + r'.*\d'
        data = re.findall(pat_data, txt)
        for i in data:
            line = i.split()
            n = int(line[1])
            m = int(line[2])
            if 'GSM' in self.filepath and n <= 1:
                continue
            shc.Cnm2d[n][m] = float(line[3])
            shc.Snm2d[n][m] = float(line[4])

        return shc

    def getCS2d(self):
        # Load sphere harmonic coefficients
        with open(self.filepath) as f:
            txt = f.read()
        pat_degree = 'max_degree *\d+'
        deg = int(re.search(pat_degree, txt).group().split()[-1])

        Cnm2d = np.zeros((deg + 1, deg + 1))
        Snm2d = np.zeros((deg + 1, deg + 1))

        key1 = 'gfc'
        key2 = 'GRCOF2'
        pat_data = key1 + '.*\d|' + key2 + '.*\d'
        data = re.findall(pat_data, txt)
        for i in data:
            line = i.split()
            n = int(line[1])
            m = int(line[2])
            Cnm2d[n][m] = float(line[3])
            Snm2d[n][m] = float(line[4])

        return Cnm2d, Snm2d


class LoadLowDegree:
    """
    This class is to load auxiliary low-degree GSM product into a dict.

    usage:
        low_degrees = LoadLowDegree(filepath).coefficients()
        low_degrees.update(LoadLowDegree(filepath).coefficients())

    parameters:
        __init__():
            filepath: filename with directory of file to be loaded.

    return:
        coefficients():
            dict, e.g. {'c20': {'200801': value, ...}, 'c30':{'200801': value, ...}}
    """

    def __init__(self, filepath):
        self.filepath = filepath
        pass

    def coefficients(self):
        """

        :return: {'c20': {'200801': value}, 'c30':{'200801': value}}
        """
        if 'TN-07' in self.filepath:
            with open(self.filepath) as f:
                txt = f.read()

            C20 = {}

            pat_data = r'\s*^\d{5}.*'
            data = re.findall(pat_data, txt, re.M)
            for i in data:
                line = i.split()
                day_begin = TimeTool.mjd2yyyymmdd(float(line[0]))[:6]  # 200801
                C20[day_begin] = float(line[2])

            return {'c20': C20}

        elif 'TN-11' in self.filepath:
            with open(self.filepath) as f:
                txt = f.read()

            C20 = {}

            pat_data = r'\s*^\d{5}.*'
            data = re.findall(pat_data, txt, re.M)
            for i in data:
                line = i.split()
                day_begin = TimeTool.mjd2yyyymmdd(float(line[0]))[:6]
                C20[day_begin] = float(line[2])

            return {'c20': C20}

        elif 'TN-13' in self.filepath:
            with open(self.filepath) as f:
                txt = f.read()

            C10 = {}
            C11 = {}
            S11 = {}

            pat_data = r'^GRCOF2.*\w'
            data = re.findall(pat_data, txt, re.M)
            for i in data:
                line = i.split()
                m = int(line[2])
                day_begin = line[7][:6]
                if m == 0:
                    C10[day_begin] = float(line[3])
                elif m == 1:
                    C11[day_begin] = float(line[3])
                    S11[day_begin] = float(line[4])

            return {'c10': C10, 'c11': C11, 's11': S11}

        elif 'TN-14' in self.filepath:
            with open(self.filepath) as f:
                txt = f.read()

            C20 = {}
            C30 = {}

            pat_data = r'\s*^\d{5}.*'
            data = re.findall(pat_data, txt, re.M)
            for i in data:
                line = i.split()
                day_begin = TimeTool.mjd2yyyymmdd(int(line[0][:5]))[:6]
                C20[day_begin] = float(line[2])
                C30[day_begin] = float(line[5])

            return {'c20': C20, 'c30': C30}

        else:
            print('check low-degree file path.')
            return {}


def demo():
    from pysrc.LoveNumber import LoveNumber, LoveNumberType
    from pysrc.Setting import DataType
    shc1 = LoadL2Product(
        '../data/L2_SH_Products/RL06/GFZ/GSM/BA01/2002/GSM-2_2002095-2002120_GRAC_GFZOP_BA01_0600').getSHC()

    shc2 = copy.deepcopy(shc1)

    LN = LoveNumber('../data/Auxiliary/')
    ln = LN.getNumber(60, LoveNumberType.Wang)

    shc1.convertTypeTo(DataType.density, ln)

    C1 = shc1.Cnm2d

    shc2.convertTypeTo(DataType.EWH, ln)
    shc1.convertTypeTo(DataType.density, ln)

    C2 = shc1.Cnm2d

    print(np.max(np.abs(C1)), np.max(np.abs(C1 - C2)))

    pass


if __name__ == '__main__':
    demo()
