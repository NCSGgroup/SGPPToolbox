import copy

import numpy as np

from pysrc.Setting import DataType
from pysrc.TimeTool import TimeTool


class Grid:
    def __init__(self, usemap: np.ndarray = None, datatype: DataType = None):
        self.grid_space = None
        self.lon, self.lat = None, None
        self.type = datatype
        self.begin_date, self.end_date = None, None  # date.date
        self.map = None

        if usemap is not None:
            self.grid_space = 180 / np.shape(usemap)[0]
            self.lat = np.arange(-90, 90, self.grid_space)
            self.lon = np.arange(-180, 180, self.grid_space)
            self.map = usemap

    def average_date(self):
        return self.begin_date + (self.end_date - self.begin_date) / 2

    def year_frac_average(self):
        return TimeTool.year_frac(self.average_date())

    def __mul__(self, other):
        new = copy.deepcopy(self)
        new.map = self.map * other
        return new

    def __sub__(self, other):
        assert type(other) in [Grid, np.ndarray]

        if type(other) is Grid:
            other = other.map

        new = copy.deepcopy(self)
        new.map = self.map - other
        return new
