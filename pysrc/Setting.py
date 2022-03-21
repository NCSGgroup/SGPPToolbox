from enum import Enum


class L2ProductType(Enum):
    GAA = 0
    GAB = 1
    GAC = 2
    GAD = 3
    GSM = 4


class L2ProductRelease(Enum):
    RL06 = 1


class L2instituteType(Enum):
    GFZ = 0
    CSR = 1
    JPL = 2


class L2MaxDegree(Enum):
    degree60 = 1
    degree96 = 2


class SAT(Enum):
    GRACE = 0
    GRACE_FO = 1


class LowDegree(Enum):
    C10 = 1
    C11 = 2
    S11 = 3
    C20 = 4
    C30 = 5


class DataType(Enum):
    geoid = 1
    EWH = 2
    density = 3


class EarthModel(Enum):
    general = 1


class GaussianFilterType(Enum):
    isotropic = 1
    anisoropic = 2
    fan = 3


class LeakageMethod(Enum):
    Wahr2006 = 1
    ForwardModeling = 2
    BufferZone = 3
    GainFactor = 4


class DecorrelatedFilterType(Enum):
    SlideWindow = 1
    StableWindow = 2
    PnMm = 3


class ProjectionType(Enum):
    PlateCarree = 1
    Mollweide = 2
    Orthographic = 3
    NorthPolarStereo = 4
    SouthPolarStereo = 5


class HydrologicModel(Enum):
    NOAH21 = 1


class GIAModel(Enum):
    ICE5G_A = 1
    ICE6G_C = 2
    ICE6G_D = 3
    Caron2018 = 4
    Caron2019 = 5


class Basin(Enum):
    Amazon = 1
    Amur = 2
    Antarctica = 3
    Aral = 4
    Brahmaputra = 5
    Caspian = 6
    Colorado = 7
    Congo = 8
    Danube = 9
    Dnieper = 10
    Euphrates = 11
    Eyre = 12
    Ganges = 13
    Greenland = 14
    Greenland_North = 15
    Greenland_South = 16
    Indus = 17
    Lena = 18
    Mackenzie = 19
    Mekong = 20
    Mississippi = 21
    Murray = 22
    Nelson = 23
    Niger = 24
    Nile = 25
    Ob = 26
    Okavango = 27
    Orange = 28
    Orinoco = 29
    Parana = 30
    Sahara = 31
    St_Lawrence = 32
    Tocantins = 33
    Yangtze = 34
    Yellow = 35
    Yenisey = 36
    Yukon = 37
    Zambeze = 38
    Ocean = 39
