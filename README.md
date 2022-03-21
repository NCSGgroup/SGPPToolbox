# Table of Contents # 

***

- Introduction
- Usage
- Maintainers

***

# Introduction #

Satellite Gravity Post-Processing(SGPP) Tool is used for GRACE(-FO) level-2 GSM products. From auto-downloading data to
the final time-series results, SGPP provides complete post-process to deal with every step. In each step of
post-process, SGPP provides optional methods and arguments for users.

# Usage #

## DataCollection.py ##
This module can be used to download GRACE(-FO) GSM products and auxiliary files.

This module currently includes following classes:

1. `RetrieveL2SH`
    
To download GRACE(-FO) GSM and GAX products from ftp://isdcftp.gfz-potsdam.de.

~~~python
download = RetrieveL2SH().config(LocalDir='../data/L2_SH_products', Institute=L2instituteType.CSR,
                                     Product=L2ProductType.GAA, Sat=SAT.GRACE)
download.byYear(2002, 2022)
~~~

2. `RetrieveLowDegree`

To download low-degree GSM data products(TN-XX files) from ftp://isdcftp.gfz-potsdam.de.

~~~python
download = RetrieveLowDegree(LocalDir='../data/LowDegreeReplace')
download.degree_one()
download.degree_C20_C30()
~~~

## LoadL2Product.py ##

This module is to load local GRACE(-FO) L2 products into a certain format for further use.

This module currently includes following classes:

1. `LoadL2Product`

To read local GSM or GAX files, and then give results in SHC format.

For more information of SHC, please see SHC.py.
 ~~~python
load = LoadL2Product('../data/L2_SH_products/RL06/CSR/GSM/BA01/2002/xxx')
shc = load.getSHC()
cs2d = load.getCS2d()   # .getCS2d() would return a tuple, (C_2d: np.ndarray, S_2d: np.ndarray)
~~~

2. `LoadLowDegree`

To read TN-XX files about low-degree GSM coefficients, and then give a dict of the results.

~~~python
load = LoadLowDegree('../data/LowDegreeReplace/xxx.txt')
low_deg = load.coefficients()   # {'c20': {'200801': value,...}, 'c30':{'200801': value,...},...}

shc.replace(low_deg, c10=True, c11=True, s11=True)
~~~

## Decorrelated.py ##

This module is to weaken the correlation error of the GRACE GSM models(i.e., even or odd degrees coefficients of the 
same order) using the method of deducting the linear fitting results. Currently, PnMm methods, stable-window 
de-correlation filter and slide-window de-correlation filter are available.

This module currently includes following classes:

1. `PnMm`

Starting with order m, to fit a polynomial between all even or odd degrees coefficients of the same order, then 
subtract the fitting results from the given coefficients.

~~~python
pnmm = PnMm(3, 5)   # the first parameter is the degree of the term in the polynomial that has the highest degree, the second one is starting order m.
shc_new = pnmm.ApplyTo(shc_initial)
~~~

2. `StableWindow`

Similar with `PnMm`, but for each coefficient, it would only contain several(a certain number) contiguous even or odd
degrees ones of the same order to fit. 

~~~python
stable_window_dec = StableWindow(3, 5, 5)   # the first two parameters are the same with `PnMm`, the third one is to define the length of the window(i.e., how many the fitting coefficients are).
shc_new = stable_window_dec.ApplyTo(shc_initial)
~~~

3. `SlideWindow`

Similar with `StableWindow`, but for different degree, the length of window would vary as the empirical formula with 
parameters A and K.

~~~python
slide_window_dec = SlideWindow(3, 5, 5, 30, 10)   # the first three parameters are the same with `StableWindow`, the forth and the fifth one are the empirical parameters.
shc_new = stable_window_dec.ApplyTo(shc_initial)
~~~

## GaussianFilter.py ##

This module is to weaken the high-degree(/order) error of the GRACE GSM models using a Gaussian kernel function. 
Currently, isotropy Gsuaaian filter, anisotropy Gsuaaian filter and Fan filter are available.

This module currently includes following classes:

1. `IsoGaussian`

As the error seems to be more serious especially in high degree and order SHCs, several kinds of filters based on 
Gaussian kernel were carried out, and one of them is isotropy Gaussian filter. Isotropy Gaussian filter could give the 
SHC a smaller weight by its filter radius while its degree becomes higher.

~~~python
isoGs = IsoGaussian(r=300, earth_model: EarthModel = EarthModel.general)    # parameter r is the filtering radius, earth_model decides the related geophysical constants.
shc_new = isoGs.ApplyTo(shc_initial)
~~~

2. `AniGaussian`

Further than isotropy Gaussian filter, anisotropy Gaussian filter not only considers the impact of degree but also the 
impact of order. Anisotropy Gaussian filter could give the SHC a smaller weight by two randius parameters and a 
order-truncation parameter while its degree/order becomes higher.

~~~python
aniGs = AniGaussian(r0=300, r1=500, trunc_m=30, earth_model: EarthModel = EarthModel.general)   # parameter r is the filtering radius, earth_model decides the related geophysical constants.
shc_new = isoGs.ApplyTo(shc_initial)
~~~

3. `Fan`

Fan filter also considers both degree and order, and it uses Gaussian kernel twice to give the SHCs weights

~~~python
aniGs = AniGaussian(r1=300, r2=500, earth_model: EarthModel = EarthModel.general)   # parameter r is the filtering radius, earth_model decides the related geophysical constants.
shc_new = isoGs.ApplyTo(shc_initial)
~~~

## Harmonic.py ##

This module is to do harmonic transform: synthesis a SHC signal into a grid one or analysis a grid signal into a 
SHC one by some necessary given parameters(i.e., latitudes, longitudes, Love number, associated Legendre function, 
et al.). 

This module currently includes following classes:

1. `Harmonic`
~~~python
lat = np.arange(-90, 90, 0.5)
lon = np.arange(-180, 180, 0.5)
Har = Harmonic(Parallel=-1).setEllipsoid(ell=EllipsoidType.gif48)

LN = LoveNumber('../data/Auxiliary/')
ln = LN.getNumber(180, LoveNumberType.Wang)
grid = Har.synthesis_for_SHC(shc, lat, lon)
~~~

~~~python
lat = np.arange(-90, 90, 0.5)
lon = np.arange(-180, 180, 0.5)
Har = Harmonic(Parallel=-1).setEllipsoid(ell=EllipsoidType.gif48)

nmax=60
Pnm = GeoMathKit.getPnm(lat, nmax, 1)

LN = LoveNumber('../data/Auxiliary/')
ln = LN.getNumber(180, LoveNumberType.Wang)
grid = analysis_for_Grid(self, Nmax=nmax, Inner=grid, Pnm=Pnm)
~~~

## Leakage.py ##

This module is to reduce the leakage and bias caused mainly by filtering. It contains two parts: To estimate and 
subtract the signal that leaked in the region of interest; To estimate a scale factor to restore the signal that leaked 
out of the region of interest(called bias).

This module currently includes following classes:
1. `Leakage`
~~~python
area = np.load('xxx.npy')
signal = [] # list of SHCs or Grids
filter_matrix = isoGs.getFilter()
leakage = Leakage().setArea(area).setSignal(signal).setGaussianFilter(filter_matrix)
years_frac, signals = leakage.getTS()
~~~

## GainFactor.py ##

This module is to reduce the leakage and bias together by estimating a scale factor via a given model(e.g., 
GLDAS NOAH hydrological model). Doing the same filtering to the models, then the proportion of the signal before and 
after filtering would be seen as the impact of filtering. Multiplying the filtered GRACE signal by this proportion can 
be seen as a way to restore the leakeage and bias of GRACE filtered signal.

This module currently includes following classes:

1. `SingleFactor`

To estimate a factor in a certain basin by given a series of models and certain filter. The factor is a fitting result 
of times series before and after filtering.

~~~python
ts = [] # list of grids
nmax = 60

sf = SingleFactor().setTimeSeries(timeseries).setNmax(nmax)
sf.setDecFilter(filtertype, parameters)
sf.setGaussianFilter(filtertype, parameters)
sf.setArea(area)
sf.run()

scale_factor = sf.k
~~~

2. `GridFactor`

Similar with `SingleFactor`, but fit signal in each grid point instead of a basin area.
~~~python
ts = [] # list of grids
nmax = 60

gf = GridFactor().setTimeSeries(timeseries).setNmax(nmax)
gf.setDecFilter(filtertype, parameters)
gf.setGaussianFilter(filtertype, parameters)

scale_factor_map = gf.getLMap
~~~

## ForwardModeling.py ##

This module is to reduce bias signal in an iteration way. It would get a result which after filtered can be similar with 
filtered GRACE signal, and that means the result can be seen as the initial GRACE signal(unfiltered but without noise).

This module currently includes following classes:

1. `ForwardModeling`
~~~python
fm = ForwardModeling().setModel(observed, initial)  # observed: filtered GRACE signal; initial: the iterative initial signal, if None, initial = observed
fm.setArea().setFilter(filter_method, *rs) # filter_method: Gaussian type 
fm.setNmax(60).setLoopTime(50)
fm.run()

signal = fm.model_tru
~~~



## TimeSeriesAnalysis.py ##

This module is to 

This module currently includes following classes:

## Plot.py ##

This module is to 

This module currently includes following classes: