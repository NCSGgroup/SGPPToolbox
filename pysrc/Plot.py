import os

import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import AxesGrid

from pysrc.Grid import Grid
from pysrc.Setting import ProjectionType

mpl.rcParams['font.family'] = 'STSong'
mpl.rcParams['font.size'] = 20
plt.rcParams['axes.unicode_minus'] = False


def _save_or_show(save: str = None):
    if save is not None:

        pic_format = save.split('.')[-1]
        if pic_format == 'eps':
            save_dpi = None
        else:
            save_dpi = 450

        path_list = save.split('/')
        path = ''
        for i in range(len(path_list) - 1):
            path += path_list[i] + '/'
        if not os.path.exists(path):
            os.makedirs(path)
        plt.savefig(save, dpi=save_dpi, bbox_inches='tight', format=pic_format)
    else:
        plt.show()

    pass


def plot_grid(*grid, area=None, save=None, projection=ProjectionType.PlateCarree):
    """

    :param projection:
    :param area: list, area to plot, [lon_start, lon_end, lat_start, lat_end]
    :param grid: list, [grid, title, vmin, vmax]
    :param save: save path and file name
    :return:
    """
    n = len(grid)
    if n == 1:
        rc = (1, 1)
    elif n == 2:
        rc = (1, 2)
    elif n == 3:
        rc = (1, 3)
    elif n == 4:
        rc = (2, 2)
    else:
        raise Exception('too many to plot!')

    lat, lon = grid[0][0].lat, grid[0][0].lon

    grids = []
    for i in range(len(grid)):
        if type(grid[i][0]) == Grid:
            grids.append(grid[i][0].map)
        else:
            grids.append(grid[i][0])

    if projection == ProjectionType.Mollweide:
        projection = ccrs.Mollweide(central_longitude=180)

    elif projection == ProjectionType.Orthographic:
        projection = ccrs.Orthographic(central_latitude=-90)

    elif projection == ProjectionType.SouthPolarStereo:
        projection = ccrs.SouthPolarStereo()

    else:
        projection = ccrs.PlateCarree(central_longitude=180)

    transform = ccrs.PlateCarree()
    axes_class = (GeoAxes, dict(map_projection=projection))
    fig1 = plt.figure(figsize=(12, 8))
    axgr = AxesGrid(fig1, 111, axes_class=axes_class,
                    nrows_ncols=rc,
                    axes_pad=1,
                    cbar_location='bottom',
                    cbar_mode='each',
                    cbar_pad=0.15,
                    cbar_size='3%',
                    label_mode='')
    for i in range(n):
        ax = axgr[i]
        lon2d, lat2d = np.meshgrid(lon, lat)
        p = ax.pcolormesh(lon2d, lat2d, grids[i], cmap="RdBu", transform=transform, vmin=grid[i][2],
                          vmax=grid[i][3])
        # p = ax.contour(lon2d, lat2d, grids[i], cmap="RdBu", transform=transform, vmin=grid[i][2],
        #                   vmax=grid[i][3])
        # ax.gridlines(draw_labels=True, dms=False, x_inline=None, y_inline=None, auto_inline=True)
        # ax.gridlines()
        ax.coastlines(resolution='110m')

        if area is not None:
            ax.set_extent(area)  # show area result

        cb = axgr.cbar_axes[i].colorbar(p, extend='both')
        cb.set_label_text(grid[i][1])

    _save_or_show(save)

    return None


def plot_time_series(year_fracs, signals: list, xticks=None, yticks=None, gridline=True, save=None,
                     figure_size=(14, 6), gap=True):
    """

    :param figure_size:
    :param gap:
    :param year_fracs:
    :param signals:  list of dict: {'signal': iterable, 'label': str, 'color': str, 'linestyle': str, 'marker': str}
    :param xticks:
    :param yticks:
    :param gridline:
    :param save:
    :return:
    """

    plt.figure(figsize=figure_size)
    plt.style.use('science')
    plt.rcParams.update({'font.size': 20})
    mpl.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    if xticks is None:
        begin_year, end_year = int(year_fracs[0]), int(year_fracs[-1]) + 1
        xticks = np.arange(begin_year, end_year + 1)
    assert year_fracs[0] >= xticks[0] and year_fracs[-1] <= xticks[-1]

    legend_flag = False
    for i in range(len(signals)):
        if not legend_flag and 'label' in signals[i].keys():
            legend_flag = True

        this_year_fracs = signals[i]['year_fracs']
        this_signal = signals[i]['signal']
        label = signals[i]['label'] if 'label' in signals[i].keys() else None
        color = signals[i]['color'] if 'color' in signals[i].keys() else None
        linestyle = signals[i]['linestyle'] if 'linestyle' in signals[i].keys() else None
        marker = signals[i]['marker'] if 'marker' in signals[i].keys() else None

        if gap:
            index1 = np.where(this_year_fracs < 2018)
            index2 = np.where(this_year_fracs > 2018)
            plt.plot(year_fracs[index1], this_signal[index1], label=label, marker=marker, color=color,
                     linestyle=linestyle)
            plt.plot(year_fracs[index2], this_signal[index2], marker=marker, color=color, linestyle=linestyle)

            plt.bar([2017.9479452054793], [9999], bottom=[-999], width=1.0136986301372417, color='grey')
        else:
            plt.plot(this_year_fracs, this_signal, label=label, marker=marker, color=color, linestyle=linestyle)

    if legend_flag:
        plt.legend()

    begin_x, end_x = xticks[0], xticks[-1]
    plt.xlim(begin_x, end_x)
    plt.xticks(xticks)

    if yticks is not None:
        begin_y, end_y = yticks[0], yticks[-1]
        plt.ylim(begin_y, end_y)
        plt.yticks(yticks)

    plt.xlabel('Time(year)')
    plt.ylabel('EWH(mm)')

    if gridline:
        plt.grid(ls=':')

    _save_or_show(save)
