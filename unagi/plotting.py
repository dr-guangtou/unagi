#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Making pretty plots."""

import numpy as np

from astropy.visualization import ZScaleInterval, \
    AsymmetricPercentileInterval

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rcParams
from matplotlib import gridspec
from matplotlib.colorbar import Colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rc('text', usetex=True)
rcParams.update({'axes.linewidth': 1.5})
rcParams.update({'xtick.direction': 'in'})
rcParams.update({'ytick.direction': 'in'})
rcParams.update({'xtick.minor.visible': 'True'})
rcParams.update({'ytick.minor.visible': 'True'})
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '8.0'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '4.0'})
rcParams.update({'xtick.minor.width': '1.5'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '8.0'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '4.0'})
rcParams.update({'ytick.minor.width': '1.5'})
rcParams.update({'axes.titlepad': '10.0'})
rcParams.update({'font.size': 25})

__all__ = ['FILTERS_COLOR', 'plot_skyobj_hist', 'map_skyobjs', 'random_cmap',
           'display_single', 'display_all']

# Default colormaps
IMG_CMAP = plt.get_cmap('viridis')
IMG_CMAP.set_bad(color='black')

FILTERS_COLOR = ['#2ca02c', '#ff7f0e', '#d62728', '#8c564b', '#7f7f7f']
FILTERS_SHORT = ['g', 'r', 'i', 'z', 'y']

def plot_skyobj_hist(X, summary, filt, prop, region=None, aper=None, fontsize=20):
    """Making 1-D summary plot of the sky objects."""
    # Range of the X-axis
    x_range = [summary['low'], summary['upp']]

    # Color for the filter
    color_use = FILTERS_COLOR[FILTERS_SHORT.index(filt)]

    # Start the figure
    fig = plt.figure(figsize=(6, 5))
    fig.subplots_adjust(left=0.01, right=0.995, bottom=0.186, top=0.996)
    ax1 = fig.add_subplot(111)

    # Grid
    ax1.grid(linestyle='--', alpha=0.3, linewidth=1.5, color='gray')

    # Vertical line to highlight 0
    ax1.axvline(
        0.0, linewidth=2.0, linestyle='-', c='gray', alpha=1.0, zorder=0)

    # Histogram
    _ = ax1.hist(X, bins='auto', density=True, histtype='stepfilled', alpha=0.4,
                 label=r'$\mathrm{Hist}$', range=x_range, color=color_use, zorder=1)

    # KDE curve
    if summary['kde'] is not None:
        x_grid = np.linspace(summary['low'], summary['upp'], 500)
        ax1.plot(x_grid, summary['kde'].evaluate(x_grid), linewidth=2.5, alpha=0.9,
                 label=r'$\mathrm{KDE}$', color=color_use, zorder=2)

    # highlight the mean and median
    ax1.axvline(summary['mean'], linewidth=2.5, linestyle='--', c=color_use,
                label=r'$\rm Mean$', zorder=3, alpha=0.9)
    ax1.axvline(summary['median'], linewidth=2.0, linestyle=':', c=color_use,
                label=r'$\rm Median$', zorder=4, alpha=0.9)

    ax1.legend(fontsize=18, loc='best')

    if prop == 'flux':
        ax1.set_xlabel(r'$\mathrm{Flux\ }[\mu\mathrm{Jy}]$', fontsize=30)
    elif prop == 'snr':
        ax1.set_xlabel(r'$\mathrm{S/N}$', fontsize=30)
    elif prop == 'mu':
        ax1.set_xlabel(r'$\mu\ [\mu\mathrm{Jy}/\mathrm{arcsec}^2]$', fontsize=30)
    else:
        raise Exception("# Wrong type of properties: flux/snr/mu")

    # Remove the Y-axis tick labels
    ax1.yaxis.set_ticklabels([])

    # Show some basic statistics
    if region is not None:
        _ = ax1.text(0.04, 0.92, region, fontsize=fontsize,
                     horizontalalignment='left', verticalalignment='center',
                     transform=ax1.transAxes)

    if aper is not None:
        _ = ax1.text(0.04, 0.83, aper, fontsize=fontsize,
                     horizontalalignment='left', verticalalignment='center',
                     transform=ax1.transAxes)

    _ = ax1.text(0.04, 0.74, r'$N:{0}$'.format(len(X)), fontsize=fontsize,
                 horizontalalignment='left', verticalalignment='center',
                 transform=ax1.transAxes)

    _ = ax1.text(0.04, 0.65, r'$\mu:{0:8.5f}$'.format(summary['mean']), fontsize=fontsize,
                 horizontalalignment='left', verticalalignment='center',
                 transform=ax1.transAxes, color='k')

    _ = ax1.text(0.04, 0.56, r'$\sigma:{0:8.5f}$'.format(summary['std']), fontsize=fontsize,
                 horizontalalignment='left', verticalalignment='center',
                 transform=ax1.transAxes, color='k')

    _ = ax1.text(0.04, 0.47, r'$\rm m:{0:8.5f}$'.format(summary['median']),
                 fontsize=fontsize, horizontalalignment='left', verticalalignment='center',
                 transform=ax1.transAxes, color='k')

    return fig

def map_skyobjs(x, y, n, mu, label=None, n_min=10, vmin=None, vmax=None,
                y_size=4, margin=0.2, fontsize=30, cbar_label=False):
    """Map the RA, Dec distributions of sky objects."""
    # Only keey the bins with enough sky objects in them
    mu[n <= n_min] = np.nan

    xy_ratio = (x.max() - x.min()) / (y.max() - y.min())

    fig = plt.figure(figsize=(xy_ratio * y_size, y_size))
    ax1 = fig.add_subplot(111)

    ax1.grid(linestyle='--', alpha=0.6)
    im = ax1.imshow(mu.T, origin='lower', extent=[x[0], x[-1], y[0], y[-1]],
                    aspect='equal', interpolation='nearest',
                    cmap=plt.get_cmap('coolwarm'), vmin=vmin, vmax=vmax)

    ax1.set_xlim(x.min() - margin, x.max() + margin)
    ax1.set_ylim(y.min() - margin, y.max() + margin)

    if label is not None:
        plt.text(0.03, 1.05, label, transform=ax1.transAxes, fontsize=38)

    # Color bar
    cb_axes = fig.add_axes([0.48, 0.90, 0.37, 0.06])
    cb = Colorbar(ax=cb_axes, mappable=im, orientation='horizontal', ticklocation='top')
    if cbar_label:
        cb.set_label(r'$\mu{\rm Jy}/\mathrm{arcsec}^2$', fontsize=25)

    _ = ax1.set_xlabel(r'$\mathrm{R.A.\ [deg]}$', fontsize=fontsize)
    _ = ax1.set_ylabel(r'$\mathrm{Dec\ [deg]}$', fontsize=fontsize)

    return fig

def random_cmap(ncolors=256, background_color='white'):
    """Random color maps.

    Generate a matplotlib colormap consisting of random (muted) colors.
    A random colormap is very useful for plotting segmentation images.

    Parameters
    ----------
    ncolors : int, optional
        The number of colors in the colormap.  The default is 256.
    random_state : int or `~numpy.random.RandomState`, optional
        The pseudo-random number generator state used for random
        sampling.  Separate function calls with the same
        ``random_state`` will generate the same colormap.

    Returns
    -------
    cmap : `matplotlib.colors.Colormap`
        The matplotlib colormap with random colors.

    Notes
    -----
    Based on: colormaps.py in photutils

    """
    prng = np.random.mtrand._rand

    h = prng.uniform(low=0.0, high=1.0, size=ncolors)
    s = prng.uniform(low=0.2, high=0.7, size=ncolors)
    v = prng.uniform(low=0.5, high=1.0, size=ncolors)

    hsv = np.dstack((h, s, v))
    rgb = np.squeeze(colors.hsv_to_rgb(hsv))

    if background_color is not None:
        if background_color not in colors.cnames:
            raise ValueError('"{0}" is not a valid background color '
                             'name'.format(background_color))
        rgb[0] = colors.hex2color(colors.cnames[background_color])

    return colors.ListedColormap(rgb)

def display_single(img,
                   pixel_scale=0.168,
                   physical_scale=None,
                   xsize=8,
                   ysize=8,
                   ax=None,
                   alpha=1.0,
                   stretch='arcsinh',
                   scale='zscale',
                   contrast=0.25,
                   no_negative=False,
                   lower_percentile=1.0,
                   upper_percentile=99.0,
                   cmap=IMG_CMAP,
                   scale_bar=True,
                   scale_bar_length=5.0,
                   scale_bar_fontsize=20,
                   scale_bar_y_offset=0.5,
                   scale_bar_color='w',
                   scale_bar_loc='left',
                   color_bar=False,
                   color_bar_loc=1,
                   color_bar_width='75%',
                   color_bar_height='5%',
                   color_bar_fontsize=18,
                   color_bar_color='w',
                   add_text=None,
                   text_fontsize=30,
                   text_color='w'):
    """Display single image.

    Parameters
    ----------
        img: np 2-D array for image

        xsize: int, default = 8
            Width of the image.

        ysize: int, default = 8
            Height of the image.

    """
    if ax is None:
        fig = plt.figure(figsize=(xsize, ysize))
        ax1 = fig.add_subplot(111)
    else:
        ax1 = ax

    # Stretch option
    if stretch.strip() == 'arcsinh':
        img_scale = np.arcsinh(img)
    elif stretch.strip() == 'log':
        if no_negative:
            img[img <= 0.0] = 1.0E-10
        img_scale = np.log(img)
    elif stretch.strip() == 'log10':
        if no_negative:
            img[img <= 0.0] = 1.0E-10
        img_scale = np.log10(img)
    elif stretch.strip() == 'linear':
        img_scale = img
    else:
        raise Exception("# Wrong stretch option.")

    # Scale option
    if scale.strip() == 'zscale':
        try:
            zmin, zmax = ZScaleInterval(contrast=contrast).get_limits(img_scale)
        except IndexError:
            # TODO: Deal with problematic image
            zmin, zmax = -1.0, 1.0
    elif scale.strip() == 'percentile':
        try:
            zmin, zmax = AsymmetricPercentileInterval(
                lower_percentile=lower_percentile,
                upper_percentile=upper_percentile).get_limits(img_scale)
        except IndexError:
            # TODO: Deal with problematic image
            zmin, zmax = -1.0, 1.0
    else:
        zmin, zmax = np.nanmin(img_scale), np.nanmax(img_scale)

    show = ax1.imshow(img_scale, origin='lower', cmap=cmap,
                      vmin=zmin, vmax=zmax, alpha=alpha)

    # Hide ticks and tick labels
    ax1.tick_params(
        labelbottom=False,
        labelleft=False,
        axis=u'both',
        which=u'both',
        length=0)

    # Put scale bar on the image
    (img_size_x, img_size_y) = img.shape
    if physical_scale is not None:
        pixel_scale *= physical_scale
    if scale_bar:
        if scale_bar_loc == 'left':
            scale_bar_x_0 = int(img_size_x * 0.04)
            scale_bar_x_1 = int(img_size_x * 0.04 +
                                (scale_bar_length / pixel_scale))
        else:
            scale_bar_x_0 = int(img_size_x * 0.95 -
                                (scale_bar_length / pixel_scale))
            scale_bar_x_1 = int(img_size_x * 0.95)

        scale_bar_y = int(img_size_y * 0.10)
        scale_bar_text_x = (scale_bar_x_0 + scale_bar_x_1) / 2
        scale_bar_text_y = (scale_bar_y * scale_bar_y_offset)
        if physical_scale is not None:
            scale_bar_text = r'$%d\ \mathrm{kpc}$' % int(scale_bar_length)
        else:
            scale_bar_text = r'$%d^{\prime\prime}$' % int(scale_bar_length)
        scale_bar_text_size = scale_bar_fontsize

        ax1.plot(
            [scale_bar_x_0, scale_bar_x_1], [scale_bar_y, scale_bar_y],
            linewidth=3,
            c=scale_bar_color,
            alpha=1.0)
        ax1.text(
            scale_bar_text_x,
            scale_bar_text_y,
            scale_bar_text,
            fontsize=scale_bar_text_size,
            horizontalalignment='center',
            color=scale_bar_color)
    if add_text is not None:
        text_x_0 = int(img_size_x*0.08)
        text_y_0 = int(img_size_y*0.80)
        ax1.text(
            text_x_0, text_y_0, r'$\mathrm{'+add_text+'}$',
            fontsize=text_fontsize, color=text_color)

    # Put a color bar on the image
    if color_bar:
        ax_cbar = inset_axes(ax1,
                             width=color_bar_width,
                             height=color_bar_height,
                             loc=color_bar_loc)
        if ax is None:
            cbar = plt.colorbar(show, ax=ax1, cax=ax_cbar,
                                orientation='horizontal')
        else:
            cbar = plt.colorbar(show, ax=ax, cax=ax_cbar,
                                orientation='horizontal')

        cbar.ax.xaxis.set_tick_params(color=color_bar_color)
        cbar.ax.yaxis.set_tick_params(color=color_bar_color)
        cbar.outline.set_edgecolor(color_bar_color)
        plt.setp(plt.getp(cbar.ax.axes, 'xticklabels'),
                 color=color_bar_color, fontsize=color_bar_fontsize)
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'),
                 color=color_bar_color, fontsize=color_bar_fontsize)

    if ax is None:
        return fig
    return ax1


def display_all(img_list, n_column=3, img_size=3., hdu_index=None, label_list=None,
                label_x=0.1, label_y=0.9, fontsize=20, **kwargs):
    """Display a list of images."""
    # Number of image to show
    n_img = len(img_list)

    if n_img <= n_column:
        n_col = n_img
        n_row = 1
    else:
        n_col = n_column
        n_row = int(np.ceil(n_img / n_column))

    fig = plt.figure(figsize=(img_size * n_col, img_size * n_row))
    fig.subplots_adjust(left=0., right=1., bottom=0., top=1., wspace=0., hspace=0.)

    gs = gridspec.GridSpec(n_row, n_col)
    gs.update(wspace=0.0, hspace=0.00)

    for ii in range(n_img):
        if hdu_index is None:
            img_show = img_list[ii]
        else:
            img_show = img_list[ii][hdu_index].data

        ax = plt.subplot(gs[ii])
        ax = display_single(img_show, ax=ax, **kwargs)

        if label_list is not None:
            if len(label_list) != n_img:
                print("# Wrong number for labels!")
            else:
                ax.text(label_x, label_y, label_list[ii], fontsize=fontsize,
                        transform=ax.transAxes, color='w')

    return fig
