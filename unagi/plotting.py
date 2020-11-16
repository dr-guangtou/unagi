#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Making pretty plots."""

import copy

import numpy as np

from astropy import wcs
from astropy.visualization import ZScaleInterval, \
    AsymmetricPercentileInterval

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
from matplotlib.patches import Ellipse
from matplotlib.colorbar import Colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from . import catalog


__all__ = ['FILTERS_COLOR', 'plot_skyobj_hist', 'map_skyobjs', 'random_cmap',
           'display_single', 'display_all', 'overplot_all', 'shape_to_ellipse',
           'cutout_show_objects', 'setup']

# Default colormaps
IMG_CMAP = copy.copy(mpl.cm.get_cmap("viridis"))
IMG_CMAP.set_bad(color='black')

FILTERS_COLOR = ['#2ca02c', '#ff7f0e', '#d62728', '#8c564b', '#7f7f7f']
FILTERS_SHORT = ['g', 'r', 'i', 'z', 'y']

def setup(style='default', fontsize=25, linewidth=1.5, latex=False):
    """Set aesthetic parameters in one step.
    """
    from matplotlib import rcParams
    # Use LaTeX font
    if latex:
        plt.rc('text', usetex=True)

    # General look of the plot
    if style == 'default':
        # The default style is just the one used by Song Huang
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
    else:
        raise KeyError("Available style: [default]")

    # Other individual parameters
    rcParams.update({'axes.linewidth': linewidth})
    rcParams.update({'font.size': fontsize})


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
                   zmin=None,
                   zmax=None,
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
        if zmin is not None:
            zmin = np.arcsinh(zmin)
        if zmax is not None:
            zmax = np.arcsinh(zmax)
    elif stretch.strip() == 'log':
        if no_negative:
            img[img <= 0.0] = 1.0E-10
        img_scale = np.log(img)
        if zmin is not None:
            zmin = np.log(zmin)
        if zmax is not None:
            zmax = np.log(zmax)
    elif stretch.strip() == 'log10':
        if no_negative:
            img[img <= 0.0] = 1.0E-10
        img_scale = np.log10(img)
        if zmin is not None:
            zmin = np.log10(zmin)
        if zmax is not None:
            zmax = np.log10(zmax)
    elif stretch.strip() == 'linear':
        img_scale = img
    else:
        raise Exception("# Wrong stretch option.")

    # Scale option
    if scale.strip() == 'zscale':
        try:
            vmin, vmax = ZScaleInterval(contrast=contrast).get_limits(img_scale)
        except IndexError:
            # TODO: Deal with problematic image
            vmin, vmax = -1.0, 1.0
    elif scale.strip() == 'percentile':
        try:
            vmin, vmax = AsymmetricPercentileInterval(
                lower_percentile=lower_percentile,
                upper_percentile=upper_percentile).get_limits(img_scale)
        except IndexError:
            # TODO: Deal with problematic image
            vmin, vmax = -1.0, 1.0
    elif scale.strip() == 'minmax':
        vmin, vmax = np.nanmin(img_scale), np.nanmax(img_scale)
    else:
        vmin, vmax = np.nanmin(img_scale), np.nanmax(img_scale)

    if zmin is not None:
        vmin = zmin
    if zmax is not None:
        vmax = zmax

    show = ax1.imshow(img_scale, origin='lower', cmap=cmap, interpolation='none',
                      vmin=vmin, vmax=vmax, alpha=alpha)

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
                cmap_list=None, label_x=0.1, label_y=0.9, fontsize=20, fontcolor='k',
                hdu_list=False, hdu_start=1, **kwargs):
    """Display a list of images."""
    if not isinstance(img_list, list):
        raise TypeError("Provide a list of image to show or use display_single()")

    # Make a numpy array list if the input is HDUList
    if hdu_list:
        img_list = [img_list[ii].data for ii in np.arange(len(img_list))[hdu_start:]]

    if cmap_list is not None:
        assert len(cmap_list) == len(img_list), "Wrong number of color maps!"

    if label_list is not None:
        assert len(label_list) == len(img_list), "Wrong number of labels!"

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
        if cmap_list is not None:
            ax = display_single(img_show, cmap=cmap_list[ii], ax=ax, **kwargs)
        else:
            ax = display_single(img_show, ax=ax, **kwargs)

        if label_list is not None:
            if len(label_list) != n_img:
                print("# Wrong number for labels!")
            else:
                ax.text(label_x, label_y, label_list[ii], fontsize=fontsize,
                        transform=ax.transAxes, color=fontcolor)

    return fig

def overplot_all(img_list, xsize=6, ysize=6, stretch='arcsinh', scale='minmax',
                 alpha=0.7, vmin=None, cmap_list=None, hdu_index=None, alpha_list=None,
                 **kwargs):
    """
    Display a list of images.
    """
    if not isinstance(img_list, list):
        raise TypeError("Provide a list of image to show or use display_single()")
    # Number of image to show
    n_img = len(img_list)

    if cmap_list is not None:
        assert len(cmap_list) == len(img_list), "Wrong number of color maps!"

    if alpha_list is not None:
        assert len(alpha_list) == len(img_list), "Wrong number of alpha list!"

    fig = plt.figure(figsize=(xsize, ysize))
    fig.subplots_adjust(left=0., right=1., bottom=0., top=1., wspace=0., hspace=0.)
    ax = plt.subplot(111)

    for ii in range(n_img):
        if hdu_index is None:
            img = img_list[ii]
        else:
            img = img_list[ii][hdu_index].data

        if cmap_list is not None:
            cmap = cmap_list[ii]
        else:
            cmap = None

        if alpha_list is not None:
            alpha_use = alpha_list[ii]
        else:
            alpha_use = alpha

        ax = display_single(
            img, xsize=xsize, ysize=ysize, stretch=stretch, alpha=alpha_use,
            ax=ax, cmap=cmap, zmin=vmin, scale=scale, **kwargs)

    # TODO: Add legend

    return fig

def shape_to_ellipse(x, y, re, ba, theta):
    """Convert parameters to Matplotlib Ellipse patch."""
    a, b = (re * 2.0), (re * ba * 2.0)

    ells = [Ellipse(xy=np.array([x[i], y[i]]),
                    width=np.array(a[i]),
                    height=np.array(b[i]),
                    angle=np.array(theta[i]))
            for i in range(len(x))]

    return ells

def cutout_show_objects(cutout, objs, show_weighted=True, show_bad=True, show_clean=False,
                        verbose=True, xsize=8, cmap='viridis', band='i',
                        show_sdssshape=False, show_mag=False, **kwargs):
    """
    Show the HSC photometry of objects on the cutout image.
    """
    # WCS of the image
    cutout_wcs = wcs.WCS(cutout[1].header)

    # Isolate the clean objects
    if show_clean:
        objs_use, clean_mask = catalog.select_clean_objects(
            objs, return_catalog=True, verbose=True)
    else:
        objs_use = objs
        clean_mask = catalog.select_clean_objects(
            objs, return_catalog=False, verbose=True)
    x_use, y_use = catalog.world_to_image(objs_use, cutout_wcs, update=False)

    # Get the stars and show the SDSS shape
    star_mask = objs_use['{}_extendedness'.format(band)] < 0.5
    cutout_star = objs_use[star_mask]
    x_star, y_star = x_use[star_mask], y_use[star_mask]
    r_star, ba_star, pa_star = catalog.moments_to_shape(
        cutout_star, shape_type='{}_sdssshape'.format(band), axis_ratio=True,
        to_pixel=True, update=False)
    if verbose:
        if cutout_star:
            print("# There are {} point sources on the cutout".format(len(cutout_star)))
        else:
            print("# No point source is found!")

    # Get the extended objects
    gal_mask = ~star_mask
    cutout_gal = objs_use[gal_mask]
    x_gal, y_gal = x_use[gal_mask], y_use[gal_mask]
    if verbose:
        if cutout_gal:
            print("# There are {} extended sources on the cutout".format(len(cutout_gal)))
        else:
            print("# No point source is found!")

    # Start to make figure
    img_shape = cutout[1].data.shape
    fig = plt.figure(figsize=(xsize, xsize * img_shape[1] / img_shape[0]))
    ax1 = fig.add_subplot(111)

    # Show the image
    ax1 = display_single(cutout[1].data, ax=ax1, contrast=0.1, cmap=cmap, **kwargs)

    # Show the stars
    if show_mag or not show_sdssshape:
        ax1.scatter(x_star, y_star, c='dodgerblue', s=100, marker='x',
                    linewidth=2, zorder=10)
    else:
        ellip_star = shape_to_ellipse(x_star, y_star, r_star, ba_star, pa_star)
        for e in ellip_star:
            ax1.add_artist(e)
            e.set_clip_box(ax1.bbox)
            e.set_alpha(1.0)
            e.set_edgecolor('red')
            e.set_facecolor('none')
            e.set_linewidth(2.0)

    # Use color of the ellipse to indicate the magnitude of the galaxy
    if show_mag:
        if '{}_cmodel_mag'.format(band) in objs_use.colnames:
            mag = np.asarray(cutout_gal['{}_cmodel_mag'.format(band)])
        elif '{}_cmodel_flux'.format(band) in objs_use.colnames:
            mag = np.asarray(
                -2.5 * np.log10(cutout_gal['{}_cmodel_flux'.format(band)]))
        else:
            raise KeyError("# No useful CModel mag available")
        # Get a color array
        color_arr = to_color_arr(mag, bottom=20., top=26.0)
        ell_cmap = plt.get_cmap('coolwarm_r')
    else:
        color_arr, ell_cmap = None, None

    # Show the galaxies
    if show_weighted:
        r_gal, ba_gal, pa_gal = catalog.moments_to_shape(
            cutout_gal, shape_type='cmodel_ellipse', axis_ratio=True,
            to_pixel=True, update=False)
        ellip_gal = shape_to_ellipse(x_gal, y_gal, r_gal, ba_gal, pa_gal)
        for ii, e in enumerate(ellip_gal):
            ax1.add_artist(e)
            e.set_clip_box(ax1.bbox)
            e.set_alpha(0.9)
            if color_arr is None:
                e.set_edgecolor('peru')
            else:
                e.set_edgecolor(ell_cmap(int(color_arr[ii])))
            e.set_facecolor('none')
            e.set_linewidth(2.0)
    else:
        r_exp, ba_exp, pa_exp = catalog.moments_to_shape(
            cutout_gal, shape_type='cmodel_exp_ellipse', axis_ratio=True,
            to_pixel=True, update=False)
        r_dev, ba_dev, pa_dev = catalog.moments_to_shape(
            cutout_gal, shape_type='cmodel_dev_ellipse', axis_ratio=True,
            to_pixel=True, update=False)
        ellip_exp = shape_to_ellipse(x_gal, y_gal, r_exp, ba_exp, pa_exp)
        ellip_dev = shape_to_ellipse(x_gal, y_gal, r_dev, ba_dev, pa_dev)
        for ii, e in enumerate(ellip_exp):
            ax1.add_artist(e)
            e.set_clip_box(ax1.bbox)
            e.set_alpha(1.0)
            if color_arr is None:
                e.set_edgecolor('peru')
            else:
                e.set_edgecolor(ell_cmap(int(color_arr[ii])))
            e.set_facecolor('none')
            e.set_linewidth(2.0)
        for ii, e in enumerate(ellip_dev):
            ax1.add_artist(e)
            e.set_clip_box(ax1.bbox)
            e.set_alpha(0.7)
            if color_arr is None:
                e.set_edgecolor('w')
            else:
                e.set_edgecolor(ell_cmap(int(color_arr[ii])))
            e.set_facecolor('none')
            e.set_linewidth(2.0)
            e.set_linestyle('--')
    ax1.scatter(x_gal, y_gal, c='peru', s=50, marker='+', linewidth=2.0)

    if show_mag:
        cax = fig.add_axes([0.16, 0.85, 0.24, 0.02])
        norm = mpl.colors.Normalize(vmin=20.0, vmax=26.0)
        cbar = mpl.colorbar.ColorbarBase(
            cax, cmap='coolwarm_r', norm=norm, orientation='horizontal')
        cbar.ax.tick_params(labelsize=15)

    # Show the non-clean objects
    if show_bad:
        x_dirty, y_dirty = catalog.world_to_image(
            objs[~clean_mask], cutout_wcs, update=False)
        ax1.scatter(x_dirty, y_dirty, facecolor='none', edgecolor='k',
                    s=100, marker='o', linewidth=1.5, zorder=10)

    ax1.set_xlim(0, img_shape[0])
    ax1.set_ylim(0, img_shape[1] - 1)

    return fig

def to_color_arr(arr, bottom=None, top=None):
    """
    Convert a data array to "color array" (between 0 and 1).
    """
    data = copy.deepcopy(arr)
    if top is not None:
        data[data >= top] = top
    if bottom is not None:
        data[data <= bottom] = bottom

    # Fix the NaN and Inf values
    data[~np.isfinite(data)] = bottom

    return (data - np.nanmin(data)) / (np.nanmax(data) - np.nanmin(data)) * 255
