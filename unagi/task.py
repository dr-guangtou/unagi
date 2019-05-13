#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from collections.abc import Iterable

import numpy as np
import astropy.units as u
from astropy import wcs
from astropy.io import fits
from astropy.visualization import make_lupton_rgb

from .hsc import Hsc
from .utils import r_phy_to_ang

__all__ = ['hsc_tricolor']

ANG_UNITS = ['arcsec', 'arcsecond', 'arcmin', 'arcminute', 'deg']
PHY_UNITS = ['pc', 'kpc', 'Mpc']

def hsc_tricolor(coord, cutout_size=10.0 * u.Unit('arcsec'), filters='gri',
                 dr='dr2', rerun='s18a_wide', redshift=None, cosmo=None,
                 prefix=None, save=False, verbose=True, rgb_order=False, hdu=1,
                 rgb_save=False, rgb_q=15, rgb_stretch=0.5, rgb_min=0):
    """
    Generate HSC 3-color picture using coadd image.
    """
    # Login to HSC archive
    archive = Hsc(dr=dr, rerun=rerun)

    # List of three filters
    filter_list = list(filters)
    if len(filter_list) is not 3:
        raise ValueError("# Need and only need three filters!")
    # Check the choices of filters
    assert np.all(
        [(f in archive.FILTER_SHORT) or (f in archive.FILTER_LIST)
         for f in filter_list]), '# Wrong filter choice!'

    # Parse the cutout image size.
    # We use central coordinate and half image size as the default format.
    if not isinstance(cutout_size, u.quantity.Quantity):
        if verbose:
            print("# Assume the cutout size is in arcsec unit.")
        ang_size = cutout_size * u.Unit('arcsec')
    else:
        cutout_unit = cutout_size.unit
        if str(cutout_unit) in ANG_UNITS:
            ang_size = cutout_size.to(u.Unit('arcsec'))
        elif str(cutout_unit) in PHY_UNITS:
            if redshift is None:
                raise ValueError("# Need to provide redshift value to use physical size!")
            elif (redshift < 0.) or (~np.isfinite(redshift)):
                raise ValueError("# Redshift value is not valid!")
            else:
                ang_size = r_phy_to_ang(cutout_size, redshift, cosmo=cosmo)
        else:
            raise ValueError("# Wrong unit for cutout size: {}".format(str(cutout_unit)))

    # Output file names
    if prefix is None:
        # TODO: More informed prefix
        preifx = 'hsc_tricolor'

    # List of fits file
    fits_list = ['_'.join([prefix, f]) + '.fits' for f in filter_list]
    # Availability of each file
    fits_save = [os.path.isfile(f) or os.path.islink(f) for f in fits_list]

    # Name of the output JPEG picture (if necessary)
    if rgb_save:
        # TODO: show filter names in the filename
        rgb_jpg = prefix + '_color.jpg'
    else:
        rgb_jpg = None

    # List of RGB data
    rgb_cube = []

    # Load the cutout images in three bands
    for ii, filt in enumerate(filter_list):
        if fits_save[ii]:
            cutout_hdu = fits.open(fits_list[ii])
        else:
            cutout_hdu = archive.get_cutout_image(
                coord, w_half=ang_size, h_half=ang_size, filt=filt)
            if save:
                _ = cutout_hdu.writeto(fits_list[ii], overwrite=True)

        if ii == 0:
            cutout_wcs = wcs.WCS(cutout_hdu[hdu].header)

        rgb_cube.append(cutout_hdu[hdu].data)

    # If the filters are not in RGB order, reverse it
    if not rgb_order:
        rgb_cube.reverse()

    cutout_rgb = make_lupton_rgb(rgb_cube[0], rgb_cube[1], rgb_cube[2],
                                 Q=rgb_q, stretch=rgb_stretch, minimum=rgb_min,
                                 filename=rgb_jpg)

    return cutout_rgb, cutout_wcs
