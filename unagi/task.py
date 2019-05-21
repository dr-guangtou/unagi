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

__all__ = ['hsc_tricolor', 'hsc_cutout']

ANG_UNITS = ['arcsec', 'arcsecond', 'arcmin', 'arcminute', 'deg']
PHY_UNITS = ['pc', 'kpc', 'Mpc']

def hsc_tricolor(coord, cutout_size=10.0 * u.Unit('arcsec'), coord_2=None,
                 filters='gri', dr='dr2', rerun='s18a_wide', redshift=None,
                 cosmo=None, prefix=None, use_saved=False, save_img=False,
                 verbose=True, rgb_order=False, hdu=1, archive=None, output_dir='./',
                 save_rgb=False, rgb_q=15, rgb_stretch=0.5, rgb_min=0):
    """
    Generate HSC 3-color picture using coadd image.
    """
    # Login to HSC archive
    if archive is None:
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
    if coord_2 is None:
        if isinstance(cutout_size, list):
            if len(cutout_size) != 2:
                raise Exception("# Cutout size should be like: [Width, Height]")
            ang_size_w = _get_cutout_size(
                cutout_size[0], redshift=redshift, cosmo=cosmo, verbose=verbose)
            ang_size_h = _get_cutout_size(
                cutout_size[1], redshift=redshift, cosmo=cosmo, verbose=verbose)
        else:
            ang_size_w = ang_size_h = _get_cutout_size(
                cutout_size[0], redshift=redshift, cosmo=cosmo, verbose=verbose)
    else:
        ang_size_w = ang_size_h = None

    # Output file names
    if prefix is None:
        if coord_2 is None:
            ra_str, dec_str = coord.to_string('decimal', precision=4).split(' ')
            size_str = "{:8.2f}{}".format(cutout_size.value, cutout_size.unit).strip()
            prefix = '{0}_{1}_{2}_{3}_{4}'.format(dr, rerun, ra_str, dec_str, size_str)
        else:
            ra1_str, dec1_str = coord.to_string('decimal', precision=4).split(' ')
            ra2_str, dec2_str = coord.to_string('decimal', precision=4).split(' ')
            prefix = '{0}_{1}_ra_{2}_{3}_dec_{4}_{5}'.format(
                dr, rerun, ra1_str, dec1_str, ra2_str, dec2_str)

    # Location of the output files
    prefix = os.path.join(output_dir, prefix)

    # List of fits file
    fits_list = ['_'.join([prefix, f]) + '.fits' for f in filter_list]

    # Availability of each file
    fits_save = [os.path.isfile(f) or os.path.islink(f) for f in fits_list]

    # Name of the output JPEG picture (if necessary)
    if save_rgb:
        # Generate a string that represents the filter combination
        color_str = ''.join(
            [f if f in archive.FILTER_SHORT else archive.FILTER_SHORT[archive.FILTER_LIST.index(f)]
             for f in filter_list])
        rgb_jpg = prefix + '_{0}.jpg'.format(color_str)
        if verbose:
            print("# RGB picture will be saved as {}".format(rgb_jpg))
    else:
        rgb_jpg = None

    # List of RGB data
    rgb_cube = []

    # Load the cutout images in three bands
    for ii, filt in enumerate(filter_list):
        if fits_save[ii] and use_saved:
            if verbose:
                print("# Read in saved FITS file: {}".format(fits_save[ii]))
            cutout_hdu = fits.open(fits_list[ii])
        else:
            if verbose:
                print("# Retrieving cutout image in filter: {}".format(filt))
            cutout_hdu = archive.get_cutout_image(
                coord, coord_2=coord_2, w_half=ang_size_w, h_half=ang_size_h, filt=filt)
            if save_img:
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

def hsc_cutout(coord, coord_2=None, cutout_size=10.0 * u.Unit('arcsec'), filters='i',
               dr='dr2', rerun='s18a_wide', redshift=None, cosmo=None, use_saved=True,
               prefix=None, verbose=True, archive=None, save_fits=True, output_dir='./',
               **kwargs):
    """
    Generate HSC cutout images.
    """
    # Login to HSC archive
    if archive is None:
        archive = Hsc(dr=dr, rerun=rerun)

    # List of three filters
    filter_list = list(filters)
    if len(filter_list) > 1 and verbose:
        print("# Will dgenerate cutouts for a list of filters:", filter_list)

    # Check the choices of filters
    assert np.all(
        [(f in archive.FILTER_SHORT) or (f in archive.FILTER_LIST)
         for f in filter_list]), '# Wrong filter choice!'

    # Parse the cutout image size.
    # We use central coordinate and half image size as the default format.
    if coord_2 is None:
        if isinstance(cutout_size, list):
            if len(cutout_size) != 2:
                raise Exception("# Cutout size should be like: [Width, Height]")
            ang_size_w = _get_cutout_size(
                cutout_size[0], redshift=redshift, cosmo=cosmo, verbose=verbose)
            ang_size_h = _get_cutout_size(
                cutout_size[1], redshift=redshift, cosmo=cosmo, verbose=verbose)
        else:
            ang_size_w = ang_size_h = _get_cutout_size(
                cutout_size[0], redshift=redshift, cosmo=cosmo, verbose=verbose)
    else:
        ang_size_w = ang_size_h = None

    # Output file names
    if prefix is None:
        if coord_2 is None:
            ra_str, dec_str = coord.to_string('decimal', precision=4).split(' ')
            size_str = "{:8.2f}{}".format(cutout_size.value, cutout_size.unit).strip()
            prefix = '{0}_{1}_{2}_{3}_{4}'.format(dr, rerun, ra_str, dec_str, size_str)
        else:
            ra1_str, dec1_str = coord.to_string('decimal', precision=4).split(' ')
            ra2_str, dec2_str = coord.to_string('decimal', precision=4).split(' ')
            prefix = '{0}_{1}_ra_{2}_{3}_dec_{4}_{5}'.format(
                dr, rerun, ra1_str, dec1_str, ra2_str, dec2_str)

    # Location of the output files
    prefix = os.path.join(output_dir, prefix)

    # List of fits file
    fits_list = ['_'.join([prefix, f]) + '.fits' for f in filter_list]

    # Availability of each file
    fits_save = [os.path.isfile(f) or os.path.islink(f) for f in fits_list]

    # Get the cutout in each band
    cutout_list = []

    for ii, filt in enumerate(filter_list):
        if fits_save[ii] and use_saved:
            if verbose:
                print("# Read in saved FITS file: {}".format(fits_save[ii]))
            cutout_hdu = fits.open(fits_list[ii])
        else:
            if verbose:
                print("# Retrieving cutout image in filter: {}".format(filt))
            if save_fits:
                cutout_hdu = archive.download_cutout(
                    coord, fits_list[ii], coord_2=coord_2, w_half=ang_size_w, h_half=ang_size_h,
                    filt=filt, **kwargs)
            else:
                cutout_hdu = archive.get_cutout_image(
                    coord, coord_2=coord_2, w_half=ang_size_w, h_half=ang_size_h,
                    filt=filt, **kwargs)

        # Append the HDU to the list
        cutout_list.append(cutout_hdu)
    
    if len(filter_list) == 1:
        return cutout_list[0]

    return cutout_list


def _get_cutout_size(cutout_size, redshift=None, cosmo=None, verbose=True):
    """Parse the input for the size of the cutout."""
    if not isinstance(cutout_size, u.quantity.Quantity):
        if verbose:
            print("# Assume the cutout size is in arcsec unit.")
        cutout_size = cutout_size * u.Unit('arcsec')
        ang_size = cutout_size
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

    return ang_size
