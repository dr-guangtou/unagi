#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil

import numpy as np
import astropy.units as u
from astropy import wcs
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.visualization import make_lupton_rgb

from . import query
from .hsc import Hsc
from .utils import r_phy_to_ang

__all__ = ['hsc_tricolor', 'hsc_cutout', 'hsc_psf',
           'hsc_cone_search', 'hsc_box_search', 'hsc_check_coverage']

ANG_UNITS = ['arcsec', 'arcsecond', 'arcmin', 'arcminute', 'deg']
PHY_UNITS = ['pc', 'kpc', 'Mpc']

def hsc_tricolor(coord, cutout_size=10.0 * u.Unit('arcsec'), coord_2=None,
                 filters='gri', dr='dr2', rerun='s18a_wide', redshift=None,
                 cosmo=None, prefix=None, use_saved=False, save_img=False,
                 verbose=True, rgb_order=False, hdu=1, archive=None, output_dir='./',
                 save_rgb=False, rgb_q=15, rgb_stretch=0.5, rgb_min=0):
    """Generate HSC 3-color picture using coadd image.

    This task generates 3-color figure of the a HSC cutout region using the
    `make_lupton_rgb` method from `astropy.visualization`. User can pass a
    `hsc.HSC` instance directly to the task or decide which `dr` and `rerun` to use.
    User can output the figure in PNG format.

    Parameters
    ----------
    coord : astropy.skycoord.SkyCoord object
        (RA, Dec) coordinate for the center of the cutout if `coord`,
        `cutout_size` pairs are used (default). If `coord_2` is provided,
        this will be the coordinate of the lower-left corner.
    cutout_size : astropy.Quantity or a pair of Quantity objects
        Defines the size of the cutout. For a single quantity, a square
        cutout will be generated. If a list of two quantities are provided,
        they will be used as the width and height of the cutout.
    coord_2 : astropy.skycoord.SkyCoord object, optional
        (RA, Dec) Coordinate of the top-right corner of the cutout image.
    filters : str, optional
        A combination of 3 filters to make the 3-color image. Default: 'gri'.
    archive : unagi.hsc.Hsc object, optional
        HSC archive object from unagi to be used. Default: None
    dr : str, optional
        Name of the data release to be used. Default: `dr2`.
    rerun : str, optional
        Name of the rerun (data reduction) to be used. Default: `s18a_wide`.
    redshift : str, optional
        Redshift of the object. Will use only when physical `cutout_size` is
        provided. Default: None
    cosmo : astropy.cosmology object, optional
        Cosmology model to be used to convert a physical size into an angular one.
        Default: None.
    prefix : str, optional
        Prefix of the output file. Default: None
    used_saved : boolen, optional
        Whether to use the saved file if available. Default: False
    rgb_order : boolen, optional
        Order of the filters. Assume it is in 'RGB' order if True. Default: False
    hdu : int, optional
        Which HDU to be used. Default: 1 (imaeg)
    output_dir : str, optional
        Directory for output picture. Default: './'
    save_image : boolen, optional
        Whether to save the image in each band in FITS file. Default: False
    save_rgb : boolen, optional
        Output figure as a JPEG file or not. Default: False
    rgb_q : int, optional
        The asinh softening parameter used in `make_rgb_lupton`. Default: 15
    rgb_stretch : float, optional
        The linear stretch of the image used in `make_rgb_lupton`. Default: 0.5
    rgb_min: float, optional
        Intensity that should be mapped to black (a scalar or array for R, G, B). Default: 0

    Returns
    -------
    cutout_rgb : ndarray
        RGB (integer, 8-bits per channel) color image as an NxNx3 numpy array.
    cutout_wcs : astropy.wcs object
        WCS information in the cutout area.

    See Also
    --------
    hsc_cutout : task that generates cutout FITS images.

    Notes
    -----

    References
    ----------
    Document for the `make_rgb_lupton` function is here [1]_.

    .. [1] http://docs.astropy.org/en/stable/api/astropy.visualization.make_lupton_rgb.html

    Examples
    --------

    >>> from unagi import hsc
    >>> from unagi.task import hsc_tricolor
    >>> pdr2 = hsc.Hsc(dr='pdr2', rerun='pdr2_dud')
    >>> coord = SkyCoord(150.09134, 2.205916, frame='icrs', unit='deg')
    >>> s_ang = 15.0 * u.arcsec
    >>> filters = 'gri'
    >>> cutout_rgb, cutout_wcs = hsc_tricolor(coord, cutout_size=s_ang, filters=filters, archive=pdr2)
    # Retrieving cutout image in filter: g
    # Retrieving cutout image in filter: r
    # Retrieving cutout image in filter: i

    """
    # Login to HSC archive
    if archive is None:
        archive = Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun
        if dr[0] == 'p':
            rerun = rerun.replace(dr + '_', '')

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
                cutout_size, redshift=redshift, cosmo=cosmo, verbose=verbose)
    else:
        ang_size_w = ang_size_h = None

    # Output file names
    if prefix is None:
        if coord_2 is None:
            ra_str, dec_str = coord.to_string('decimal', precision=4).split(' ')
            if isinstance(cutout_size, list):
                w_half_str = "{:8.2f}{}".format(ang_size_w.value, ang_size_w.unit).strip()
                h_half_str = "{:8.2f}{}".format(ang_size_h.value, ang_size_h.unit).strip()
                size_str = '{0}_{1}'.format(w_half_str, h_half_str)
            else:
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
               dr='dr2', rerun='s18a_wide', redshift=None, cosmo=None, img_type='coadd',
               prefix=None, verbose=True, archive=None, save_output=True, use_saved=False,
               output_dir='./', **kwargs):
    """
    Generate HSC cutout images.
    """
    # Login to HSC archive
    if archive is None:
        archive = Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun
        if dr[0] == 'p':
            rerun = rerun.replace(dr + '_', '')

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
                cutout_size, redshift=redshift, cosmo=cosmo, verbose=verbose)
    else:
        ang_size_w = ang_size_h = None

    # Output file names
    if prefix is None:
        if coord_2 is None:
            ra_str, dec_str = coord.to_string('decimal', precision=4).split(' ')
            if isinstance(cutout_size, list):
                w_half_str = "{:8.2f}{}".format(ang_size_w.value, ang_size_w.unit).strip()
                h_half_str = "{:8.2f}{}".format(ang_size_h.value, ang_size_h.unit).strip()
                size_str = '{0}_{1}'.format(w_half_str, h_half_str)
            else:
                size_str = "{:8.2f}{}".format(cutout_size.value, cutout_size.unit).strip()
            prefix = '{0}_{1}_{2}_{3}_{4}'.format(dr, rerun, ra_str, dec_str, size_str)
        else:
            ra1_str, dec1_str = coord.to_string('decimal', precision=4).split(' ')
            ra2_str, dec2_str = coord_2.to_string('decimal', precision=4).split(' ')
            prefix = '{0}_{1}_ra_{2}_{3}_dec_{4}_{5}'.format(
                dr, rerun, ra1_str, dec1_str, ra2_str, dec2_str)

    # Location of the output files
    prefix = os.path.join(output_dir, prefix)

    # List of fits file
    if img_type == 'coadd':
        output_list = ['_'.join([prefix, f]) + '.fits' for f in filter_list]
    elif img_type == 'warp':
        output_list = ['_'.join([prefix, f]) + '.tar' for f in filter_list]
    else:
        raise Exception("# Wrong image type: coadd or warp !")

    # Availability of each file
    file_available = [os.path.isfile(f) or os.path.islink(f) for f in output_list]

    # A "datacube" mode for people who want to compile just multi-band images
    # together in one file. This only applies to `coadd` data type, and only works
    # when multiple filters are selected.

    # Get the cutout in each band
    cutout_list = []
    for ii, filt in enumerate(filter_list):
        if file_available[ii] and use_saved:
            if img_type == 'coadd':
                if verbose:
                    print("# Read in saved FITS file: {}".format(output_list[ii]))
                cutout_hdu = fits.open(output_list[ii])
            else:
                if verbose:
                    print("# Read in saved TAR file: {}".format(output_list[ii]))
                # TODO
                raise NotImplementedError("# Not yet...")
        else:
            if verbose:
                if img_type == 'coadd':
                    print("# Retrieving cutout image in filter: {}".format(filt))
                else:
                    print("# Retrieving warped images in filter: {}".format(filt))

            # Get the FITS data or the URL of compressed tarball
            cutout_hdu = archive.get_cutout_image(
                coord, coord_2=coord_2, w_half=ang_size_w, h_half=ang_size_h,
                filt=filt, img_type=img_type, **kwargs)

            if img_type == 'coadd' and save_output:
                # Download FITS file for coadd image.
                _ = cutout_hdu.writeto(output_list[ii], overwrite=True)

            if img_type == 'warp':
                # Download the tarball for warpped images.
                _ = shutil.move(
                    download_file(cutout_hdu, show_progress=False), output_list[ii])

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

def hsc_psf(coord, centered=True, filters='i', dr='dr2', rerun='s18a_wide',
            img_type='coadd', prefix=None, verbose=True, archive=None, save_output=True,
            use_saved=False, output_dir='./'):
    """
    Generate HSC PSF models.
    """
    # Login to HSC archive
    if archive is None:
        archive = Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun
        if dr[0] == 'p':
            rerun = rerun.replace(dr + '_', '')

    # List of three filters
    filter_list = list(filters)
    if len(filter_list) > 1 and verbose:
        print("# Will dgenerate cutouts for a list of filters:", filter_list)

    # Check the choices of filters
    assert np.all(
        [(f in archive.FILTER_SHORT) or (f in archive.FILTER_LIST)
         for f in filter_list]), '# Wrong filter choice!'

    # Output file names
    if prefix is None:
        ra_str, dec_str = coord.to_string('decimal', precision=4).split(' ')
        prefix = '{0}_{1}_{2}_{3}_psf'.format(dr, rerun, ra_str, dec_str)

    # Location of the output files
    prefix = os.path.join(output_dir, prefix)

    # List of fits file
    if img_type == 'coadd':
        output_list = ['_'.join([prefix, f]) + '.fits' for f in filter_list]
    elif img_type == 'warp':
        output_list = ['_'.join([prefix, f]) + '.tar' for f in filter_list]
    else:
        raise Exception("# Wrong image type: coadd or warp !")

    # Availability of each file
    file_available = [os.path.isfile(f) or os.path.islink(f) for f in output_list]

    # Get the cutout in each band
    psf_list = []

    for ii, filt in enumerate(filter_list):
        if file_available[ii] and use_saved:
            if img_type == 'coadd':
                if verbose:
                    print("# Read in saved FITS file: {}".format(output_list[ii]))
                psf_hdu = fits.open(output_list[ii])
            else:
                if verbose:
                    print("# Read in saved TAR file: {}".format(output_list[ii]))
                # TODO: read the files in the tarball
                raise NotImplementedError("# Not yet...")
        else:
            if verbose:
                if img_type == 'coadd':
                    print("# Retrieving coadd PSF model in filter: {}".format(filt))
                else:
                    print("# Retrieving warped PSF model in filter: {}".format(filt))

            # Get the FITS data or the URL of compressed tarball
            psf_hdu = archive.get_psf_model(
                coord, filt=filt, img_type=img_type, centered=centered)

            if img_type == 'coadd' and save_output:
                # Download FITS file for coadd image.
                _ = psf_hdu.writeto(output_list[ii], overwrite=True)

            if img_type == 'warp':
                # Download the tarball for warpped images.
                _ = shutil.move(
                    download_file(psf_hdu, show_progress=False), output_list[ii])

        # Append the HDU to the list
        psf_list.append(psf_hdu)

    if len(filter_list) == 1:
        return psf_list[0]

    return psf_list

def hsc_cone_search(coord, radius=10.0 * u.Unit('arcsec'), redshift=None,
                    archive=None, dr='pdr2', rerun='pdr2_wide', cosmo=None,
                    verbose=True, **kwargs):
    """
    Search for objects within a cone area.
    """
    if archive is None:
        archive = Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun
        if dr[0] == 'p':
            rerun = rerun.replace(dr + '_', '')

    # We use central coordinate and half image size as the default format.
    ra, dec = coord.ra.value, coord.dec.value
    rad_arcsec = _get_cutout_size(
        radius, redshift=redshift, cosmo=cosmo, verbose=verbose).to(u.Unit('arcsec'))

    objects = archive.sql_query(
        query.cone_search(ra, dec, rad_arcsec, archive=archive, **kwargs), verbose=True)

    return objects

def hsc_box_search(coord, box_size=10.0 * u.Unit('arcsec'), coord_2=None, redshift=None,
                   archive=None, dr='pdr2', rerun='pdr2_wide', cosmo=None,
                   verbose=True, **kwargs):
    """
    Search for objects within a box area.
    """
    # Login to HSC archive
    if archive is None:
        archive = Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun
        if dr[0] == 'p':
            rerun = rerun.replace(dr + '_', '')

    # We use central coordinate and half image size as the default format.
    if coord_2 is None:
        if isinstance(box_size, list):
            if len(box_size) != 2:
                raise Exception("# Cutout size should be like: [Width, Height]")
            ang_size_w = _get_cutout_size(
                box_size[0], redshift=redshift, cosmo=cosmo, verbose=verbose)
            ang_size_h = _get_cutout_size(
                box_size[1], redshift=redshift, cosmo=cosmo, verbose=verbose)
        else:
            ang_size_w = ang_size_h = _get_cutout_size(
                box_size, redshift=redshift, cosmo=cosmo, verbose=verbose)
        ra_size = ang_size_w.to(u.Unit('deg'))
        dec_size = ang_size_h.to(u.Unit('deg'))
        ra1, ra2 = coord.ra.value - ra_size.value, coord.ra.value + ra_size.value
        dec1, dec2 = coord.dec.value - dec_size.value, coord.dec.value + dec_size.value
    else:
        ra1, dec1 = coord.ra.value, coord.dec.value
        ra2, dec2 = coord_2.ra.value, coord_2.dec.value

    objects = archive.sql_query(
        query.box_search(ra1, ra2, dec1, dec2, archive=archive, **kwargs), verbose=True)

    return objects

def hsc_check_coverage(coord, dr='pdr2', rerun='pdr2_wide', archive=None, verbose=False,
                       return_filter=False):
    """
    Check if the coordinate is covered by HSC footprint.

    TODO: This is not very fast, will take at least a few seconds for one object.
          And it does not guarentee that the location has data.
    """
    # Login to HSC archive
    if archive is None:
        archive = Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun

    sql_str = query.PATCH_CONTAIN.format(rerun, coord.ra.value, coord.dec.value)

    coverage = archive.sql_query(sql_str, verbose=False)
    filter_list = list(np.unique(coverage['filter01']))

    if verbose:
        if filter_list:
            print("# Covered by {}-band".format(len(filter_list)))
            print(filter_list)
        else:
            print("# Not covered")

    if return_filter:
        return filter_list
    return coverage
