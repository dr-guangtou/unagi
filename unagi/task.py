#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import tempfile
from collections.abc import Iterable

import tarfile
import numpy as np
import astropy.units as u
from astropy import wcs
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.visualization import make_lupton_rgb
from multiprocessing import Pool
from functools import partial
import glob

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
    """
    Generate HSC 3-color picture using coadd image.
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

def _download_cutouts(args, url=None, tmp_dir=None, output_dir=None, session=None):
    filename, ids, filter = args

    # Check if the file is archive is already download, if not, download it
    resp = session.post(url,
                        files={'list': open(filename, 'rb')},
                        stream=True)

    # Checking that access worked
    assert(resp.status_code == 200)
    output_filename = resp.headers['Content-Disposition'].split('"')[-2]

    # Proceed to download the data
    with open(os.path.join(tmp_dir, output_filename), 'wb') as f:
        for chunk in resp.iter_content(chunk_size=1024):
            f.write(chunk)

    # Untar the archive
    with tarfile.TarFile(os.path.join(tmp_dir, output_filename), "r") as tarball:
        tarball.extractall(tmp_dir)

    # Recover path to output dir
    output_path = os.path.join(tmp_dir, output_filename.split('.tar')[0])

    # Rename files based on the object ids
    fnames = glob.glob(output_path+'/*.fits')
    for fname in fnames:
        indx = int(fname.split(output_path+'/')[1].split('-')[0]) - 2
        output_filename = os.path.join(output_dir, filter, '%d.fits'%ids[indx])
        shutil.move(fname, output_filename)

    # Remove the tar file and temporary download folder
    os.remove(os.path.join(tmp_dir, output_filename))
    os.rmdir(output_path)

    # Return true that everything is ok
    return True

def hsc_cutout_bulk_download(table, cutout_size=10.0 * u.Unit('arcsec'),
                             filters='i', dr='dr2', rerun='s18a_wide', img_type='coadd',
                             verbose=True, archive=None, use_saved=False,
                             image=True, variance=False, mask=False, nproc=1,
                             tmp_dir=None, output_dir='./', **kwargs):
    """
    Generate HSC cutout images in bulk.

    table: astropy table
        Astropy table of with at least object_id, and (ra, dec) in deg.
    """
    # Login to HSC archive
    if archive is None:
        archive = Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun
        if dr[0] == 'p':
            rerun = rerun.replace(dr + '_', '')

    # Get temporary directory for dowloading and staging
    if tmp_dir is None:
        tmp_dir = tempfile.mkdtemp()

    # List of three filters
    filter_list = list(filters)

    # Check the choices of filters
    assert np.all(
        [(f in archive.FILTER_SHORT) or (f in archive.FILTER_LIST)
         for f in filter_list]), '# Wrong filter choice!'

    # Parse the cutout image size.
    # We use central coordinate and half image size as the default format.
    if isinstance(cutout_size, list):
        if len(cutout_size) != 2:
            raise Exception("# Cutout size should be like: [Width, Height]")
        ang_size_w = _get_cutout_size(
            cutout_size[0], verbose=verbose)
        ang_size_h = _get_cutout_size(
            cutout_size[1], verbose=verbose)
    else:
        ang_size_w = ang_size_h = _get_cutout_size(cutout_size, verbose=verbose)

    # Compute the number of batches to download
    # There is a hard limit of 1000 cutouts at a time
    batch_size = 1000
    n_batches = len(table) // batch_size
    if len(table) % batch_size > 0:
        n_batches = n_batches + 1

    # Step 1: Create all download files, for all filters
    batch_files = []
    for batch_index in range(n_batches):
        list_table = table[['ra', 'dec', 'object_id']][batch_index*batch_size:(batch_index+1)*batch_size]
        list_table['sw'] = str(ang_size_w.value)+'asec'
        list_table['sh'] = str(ang_size_h.value)+'asec'
        list_table['rerun'] = archive.rerun
        list_table['filter'] = archive._check_filter(filters[0])
        list_table['type'] = img_type
        list_table['image'] = 'true' if image else 'false'
        list_table['variance'] = 'true' if variance else 'false'
        list_table['mask'] = 'true' if mask else 'false'
        list_table['#?'] = ''
        # Saving object ids corresponding to the downloaded objects
        ids = list_table['object_id']
        list_table = list_table[['#?', 'ra', 'dec', 'sw', 'sh', 'filter', 'rerun', 'image', 'variance', 'mask', 'type']]

        # Generate the table for requested filter
        for f in filters:
            list_table['filter'] = archive._check_filter(f)
            filename = os.path.join(tmp_dir, ('batch_%s_%d')%(f, batch_index))
            list_table.write(filename, format='ascii.tab')
            batch_files.append((filename, ids, archive._check_filter(f)))

    # Step 2: Download fits files
    # Create output directories
    for f in filters:
        directory = os.path.join(output_dir, archive._check_filter(f))
        if not os.path.exists(directory):
            os.makedirs(directory)
    download_cutouts = partial(_download_cutouts,
                               url=archive.archive.img_url,
                               tmp_dir=tmp_dir,
                               output_dir=output_dir,
                               session=archive.session)
    with Pool(nproc) as pool:
        res = pool.map(download_cutouts, batch_files)

    return res

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
