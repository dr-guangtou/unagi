#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions about using HSC catalogs."""

import numpy as np

import astropy.units as u
from astropy.table import Column

__all__ = ['moments_to_shape', 'abmag_to_image', 'world_to_image', 'select_clean_objects',
           'objs_to_galsim']

# Flux unit in HSC catalog
FLUX_UNIT_S16A = (u.erg / u.s / u.Hz / u.cm ** 2)
FLUX_UNIT_S17A = (u.erg / u.s / u.Hz / u.cm ** 2)
FLUX_UNIT_S18A = (u.erg / u.s / u.Hz / u.cm ** 2)
FLUX_UNIT_PDR1 = (u.Jansky * 1.E-9)
FLUX_UNIT_PDR2 = (u.Jansky * 1.E-9)

def abmag_to_image(abmag):
    """
    Convert AB magnitude into HSC image flux unit.
    """
    return 10.0 ** ((27.0 - abmag) / 2.5)

def world_to_image(catalog, wcs, ra='ra', dec='dec', update=True):
    """
    Get the X, Y coordinate on the image.
    """
    xy_arr = np.asarray(
        [wcs.all_world2pix(catalog[ra][ii], catalog[dec][ii], 0)
         for ii in np.arange(len(catalog))])
    x_arr, y_arr = xy_arr[:, 0], xy_arr[:, 1]

    if update:
        catalog.add_column(Column(data=x_arr, name='x'))
        catalog.add_column(Column(data=y_arr, name='y'))
        return catalog
    return x_arr, y_arr

def moments_to_shape(catalog, shape_type='i_sdssshape', axis_ratio=False,
                     radian=False, update=True, to_pixel=False):
    """
    Convert the 2nd moments into elliptical shape: radius, ellipticity, position angle.
    """
    try:
        xx = catalog["{}_11".format(shape_type)]
        yy = catalog["{}_22".format(shape_type)]
        xy = catalog["{}_12".format(shape_type)]
    except KeyError:
        print("Wrong column name!")
        raise

    e1 = (xx - yy) / (xx + yy)
    e2 = (2.0 * xy / (xx + yy))
    # Get the r50 or determinant radius
    rad = np.sqrt(xx + yy)
    rad = rad / 0.168 if to_pixel else rad
    # Ellipticity or axis ratio
    ell = np.sqrt(e1 ** 2.0 + e2 ** 2.0)
    ell = 1.0 - ell if axis_ratio else ell
    # Position angle in degree or radian
    theta = (-0.5 * np.arctan2(e2, e1))
    theta = (theta * 180. / np.pi) if not radian else theta

    if update:
        rad_col = "{}_r".format(shape_type)
        theta_col = "{}_theta".format(shape_type)
        if axis_ratio:
            ell_col = "{}_ba".format(shape_type)
        else:
            ell_col = "{}_e".format(shape_type)
        catalog.add_column(Column(data=rad, name=rad_col))
        catalog.add_column(Column(data=ell, name=ell_col))
        catalog.add_column(Column(data=theta, name=theta_col))
        return catalog
    return rad, ell, theta

def select_clean_objects(catalog, check_flag='gri', check_psf='i', check_cmodel='i',
                         return_catalog=False, verbose=False):
    """
    Select the "clean" objects.
    """
    clean_mask = np.ones(len(catalog)).astype(np.bool)

    # Check data quality
    if check_flag is not None:
        for f in check_flag:
            mask_of_this_band = (
                np.isfinite(catalog['{}_extendedness'.format(f)]) &
                ~catalog['{}_flag_edge'.format(f)] &
                ~catalog['{}_flag_saturated_cen'.format(f)] &
                ~catalog['{}_flag_interpolated_cen'.format(f)]
                )
            clean_mask = clean_mask & mask_of_this_band

    # Check PSF flux/magnitude
    if check_psf is not None and check_psf in 'grizy':
        psf_mag = '{}_psf_mag'.format(check_psf)
        psf_flux = '{}_psf_flux'.format(check_psf)
        if psf_mag in catalog.colnames:
            psf_mask = (np.isfinite(catalog['{}_psf_mag'.format(check_psf)]) &
                        (catalog['{}_psf_mag'.format(check_psf)] > 0))
        elif psf_flux in catalog.colnames:
            psf_mask = (np.isfinite(catalog['{}_psf_flux'.format(check_psf)]) &
                        (catalog['{}_psf_flux'.format(check_psf)] > 0))
        else:
            raise KeyError("# PSF flux/mag not available!")
        clean_mask = clean_mask & psf_mask

    # Check PSF flux/magnitude
    if check_cmodel is not None and check_cmodel in 'grizy':
        cmodel_mag = '{}_cmodel_mag'.format(check_cmodel)
        cmodel_flux = '{}_cmodel_flux'.format(check_cmodel)
        if cmodel_mag in catalog.colnames:
            cmodel_mask = (
                np.isfinite(catalog['{}_cmodel_mag'.format(check_cmodel)]) &
                (catalog['{}_cmodel_mag'.format(check_cmodel)] > 0))
        elif cmodel_flux in catalog.colnames:
            cmodel_mask = (
                np.isfinite(catalog['{}_cmodel_flux'.format(check_cmodel)]) &
                (catalog['{}_cmodel_flux'.format(check_cmodel)] > 0))
        else:
            raise KeyError("# CModel flux/mag not available!")
        clean_mask = clean_mask & cmodel_mask

    if verbose:
        print("# {}/{} objects are clean.".format(clean_mask.sum(), len(catalog)))

    if return_catalog:
        return catalog[clean_mask], clean_mask
    return clean_mask

def objs_to_galsim(catalog):
    """
    Evaluate HSC PSF or CModel objects as a GalSim 2-D model.
    """
    try:
        import galsim
    except ImportError:
        raise Exception("# Please install GalSim first!")
    pass
