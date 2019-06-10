#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions about using HSC catalogs."""

import numpy as np

import astropy.units as u
from astropy.table import Column

__all__ = ['moments_to_shape', 'abmag_to_image', 'world_to_image', 'select_clean_objects']

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

def select_clean_objects(catalog, return_catalog=False, verbose=False):
    """
    Select the "clean" objects.
    """
    clean_mask = (
        np.isfinite(catalog['i_extendedness']) &
        ~catalog['g_flag_edge'] &
        ~catalog['r_flag_edge'] &
        ~catalog['i_flag_edge'] &
        ~catalog['z_flag_edge'] &
        ~catalog['y_flag_edge'] &
        ~catalog['g_flag_saturated_cen'] &
        ~catalog['r_flag_saturated_cen'] &
        ~catalog['i_flag_saturated_cen'] &
        ~catalog['z_flag_saturated_cen'] &
        ~catalog['y_flag_saturated_cen'] &
        ~catalog['g_flag_interpolated_cen'] &
        ~catalog['r_flag_interpolated_cen'] &
        ~catalog['i_flag_interpolated_cen'] &
        ~catalog['z_flag_interpolated_cen'] &
        ~catalog['y_flag_interpolated_cen']
    )

    if verbose:
        print("# {}/{} objects are clean.".format(clean_mask.sum(), len(catalog)))

    if return_catalog:
        return catalog[clean_mask], clean_mask
    return clean_mask
