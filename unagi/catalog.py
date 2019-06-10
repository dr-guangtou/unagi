#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions about using HSC catalogs."""

import numpy as np

import astropy.units as u
from astropy.table import Column

__all__ = ['moments_to_shape', 'abmag_to_image', 'world_to_image']

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

def moments_to_shape(xx, yy, xy, axis_ratio=False, radian=False):
    """
    Convert the 2nd moments into elliptical shape: radius, ellipticity, position angle.
    """
    e1 = (xx - yy) / (xx + yy)
    e2 = (2.0 * xy / (xx + yy))
    rad = np.sqrt(xx + yy)
    ell = np.sqrt(e1 ** 2.0 + e2 ** 2.0)
    theta = (0.5 * np.arctan2(e2, e1))
    if not radian:
        theta *= (180.0 / np.pi)

    if axis_ratio:
        return rad, 1.0 - ell, theta

    return rad, ell, theta
