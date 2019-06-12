#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions about using HSC catalogs."""

import numpy as np

import astropy.units as u
from astropy.table import Column

__all__ = ['moments_to_shape', 'abmag_to_image', 'world_to_image', 'select_clean_objects',
           'objects_to_galsim', 'mag_to_flux']

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

def mag_to_flux(catalog, mag_col, zeropoint=27.0, update=True):
    """
    Convert AB magnitude into HSC image flux unit.
    """
    flux = 10.0 ** ((zeropoint - catalog[mag_col]) / 2.5)
    if update:
        flux_col = mag_col.replace('mag', 'flux')
        if flux_col in catalog.colnames:
            catalog.remove_column(flux_col)
        catalog.add_column(data=flux, name=flux_col)
        return catalog
    return flux

def world_to_image(catalog, wcs, ra='ra', dec='dec', update=True):
    """
    Get the X, Y coordinate on the image.
    """
    xy_arr = np.asarray(
        [wcs.all_world2pix(catalog[ra][ii], catalog[dec][ii], 0)
         for ii in np.arange(len(catalog))])
    x_arr, y_arr = xy_arr[:, 0], xy_arr[:, 1]

    if update:
        if 'x' in catalog.colnames:
            catalog.remove_column('x')
        if 'y' in catalog.colnames:
            catalog.remove_column('y')
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
        if rad_col in catalog.colnames:
            catalog.remove_column(rad_col)
        catalog.add_column(Column(data=rad, name=rad_col))
        if ell_col in catalog.colnames:
            catalog.remove_column(ell_col)
        catalog.add_column(Column(data=ell, name=ell_col))
        if theta_col in catalog.colnames:
            catalog.remove_column(theta_col)
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

def objects_to_galsim(img, objects, psf_model=None, extended='i_extendedness',
                      psf_mag='psf_mag', gal_mag='cmodel_mag',
                      exp_mag='cmodel_exp_mag', dev_mag='cmodel_dev_mag'):
    """
    Convert the HSC objects into GalSim model image.
    """
    try:
        import galsim
    except ImportError:
        raise Exception("# Please install GalSim first!")

    # Shape of the image
    img_shape = img.shape

    # Empty image
    img_empty = np.zeros(img_shape)

    # This is the "canvas" of the model
    img_model = galsim.ImageF(img_shape[1], img_shape[0])

    # The "true" center of the image is allowed to be halfway between two pixels, as is the
    # case for even-sized images.  full_image.center is an integer position,
    # which would be 1/2 pixel up and to the right of the true center in this case.
    img_center = img_model.true_center

    # Pass the WCS from the header of the image to GalSim
    # img_wcs = galsim.AstropyWCS(wcs=cutout_wcs)

    # Pass the PSF image to GalSim
    if psf_model is not None:
        psf_obj = galsim.InterpolatedImage(
            galsim.image.Image(psf_model), scale=1.0)
    else:
        psf_obj = None

    # Go through the list of objects, there appears to be no other way than a for loop

    for ii, obj in enumerate(objects):
        # Position of the image
        img_position = galsim.PositionD(obj['x'], obj['y'])

        # Get the offset of the object on the canvas
        obj_offset = galsim.PositionD(
            img_position.x - img_center.x + 1.,
            img_position.y - img_center.y + 1.)

        # Seperate the type of the object
        if obj[extended] < 0.5:
            if np.isfinite(obj[psf_mag]) & (obj[psf_mag] > 0):
                # Generate a star
                psf_flux = abmag_to_image(obj[psf_mag])
                psf = psf_obj.withFlux(psf_flux)
                # Add it to the empty image
                img_empty += psf.drawImage(img_model, offset=obj_offset).array
            else:
                print("# Cannot generate star {} with magnitude {}".format(
                    ii, obj[psf_mag]))
        else:
            if (np.isfinite(obj[exp_mag]) & (obj[exp_mag] > 0) &
                    np.isfinite(obj[dev_mag]) & (obj[dev_mag] > 0) &
                    np.isfinite(obj[gal_mag]) & (obj[gal_mag] > 0)):
                try:
                    # The exponential component
                    flux_exp = abmag_to_image(obj[exp_mag])
                    shape_exp = galsim.Shear(
                        q=obj['cmodel_exp_ellipse_ba'],
                        beta=obj['cmodel_exp_ellipse_theta'] * galsim.degrees)
                    comp_exp = galsim.Exponential(
                        half_light_radius=obj['cmodel_exp_ellipse_r'], flux=flux_exp)
                    comp_exp = comp_exp.shear(shape_exp)

                    # The De Vacouleurs component
                    flux_dev = abmag_to_image(obj[dev_mag])
                    shape_dev = galsim.Shear(
                        q=obj['cmodel_dev_ellipse_ba'],
                        beta=obj['cmodel_dev_ellipse_theta'] * galsim.degrees)
                    comp_dev = galsim.DeVaucouleurs(
                        half_light_radius=obj['cmodel_dev_ellipse_r'], flux=flux_dev,
                        trunc=(6.0 * obj['cmodel_dev_ellipse_r']))
                    comp_dev = comp_dev.shear(shape_dev)

                    # Combine the two component
                    cmodel = galsim.Add([comp_exp, comp_dev])
                    if psf_obj is not None:
                        # Convolution with PSF
                        cmodel = galsim.Convolve([cmodel, psf_obj])

                    # Add it to the empty image
                    img_empty += cmodel.drawImage(
                        img_model, method='no_pixel', offset=obj_offset).array
                except Exception:
                    print("# Cannot generate galaxy {} with magnitude {}".format(
                        ii, obj[gal_mag]))
            else:
                print("# Problematic galaxy {} with magnitude {}".format(
                    ii, obj[gal_mag]))

    return img_empty
