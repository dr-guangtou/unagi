#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SQL search related functions"""

from . import hsc

__all__ = ['HELP_BASIC', 'COLUMNS_CONTAIN', 'TABLE_SCHEMA', 'PATCH_CONTAIN',
           'DR1_CLEAN', 'DR2_CLEAN', 'basic_meas_photometry',
           'basic_forced_photometry', 'column_dict_to_str', 'join_table_by_id',
           'box_search', 'cone_search']

HELP_BASIC = "SELECT * FROM help('{0}');"

TABLE_SCHEMA = "SELECT * FROM help('{0}.{1}');"

COLUMNS_CONTAIN = "SELECT * FROM help('{0}.%{1}%');"

PATCH_CONTAIN = """
    --- Find coadded patch images
    SELECT
        mosaic.tract,
        mosaic.patch,
        mosaic.filter01
    FROM
        {0}.mosaic JOIN public.skymap USING (skymap_id)
    WHERE
        patch_contains(patch_area, wcs, {1}, {2})
    ;
    """

DR2_CLEAN = [
    'g_pixelflags_edge', 'r_pixelflags_edge', 'i_pixelflags_edge',
    'z_pixelflags_edge', 'z_pixelflags_edge',
    'g_pixelflags_interpolatedcenter', 'r_pixelflags_interpolatedcenter',
    'i_pixelflags_interpolatedcenter', 'z_pixelflags_interpolatedcenter',
    'z_pixelflags_interpolatedcenter',
    'g_pixelflags_saturatedcenter', 'r_pixelflags_saturatedcenter',
    'i_pixelflags_saturatedcenter', 'z_pixelflags_saturatedcenter',
    'z_pixelflags_saturatedcenter',
    'g_pixelflags_crcenter', 'r_pixelflags_crcenter',
    'i_pixelflags_crcenter', 'z_pixelflags_crcenter',
    'z_pixelflags_crcenter'
    ]

DR1_CLEAN = [
    'gflags_pixel_edge', 'rflags_pixel_edge', 'iflags_pixel_edge',
    'zflags_pixel_edge', 'yflags_pixel_edge',
    'gflags_pixel_interpolated_center', 'rflags_pixel_interpolated_center',
    'iflags_pixel_interpolated_center', 'zflags_pixel_interpolated_center',
    'yflags_pixel_interpolated_center',
    'gflags_pixel_saturated_center', 'rflags_pixel_saturated_center',
    'iflags_pixel_saturated_center', 'zflags_pixel_saturated_center',
    'yflags_pixel_saturated_center',
    'gflags_pixel_cr_center', 'rflags_pixel_cr_center',
    'iflags_pixel_cr_center', 'zflags_pixel_cr_center',
    'yflags_pixel_cr_center'
]

def basic_meas_photometry(rerun, band):
    """
    Return a dict of column names for basic independent photometric measurements.
    """
    if ('pdr2' in rerun) or ('s18a' in rerun):
        meas_dict = {
            'object_id': 'meas.object_id',
            'ra': 'meas.{}_ra'.format(band),
            'dec': 'meas.{}_dec'.format(band),
            'cmodel_mag': 'meas.{}_cmodel_mag'.format(band),
            'cmodel_mag_err': 'meas.{}_cmodel_magsigma'.format(band),
            'cmodel_ellipse_11': 'meas.{}_cmodel_ellipse_11'.format(band),
            'cmodel_ellipse_22': 'meas.{}_cmodel_ellipse_22'.format(band),
            'cmodel_ellipse_12': 'meas.{}_cmodel_ellipse_12'.format(band),
            'cmodel_exp_mag': 'meas.{}_cmodel_exp_mag'.format(band),
            'cmodel_exp_mag_err': 'meas.{}_cmodel_exp_magsigma'.format(band),
            'cmodel_exp_ellipse_11': 'meas.{}_cmodel_exp_ellipse_11'.format(band),
            'cmodel_exp_ellipse_22': 'meas.{}_cmodel_exp_ellipse_22'.format(band),
            'cmodel_exp_ellipse_12': 'meas.{}_cmodel_exp_ellipse_12'.format(band),
            'cmodel_dev_mag': 'meas.{}_cmodel_dev_mag'.format(band),
            'cmodel_dev_mag_err': 'meas.{}_cmodel_dev_magsigma'.format(band),
            'cmodel_dev_ellipse_11': 'meas.{}_cmodel_dev_ellipse_11'.format(band),
            'cmodel_dev_ellipse_22': 'meas.{}_cmodel_dev_ellipse_22'.format(band),
            'cmodel_dev_ellipse_12': 'meas.{}_cmodel_dev_ellipse_12'.format(band),
            'cmodel_fracdev': 'meas.{}_cmodel_fracdev'.format(band),
            'psf_mag': 'meas2.{}_psfflux_mag'.format(band),
            'psf_mag_err': 'meas2.{}_psfflux_magsigma'.format(band),
        }
    elif 's17a' in rerun:
        meas_dict = {
            'object_id': 'meas.object_id',
            'ra': 'meas.{}_ra'.format(band),
            'dec': 'meas.{}_dec'.format(band),
            'cmodel_mag': 'meas.{}_cmodel_mag'.format(band),
            'cmodel_mag_err': 'meas.{}_cmodel_magsigma'.format(band),
            'cmodel_ellipse_11': 'meas.{}_cmodel_ellipse_11'.format(band),
            'cmodel_ellipse_22': 'meas.{}_cmodel_ellipse_22'.format(band),
            'cmodel_ellipse_12': 'meas.{}_cmodel_ellipse_12'.format(band),
            'cmodel_exp_mag': 'meas.{}_cmodel_exp_mag'.format(band),
            'cmodel_exp_mag_err': 'meas.{}_cmodel_exp_magsigma'.format(band),
            'cmodel_exp_ellipse_11': 'meas.{}_cmodel_exp_ellipse_11'.format(band),
            'cmodel_exp_ellipse_22': 'meas.{}_cmodel_exp_ellipse_22'.format(band),
            'cmodel_exp_ellipse_12': 'meas.{}_cmodel_exp_ellipse_12'.format(band),
            'cmodel_dev_mag': 'meas.{}_cmodel_dev_mag'.format(band),
            'cmodel_dev_mag_err': 'meas.{}_cmodel_dev_magsigma'.format(band),
            'cmodel_dev_ellipse_11': 'meas.{}_cmodel_dev_ellipse_11'.format(band),
            'cmodel_dev_ellipse_22': 'meas.{}_cmodel_dev_ellipse_22'.format(band),
            'cmodel_dev_ellipse_12': 'meas.{}_cmodel_dev_ellipse_12'.format(band),
            'cmodel_fracdev': 'meas.{}_cmodel_fracdev'.format(band),
            'psf_mag': 'meas.{}_psfflux_mag'.format(band),
            'psf_mag_err': 'meas.{}_psfflux_magsigma'.format(band),
        }
    elif ('pdr1' in rerun) or ('s16a' in rerun):
        meas_dict = {
            'object_id': 'meas.object_id',
            'ra': 'meas.{}ra'.format(band),
            'dec': 'meas.{}dec'.format(band),
            'cmodel_mag': 'meas.{}cmodel_mag'.format(band),
            'cmodel_mag_err': 'meas.{}cmodel_mag_err'.format(band),
            'cmodel_ellipse_11': 'meas.{}cmodel_ellipse_11'.format(band),
            'cmodel_ellipse_22': 'meas.{}cmodel_ellipse_22'.format(band),
            'cmodel_ellipse_12': 'meas.{}cmodel_ellipse_12'.format(band),
            'cmodel_exp_mag': 'meas.{}cmodel_exp_mag'.format(band),
            'cmodel_exp_mag_err': 'meas.{}cmodel_exp_mag_err'.format(band),
            'cmodel_exp_ellipse_11': 'meas.{}cmodel_exp_ellipse_11'.format(band),
            'cmodel_exp_ellipse_22': 'meas.{}cmodel_exp_ellipse_22'.format(band),
            'cmodel_exp_ellipse_12': 'meas.{}cmodel_exp_ellipse_12'.format(band),
            'cmodel_dev_mag': 'meas.{}cmodel_dev_mag'.format(band),
            'cmodel_dev_mag_err': 'meas.{}cmodel_dev_mag_err'.format(band),
            'cmodel_dev_ellipse_11': 'meas.{}cmodel_dev_ellipse_11'.format(band),
            'cmodel_dev_ellipse_22': 'meas.{}cmodel_dev_ellipse_22'.format(band),
            'cmodel_dev_ellipse_12': 'meas.{}cmodel_dev_ellipse_12'.format(band),
            'cmodel_fracdev': 'meas.{}cmodel_fracdev'.format(band),
            'psf_mag': 'meas.{}mag_psf'.format(band),
            'psf_mag_err': 'meas.{}mag_psf_err'.format(band),
        }
    else:
        raise NameError("Wrong rerun name")

    return meas_dict

def basic_forced_photometry(rerun, psf=True, cmodel=True, aper=False,
                            shape=False, flux=False, aper_type='3_20'):
    """
    Return a dict of column names for basic photometric measurements.
    """
    # Basic information
    basic_dict = {
        'object_id': 'forced.object_id', 'ra': 'forced.ra', 'dec': 'forced.dec',
        'tract': 'forced.tract', 'patch': 'forced.patch',
        'a_g': 'forced.a_g', 'a_r': 'forced.a_r', 'a_i': 'forced.a_i',
        'a_z': 'forced.a_z', 'a_y': 'forced.a_y',
    }

    if ('pdr2' in rerun) or ('s18a' in rerun):
        # This is for the columns in PDR2 rerun and S18A
        # Flag
        meta_dict = {
            'merge_peak_sky': 'forced.merge_peak_sky',
            'g_inputcount': 'forced.g_inputcount_value',
            'r_inputcount': 'forced.r_inputcount_value',
            'i_inputcount': 'forced.i_inputcount_value',
            'z_inputcount': 'forced.z_inputcount_value',
            'y_inputcount': 'forced.y_inputcount_value',
            'g_flag_edge': 'forced.g_pixelflags_edge',
            'r_flag_edge': 'forced.r_pixelflags_edge',
            'i_flag_edge': 'forced.i_pixelflags_edge',
            'z_flag_edge': 'forced.z_pixelflags_edge',
            'y_flag_edge': 'forced.y_pixelflags_edge',
            'g_flag_saturated': 'forced.g_pixelflags_saturated',
            'r_flag_saturated': 'forced.r_pixelflags_saturated',
            'i_flag_saturated': 'forced.i_pixelflags_saturated',
            'z_flag_saturated': 'forced.z_pixelflags_saturated',
            'y_flag_saturated': 'forced.y_pixelflags_saturated',
            'g_flag_interpolated': 'forced.g_pixelflags_interpolated',
            'r_flag_interpolated': 'forced.r_pixelflags_interpolated',
            'i_flag_interpolated': 'forced.i_pixelflags_interpolated',
            'z_flag_interpolated': 'forced.z_pixelflags_interpolated',
            'y_flag_interpolated': 'forced.y_pixelflags_interpolated',
            'g_flag_saturated_cen': 'forced.g_pixelflags_saturatedcenter',
            'r_flag_saturated_cen': 'forced.r_pixelflags_saturatedcenter',
            'i_flag_saturated_cen': 'forced.i_pixelflags_saturatedcenter',
            'z_flag_saturated_cen': 'forced.z_pixelflags_saturatedcenter',
            'y_flag_saturated_cen': 'forced.y_pixelflags_saturatedcenter',
            'g_flag_interpolated_cen': 'forced.g_pixelflags_interpolatedcenter',
            'r_flag_interpolated_cen': 'forced.r_pixelflags_interpolatedcenter',
            'i_flag_interpolated_cen': 'forced.i_pixelflags_interpolatedcenter',
            'z_flag_interpolated_cen': 'forced.z_pixelflags_interpolatedcenter',
            'y_flag_interpolated_cen': 'forced.y_pixelflags_interpolatedcenter',
            'g_extendedness': 'forced.g_extendedness_value',
            'r_extendedness': 'forced.r_extendedness_value',
            'i_extendedness': 'forced.i_extendedness_value',
            'z_extendedness': 'forced.z_extendedness_value',
            'y_extendedness': 'forced.y_extendedness_value'
        }

        # CModel photometry
        if cmodel:
            cmodel_flag = {
                'g_cmodel_flag': 'forced.g_cmodel_flag',
                'r_cmodel_flag': 'forced.r_cmodel_flag',
                'i_cmodel_flag': 'forced.i_cmodel_flag',
                'z_cmodel_flag': 'forced.z_cmodel_flag',
                'y_cmodel_flag': 'forced.y_cmodel_flag'
            }
            if flux:
                cmodel_dict = {
                    'g_cmodel_flux': 'forced.g_cmodel_flux',
                    'r_cmodel_flux': 'forced.r_cmodel_flux',
                    'i_cmodel_flux': 'forced.i_cmodel_flux',
                    'z_cmodel_flux': 'forced.z_cmodel_flux',
                    'y_cmodel_flux': 'forced.y_cmodel_flux',
                    'g_cmodel_flux_err': 'forced.g_cmodel_fluxsigma',
                    'r_cmodel_flux_err': 'forced.r_cmodel_fluxsigma',
                    'i_cmodel_flux_err': 'forced.i_cmodel_fluxsigma',
                    'z_cmodel_flux_err': 'forced.z_cmodel_fluxsigma',
                    'y_cmodel_flux_err': 'forced.y_cmodel_fluxsigma'
                }
            else:
                cmodel_dict = {
                    'g_cmodel_mag': 'forced.g_cmodel_mag',
                    'r_cmodel_mag': 'forced.r_cmodel_mag',
                    'i_cmodel_mag': 'forced.i_cmodel_mag',
                    'z_cmodel_mag': 'forced.z_cmodel_mag',
                    'y_cmodel_mag': 'forced.y_cmodel_mag',
                    'g_cmodel_mag_err': 'forced.g_cmodel_magsigma',
                    'r_cmodel_mag_err': 'forced.r_cmodel_magsigma',
                    'i_cmodel_mag_err': 'forced.i_cmodel_magsigma',
                    'z_cmodel_mag_err': 'forced.z_cmodel_magsigma',
                    'y_cmodel_mag_err': 'forced.y_cmodel_magsigma'
                }
            # Put the CModel flag
            cmodel_dict.update(cmodel_flag)
        else:
            cmodel_dict = {}

        # PSF photometry
        if psf:
            psf_flag = {
                'g_psf_flag': 'forced2.g_psfflux_flag',
                'r_psf_flag': 'forced2.r_psfflux_flag',
                'i_psf_flag': 'forced2.i_psfflux_flag',
                'z_psf_flag': 'forced2.z_psfflux_flag',
                'y_psf_flag': 'forced2.y_psfflux_flag'
            }
            if flux:
                psf_dict = {
                    'g_psf_flux': 'forced2.g_psfflux_flux',
                    'r_psf_flux': 'forced2.r_psfflux_flux',
                    'i_psf_flux': 'forced2.i_psfflux_flux',
                    'z_psf_flux': 'forced2.z_psfflux_flux',
                    'y_psf_flux': 'forced2.y_psfflux_flux',
                    'g_psf_flux_err': 'forced2.g_psfflux_fluxsigma',
                    'r_psf_flux_err': 'forced2.r_psfflux_fluxsigma',
                    'i_psf_flux_err': 'forced2.i_psfflux_fluxsigma',
                    'z_psf_flux_err': 'forced2.z_psfflux_fluxsigma',
                    'y_psf_flux_err': 'forced2.y_psfflux_fluxsigma'
                }
            else:
                psf_dict = {
                    'g_psf_mag': 'forced2.g_psfflux_mag',
                    'r_psf_mag': 'forced2.r_psfflux_mag',
                    'i_psf_mag': 'forced2.i_psfflux_mag',
                    'z_psf_mag': 'forced2.z_psfflux_mag',
                    'y_psf_mag': 'forced2.y_psfflux_mag',
                    'g_psf_mag_err': 'forced2.g_psfflux_magsigma',
                    'r_psf_mag_err': 'forced2.r_psfflux_magsigma',
                    'i_psf_mag_err': 'forced2.i_psfflux_magsigma',
                    'z_psf_mag_err': 'forced2.z_psfflux_magsigma',
                    'y_psf_mag_err': 'forced2.y_psfflux_magsigma'
                }
            # Put the PSF flag
            psf_dict.update(psf_flag)
        else:
            psf_dict = {}

        # Aperture photometry with matched PSF
        if aper:
            # Flag for aperture photometry
            aper_flag = {
                'g_aper_flag': 'forced4.g_convolvedflux_{}_flag'.format(aper_type),
                'r_aper_flag': 'forced4.r_convolvedflux_{}_flag'.format(aper_type),
                'i_aper_flag': 'forced4.i_convolvedflux_{}_flag'.format(aper_type),
                'z_aper_flag': 'forced4.z_convolvedflux_{}_flag'.format(aper_type),
                'y_aper_flag': 'forced4.y_convolvedflux_{}_flag'.format(aper_type),
            }
            if flux:
                aper_dict = {
                    'g_aper_flux': 'forced4.g_convolvedflux_{}_flux'.format(aper_type),
                    'r_aper_flux': 'forced4.r_convolvedflux_{}_flux'.format(aper_type),
                    'i_aper_flux': 'forced4.i_convolvedflux_{}_flux'.format(aper_type),
                    'z_aper_flux': 'forced4.z_convolvedflux_{}_flux'.format(aper_type),
                    'y_aper_flux': 'forced4.y_convolvedflux_{}_flux'.format(aper_type),
                    'g_aper_flux_err': 'forced4.g_convolvedflux_{}_fluxsigma'.format(aper_type),
                    'r_aper_flux_err': 'forced4.r_convolvedflux_{}_fluxsigma'.format(aper_type),
                    'i_aper_flux_err': 'forced4.i_convolvedflux_{}_fluxsigma'.format(aper_type),
                    'z_aper_flux_err': 'forced4.z_convolvedflux_{}_fluxsigma'.format(aper_type),
                    'y_aper_flux_err': 'forced4.y_convolvedflux_{}_fluxsigma'.format(aper_type),
                }
            else:
                aper_dict = {
                    'g_aper_mag': 'forced4.g_convolvedflux_{}_mag'.format(aper_type),
                    'r_aper_mag': 'forced4.r_convolvedflux_{}_mag'.format(aper_type),
                    'i_aper_mag': 'forced4.i_convolvedflux_{}_mag'.format(aper_type),
                    'z_aper_mag': 'forced4.z_convolvedflux_{}_mag'.format(aper_type),
                    'y_aper_mag': 'forced4.y_convolvedflux_{}_mag'.format(aper_type),
                    'g_aper_mag_err': 'forced4.g_convolvedflux_{}_magsigma'.format(aper_type),
                    'r_aper_mag_err': 'forced4.r_convolvedflux_{}_magsigma'.format(aper_type),
                    'i_aper_mag_err': 'forced4.i_convolvedflux_{}_magsigma'.format(aper_type),
                    'z_aper_mag_err': 'forced4.z_convolvedflux_{}_magsigma'.format(aper_type),
                    'y_aper_mag_err': 'forced4.y_convolvedflux_{}_magsigma'.format(aper_type),
                }
            aper_dict.update(aper_flag)
        else:
            aper_dict = {}

        # Shape of the object
        if shape:
            shape_dict = {
                'g_sdssshape_11': 'forced2.g_sdssshape_shape11',
                'g_sdssshape_22': 'forced2.g_sdssshape_shape22',
                'g_sdssshape_12': 'forced2.g_sdssshape_shape12',
                'r_sdssshape_11': 'forced2.r_sdssshape_shape11',
                'r_sdssshape_22': 'forced2.r_sdssshape_shape22',
                'r_sdssshape_12': 'forced2.r_sdssshape_shape12',
                'i_sdssshape_11': 'forced2.i_sdssshape_shape11',
                'i_sdssshape_22': 'forced2.i_sdssshape_shape22',
                'i_sdssshape_12': 'forced2.i_sdssshape_shape12',
                'z_sdssshape_11': 'forced2.z_sdssshape_shape11',
                'z_sdssshape_22': 'forced2.z_sdssshape_shape22',
                'z_sdssshape_12': 'forced2.z_sdssshape_shape12',
                'y_sdssshape_11': 'forced2.y_sdssshape_shape11',
                'y_sdssshape_22': 'forced2.y_sdssshape_shape22',
                'y_sdssshape_12': 'forced2.y_sdssshape_shape12'
            }
        else:
            shape_dict = {}
    elif 's17a' in rerun:
        # This is for the columns in S17A release
        # Flag
        meta_dict = {
            'merge_peak_sky': 'forced.merge_peak_sky',
            'g_inputcount': 'forced.g_inputcount_value',
            'r_inputcount': 'forced.r_inputcount_value',
            'i_inputcount': 'forced.i_inputcount_value',
            'z_inputcount': 'forced.z_inputcount_value',
            'y_inputcount': 'forced.y_inputcount_value',
            'g_flag_edge': 'forced.g_pixelflags_edge',
            'r_flag_edge': 'forced.r_pixelflags_edge',
            'i_flag_edge': 'forced.i_pixelflags_edge',
            'z_flag_edge': 'forced.z_pixelflags_edge',
            'y_flag_edge': 'forced.y_pixelflags_edge',
            'g_flag_saturated': 'forced.g_pixelflags_saturated',
            'r_flag_saturated': 'forced.r_pixelflags_saturated',
            'i_flag_saturated': 'forced.i_pixelflags_saturated',
            'z_flag_saturated': 'forced.z_pixelflags_saturated',
            'y_flag_saturated': 'forced.y_pixelflags_saturated',
            'g_flag_interpolated': 'forced.g_pixelflags_interpolated',
            'r_flag_interpolated': 'forced.r_pixelflags_interpolated',
            'i_flag_interpolated': 'forced.i_pixelflags_interpolated',
            'z_flag_interpolated': 'forced.z_pixelflags_interpolated',
            'y_flag_interpolated': 'forced.y_pixelflags_interpolated',
            'g_flag_saturated_cen': 'forced.g_pixelflags_saturatedcenter',
            'r_flag_saturated_cen': 'forced.r_pixelflags_saturatedcenter',
            'i_flag_saturated_cen': 'forced.i_pixelflags_saturatedcenter',
            'z_flag_saturated_cen': 'forced.z_pixelflags_saturatedcenter',
            'y_flag_saturated_cen': 'forced.y_pixelflags_saturatedcenter',
            'g_flag_interpolated_cen': 'forced.g_pixelflags_interpolatedcenter',
            'r_flag_interpolated_cen': 'forced.r_pixelflags_interpolatedcenter',
            'i_flag_interpolated_cen': 'forced.i_pixelflags_interpolatedcenter',
            'z_flag_interpolated_cen': 'forced.z_pixelflags_interpolatedcenter',
            'y_flag_interpolated_cen': 'forced.y_pixelflags_interpolatedcenter',
            'g_extendedness': 'forced.g_extendedness_value',
            'r_extendedness': 'forced.r_extendedness_value',
            'i_extendedness': 'forced.i_extendedness_value',
            'z_extendedness': 'forced.z_extendedness_value',
            'y_extendedness': 'forced.y_extendedness_value'
        }

        # CModel photometry
        if cmodel:
            cmodel_flag = {
                'g_cmodel_flag': 'forced.g_cmodel_flag',
                'r_cmodel_flag': 'forced.r_cmodel_flag',
                'i_cmodel_flag': 'forced.i_cmodel_flag',
                'z_cmodel_flag': 'forced.z_cmodel_flag',
                'y_cmodel_flag': 'forced.y_cmodel_flag'
            }
            if flux:
                cmodel_dict = {
                    'g_cmodel_flux': 'forced.g_cmodel_flux',
                    'r_cmodel_flux': 'forced.r_cmodel_flux',
                    'i_cmodel_flux': 'forced.i_cmodel_flux',
                    'z_cmodel_flux': 'forced.z_cmodel_flux',
                    'y_cmodel_flux': 'forced.y_cmodel_flux',
                    'g_cmodel_flux_err': 'forced.g_cmodel_fluxsigma',
                    'r_cmodel_flux_err': 'forced.r_cmodel_fluxsigma',
                    'i_cmodel_flux_err': 'forced.i_cmodel_fluxsigma',
                    'z_cmodel_flux_err': 'forced.z_cmodel_fluxsigma',
                    'y_cmodel_flux_err': 'forced.y_cmodel_fluxsigma'
                }
            else:
                cmodel_dict = {
                    'g_cmodel_mag': 'forced.g_cmodel_mag',
                    'r_cmodel_mag': 'forced.r_cmodel_mag',
                    'i_cmodel_mag': 'forced.i_cmodel_mag',
                    'z_cmodel_mag': 'forced.z_cmodel_mag',
                    'y_cmodel_mag': 'forced.y_cmodel_mag',
                    'g_cmodel_mag_err': 'forced.g_cmodel_magsigma',
                    'r_cmodel_mag_err': 'forced.r_cmodel_magsigma',
                    'i_cmodel_mag_err': 'forced.i_cmodel_magsigma',
                    'z_cmodel_mag_err': 'forced.z_cmodel_magsigma',
                    'y_cmodel_mag_err': 'forced.y_cmodel_magsigma'
                }
            # Put the CModel flag
            cmodel_dict.update(cmodel_flag)
        else:
            cmodel_dict = {}

        # PSF photometry
        if psf:
            psf_flag = {
                'g_psf_flag': 'forced.g_psfflux_flag',
                'r_psf_flag': 'forced.r_psfflux_flag',
                'i_psf_flag': 'forced.i_psfflux_flag',
                'z_psf_flag': 'forced.z_psfflux_flag',
                'y_psf_flag': 'forced.y_psfflux_flag'
            }
            if flux:
                psf_dict = {
                    'g_psf_flux': 'forced.g_psfflux_flux',
                    'r_psf_flux': 'forced.r_psfflux_flux',
                    'i_psf_flux': 'forced.i_psfflux_flux',
                    'z_psf_flux': 'forced.z_psfflux_flux',
                    'y_psf_flux': 'forced.y_psfflux_flux',
                    'g_psf_flux_err': 'forced.g_psfflux_fluxsigma',
                    'r_psf_flux_err': 'forced.r_psfflux_fluxsigma',
                    'i_psf_flux_err': 'forced.i_psfflux_fluxsigma',
                    'z_psf_flux_err': 'forced.z_psfflux_fluxsigma',
                    'y_psf_flux_err': 'forced.y_psfflux_fluxsigma'
                }
            else:
                psf_dict = {
                    'g_psf_mag': 'forced.g_psfflux_mag',
                    'r_psf_mag': 'forced.r_psfflux_mag',
                    'i_psf_mag': 'forced.i_psfflux_mag',
                    'z_psf_mag': 'forced.z_psfflux_mag',
                    'y_psf_mag': 'forced.y_psfflux_mag',
                    'g_psf_mag_err': 'forced.g_psfflux_magsigma',
                    'r_psf_mag_err': 'forced.r_psfflux_magsigma',
                    'i_psf_mag_err': 'forced.i_psfflux_magsigma',
                    'z_psf_mag_err': 'forced.z_psfflux_magsigma',
                    'y_psf_mag_err': 'forced.y_psfflux_magsigma'
                }
            # Put the PSF flag
            psf_dict.update(psf_flag)
        else:
            psf_dict = {}

        # Aperture photometry with matched PSF
        if aper:
            # Flag for aperture photometry
            aper_flag = {
                'g_aper_flag': 'forced3.g_convolvedflux_{}_flag'.format(aper_type),
                'r_aper_flag': 'forced3.r_convolvedflux_{}_flag'.format(aper_type),
                'i_aper_flag': 'forced3.i_convolvedflux_{}_flag'.format(aper_type),
                'z_aper_flag': 'forced3.z_convolvedflux_{}_flag'.format(aper_type),
                'y_aper_flag': 'forced3.y_convolvedflux_{}_flag'.format(aper_type),
            }
            if flux:
                aper_dict = {
                    'g_aper_flux': 'forced3.g_convolvedflux_{}_flux'.format(aper_type),
                    'r_aper_flux': 'forced3.r_convolvedflux_{}_flux'.format(aper_type),
                    'i_aper_flux': 'forced3.i_convolvedflux_{}_flux'.format(aper_type),
                    'z_aper_flux': 'forced3.z_convolvedflux_{}_flux'.format(aper_type),
                    'y_aper_flux': 'forced3.y_convolvedflux_{}_flux'.format(aper_type),
                    'g_aper_flux_err': 'forced3.g_convolvedflux_{}_fluxsigma'.format(aper_type),
                    'r_aper_flux_err': 'forced3.r_convolvedflux_{}_fluxsigma'.format(aper_type),
                    'i_aper_flux_err': 'forced3.i_convolvedflux_{}_fluxsigma'.format(aper_type),
                    'z_aper_flux_err': 'forced3.z_convolvedflux_{}_fluxsigma'.format(aper_type),
                    'y_aper_flux_err': 'forced3.y_convolvedflux_{}_fluxsigma'.format(aper_type),
                }
            else:
                aper_dict = {
                    'g_aper_mag': 'forced3.g_convolvedflux_{}_mag'.format(aper_type),
                    'r_aper_mag': 'forced3.r_convolvedflux_{}_mag'.format(aper_type),
                    'i_aper_mag': 'forced3.i_convolvedflux_{}_mag'.format(aper_type),
                    'z_aper_mag': 'forced3.z_convolvedflux_{}_mag'.format(aper_type),
                    'y_aper_mag': 'forced3.y_convolvedflux_{}_mag'.format(aper_type),
                    'g_aper_mag_err': 'forced3.g_convolvedflux_{}_magsigma'.format(aper_type),
                    'r_aper_mag_err': 'forced3.r_convolvedflux_{}_magsigma'.format(aper_type),
                    'i_aper_mag_err': 'forced3.i_convolvedflux_{}_magsigma'.format(aper_type),
                    'z_aper_mag_err': 'forced3.z_convolvedflux_{}_magsigma'.format(aper_type),
                    'y_aper_mag_err': 'forced3.y_convolvedflux_{}_magsigma'.format(aper_type),
                }
            aper_dict.update(aper_flag)
        else:
            aper_dict = {}

        # Shape of the object
        if shape:
            print("# SDSS Shape is not available for forced photometry in S17A release.")
        shape_dict = {}
    elif ('pdr1' in rerun) or ('s16a' in rerun):
        # This is for the columns in PDR1 and S16A release
        # Flag
        meta_dict = {
            'merge_peak_sky': 'forced.merge_peak_sky',
            'g_inputcount': 'forced.gcountinputs',
            'r_inputcount': 'forced.rcountinputs',
            'i_inputcount': 'forced.icountinputs',
            'z_inputcount': 'forced.zcountinputs',
            'y_inputcount': 'forced.ycountinputs',
            'g_flag_edge': 'forced.gflags_pixel_edge',
            'r_flag_edge': 'forced.rflags_pixel_edge',
            'i_flag_edge': 'forced.iflags_pixel_edge',
            'z_flag_edge': 'forced.zflags_pixel_edge',
            'y_flag_edge': 'forced.yflags_pixel_edge',
            'g_flag_saturated': 'forced.gflags_pixel_saturated_any',
            'r_flag_saturated': 'forced.rflags_pixel_saturated_any',
            'i_flag_saturated': 'forced.iflags_pixel_saturated_any',
            'z_flag_saturated': 'forced.zflags_pixel_saturated_any',
            'y_flag_saturated': 'forced.yflags_pixel_saturated_any',
            'g_flag_interpolated': 'forced.gflags_pixel_interpolated_any',
            'r_flag_interpolated': 'forced.rflags_pixel_interpolated_any',
            'i_flag_interpolated': 'forced.iflags_pixel_interpolated_any',
            'z_flag_interpolated': 'forced.zflags_pixel_interpolated_any',
            'y_flag_interpolated': 'forced.yflags_pixel_interpolated_any',
            'g_flag_saturated_cen': 'forced.gflags_pixel_saturated_center',
            'r_flag_saturated_cen': 'forced.rflags_pixel_saturated_center',
            'i_flag_saturated_cen': 'forced.iflags_pixel_saturated_center',
            'z_flag_saturated_cen': 'forced.zflags_pixel_saturated_center',
            'y_flag_saturated_cen': 'forced.yflags_pixel_saturated_center',
            'g_flag_interpolated_cen': 'forced.gflags_pixel_interpolated_center',
            'r_flag_interpolated_cen': 'forced.rflags_pixel_interpolated_center',
            'i_flag_interpolated_cen': 'forced.iflags_pixel_interpolated_center',
            'z_flag_interpolated_cen': 'forced.zflags_pixel_interpolated_center',
            'y_flag_interpolated_cen': 'forced.yflags_pixel_interpolated_center',
            'g_extendedness': 'forced.gclassification_extendedness',
            'r_extendedness': 'forced.rclassification_extendedness',
            'i_extendedness': 'forced.iclassification_extendedness',
            'z_extendedness': 'forced.zclassification_extendedness',
            'y_extendedness': 'forced.yclassification_extendedness'
        }

        # CModel photometry
        if cmodel:
            cmodel_flag = {
                'g_cmodel_flag': 'forced.gcmodel_flux_flags',
                'r_cmodel_flag': 'forced.rcmodel_flux_flags',
                'i_cmodel_flag': 'forced.icmodel_flux_flags',
                'z_cmodel_flag': 'forced.zcmodel_flux_flags',
                'y_cmodel_flag': 'forced.ycmodel_flux_flags'
            }
            if flux:
                cmodel_dict = {
                    'g_cmodel_flux': 'forced.gcmodel_flux',
                    'r_cmodel_flux': 'forced.rcmodel_flux',
                    'i_cmodel_flux': 'forced.icmodel_flux',
                    'z_cmodel_flux': 'forced.zcmodel_flux',
                    'y_cmodel_flux': 'forced.ycmodel_flux',
                    'g_cmodel_flux_err': 'forced.gcmodel_flux_err',
                    'r_cmodel_flux_err': 'forced.rcmodel_flux_err',
                    'i_cmodel_flux_err': 'forced.icmodel_flux_err',
                    'z_cmodel_flux_err': 'forced.zcmodel_flux_err',
                    'y_cmodel_flux_err': 'forced.ycmodel_flux_err'
                }
            else:
                cmodel_dict = {
                    'g_cmodel_mag': 'forced.gcmodel_mag',
                    'r_cmodel_mag': 'forced.rcmodel_mag',
                    'i_cmodel_mag': 'forced.icmodel_mag',
                    'z_cmodel_mag': 'forced.zcmodel_mag',
                    'y_cmodel_mag': 'forced.ycmodel_mag',
                    'g_cmodel_mag_err': 'forced.gcmodel_mag_err',
                    'r_cmodel_mag_err': 'forced.rcmodel_mag_err',
                    'i_cmodel_mag_err': 'forced.icmodel_mag_err',
                    'z_cmodel_mag_err': 'forced.zcmodel_mag_err',
                    'y_cmodel_mag_err': 'forced.ycmodel_mag_err'
                }
            # Put the CModel flag
            cmodel_dict.update(cmodel_flag)
        else:
            cmodel_dict = {}

        # PSF photometry
        if psf:
            psf_flag = {
                'g_psf_flag': 'forced.gflux_psf_flags',
                'r_psf_flag': 'forced.rflux_psf_flags',
                'i_psf_flag': 'forced.iflux_psf_flags',
                'z_psf_flag': 'forced.zflux_psf_flags',
                'y_psf_flag': 'forced.yflux_psf_flags'
            }
            if flux:
                psf_dict = {
                    'g_psf_flux': 'forced.gflux_psf',
                    'r_psf_flux': 'forced.rflux_psf',
                    'i_psf_flux': 'forced.iflux_psf',
                    'z_psf_flux': 'forced.zflux_psf',
                    'y_psf_flux': 'forced.yflux_psf',
                    'g_psf_flux_err': 'forced.gflux_psf_err',
                    'r_psf_flux_err': 'forced.rflux_psf_err',
                    'i_psf_flux_err': 'forced.iflux_psf_err',
                    'z_psf_flux_err': 'forced.zflux_psf_err',
                    'y_psf_flux_err': 'forced.yflux_psf_err'
                }
            else:
                psf_dict = {
                    'g_psf_mag': 'forced.gmag_psf',
                    'r_psf_mag': 'forced.rmag_psf',
                    'i_psf_mag': 'forced.imag_psf',
                    'z_psf_mag': 'forced.zmag_psf',
                    'y_psf_mag': 'forced.ymag_psf',
                    'g_psf_mag_err': 'forced.gmag_psf_err',
                    'r_psf_mag_err': 'forced.rmag_psf_err',
                    'i_psf_mag_err': 'forced.imag_psf_err',
                    'z_psf_mag_err': 'forced.zmag_psf_err',
                    'y_psf_mag_err': 'forced.ymag_psf_err'
                }
            # Put the PSF flag
            psf_dict.update(psf_flag)
        else:
            psf_dict = {}

        # Aperture photometry with matched PSF
        if aper:
            print("# No PSF matched aperture photometry is available in PDR1")
            aper_dict = {}

        # Shape of the object
        if shape:
            shape_dict = {
                'g_sdssshape_11': 'forced.gshape_sdss_11',
                'g_sdssshape_22': 'forced.gshape_sdss_22',
                'g_sdssshape_12': 'forced.gshape_sdss_12',
                'r_sdssshape_11': 'forced.rshape_sdss_11',
                'r_sdssshape_22': 'forced.rshape_sdss_22',
                'r_sdssshape_12': 'forced.rshape_sdss_12',
                'i_sdssshape_11': 'forced.ishape_sdss_11',
                'i_sdssshape_22': 'forced.ishape_sdss_22',
                'i_sdssshape_12': 'forced.ishape_sdss_12',
                'z_sdssshape_11': 'forced.zshape_sdss_11',
                'z_sdssshape_22': 'forced.zshape_sdss_22',
                'z_sdssshape_12': 'forced.zshape_sdss_12',
                'y_sdssshape_11': 'forced.yshape_sdss_11',
                'y_sdssshape_22': 'forced.yshape_sdss_22',
                'y_sdssshape_12': 'forced.yshape_sdss_12'
            }
        else:
            shape_dict = {}
    else:
        raise NameError("Wrong rerun name")

    # Combine all the dicts together
    basic_dict.update(cmodel_dict)
    basic_dict.update(psf_dict)
    basic_dict.update(aper_dict)
    basic_dict.update(shape_dict)
    basic_dict.update(meta_dict)

    return basic_dict

def column_dict_to_str(columns, add_select=True):
    """
    Convert the dictionary of columns into a SQL search SQL.
    """
    col_str = ', '.join(["{0} AS {1}".format(v, k) for k, v in columns.items()])
    if add_select:
        return 'SELECT ' + col_str
    return col_str

def join_table_by_id(rerun, tables):
    """
    Convert a list of tables into the "FROM" part of SQL search string.
    """
    table_list = ["{0}.{1}".format(rerun, t) for t in tables]
    from_str_0 = 'FROM {0} '
    from_str_1 = 'LEFT JOIN {0} USING (object_id)'

    for ii, tab in enumerate(table_list):
        if ii == 0:
            table_list[ii] = from_str_0.format(tab)
        else:
            table_list[ii] = from_str_1.format(tab)

    from_str = ' '.join(table_list)

    return from_str

def sql_clean_objects(rerun):
    """
    Return a "WHERE" string to select "clean" objects.
    """
    if 'pdr2' in rerun or 's18a' in rerun or 's17a' in rerun:
        return "AND NOT " + " AND NOT ".join(DR2_CLEAN)
    elif 'pdr1' in rerun or 's16a' in rerun:
        return "AND NOT " + " AND NOT ".join(DR1_CLEAN)
    else:
        # TODO: need to support other reruns
        raise NameError("Wrong rerun name")

def box_search(ra1, ra2, dec1, dec2, primary=True, clean=False, dr='pdr2', rerun='pdr2_wide',
               archive=None, psf=True, cmodel=True, aper=False, meas=None,
               shape=False, flux=False, aper_type='3_20', where_list=None):
    """
    Get the SQL template for box search.
    """
    # Login to HSC archive
    if archive is None:
        archive = hsc.Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun

    # The "SELECT" part of the SQL search
    column_dict = basic_forced_photometry(
        rerun, psf=psf, cmodel=cmodel, aper=aper, shape=shape,
        flux=flux, aper_type=aper_type)
    # Only support wide filters for now
    if meas and meas.strip() in 'grizy':
        column_dict.update(basic_meas_photometry(rerun, meas.strip()))
    select_str = column_dict_to_str(column_dict)

    # The "FROM" part of the SQL search
    tables = ['forced']
    if ('pdr2' in rerun) or ('s18a' in rerun):
        if psf or shape:
            tables.append('forced2')
        if aper:
            tables.append('forced4')
        if meas:
            tables.append('meas')
            tables.append('meas2')
    elif 's17a' in rerun:
        if aper:
            tables.append('forced3')
        if meas:
            tables.append('meas')
    elif ('pdr1' in rerun) or ('s16a' in rerun):
        # Only forced catalog is used
        if meas:
            tables.append('meas')
    else:
        raise NameError("Wrong rerun name")
    from_str = join_table_by_id(rerun, tables)

    # The "WHERE" part of the SQL search
    # boxSearch(): Returns true if coord is in a box [ra1, ra2] × [dec1, dec2].
    # (Units are degrees). Note that boxSearch(coord, 350, 370, dec1, dec2) is different
    # from boxSearch(coord, 350, 10, dec1, dec2). In the former, ra ∈ [350, 360] ∪ [0, 10];
    # while in the latter, ra ∈ [10, 350].
    where_str = "WHERE boxSearch(coord, {0}, {1}, {2}, {3})".format(ra1, ra2, dec1, dec2)
    if primary:
        where_str += ' AND isprimary'
    if clean:
        where_str += sql_clean_objects(rerun)
    if where_list:
        where_str += " AND ".join(where_list)

    return ' '.join([select_str, from_str, where_str])

def cone_search(ra, dec, rad, primary=True, clean=False, dr='pdr2', rerun='pdr2_wide',
                archive=None, psf=True, cmodel=True, aper=False, meas=None,
                shape=False, flux=False, aper_type='3_20', where_list=None):
    """
    Get the SQL template for cone search.
    """
    # Login to HSC archive
    if archive is None:
        archive = hsc.Hsc(dr=dr, rerun=rerun)
    else:
        dr = archive.dr
        rerun = archive.rerun

    # The "SELECT" part of the SQL search
    column_dict = basic_forced_photometry(
        rerun, psf=psf, cmodel=cmodel, aper=aper, shape=shape,
        flux=flux, aper_type=aper_type)
    # Only support wide filters for now
    if meas and meas.strip() in 'grizy':
        column_dict.update(basic_meas_photometry(rerun, meas.strip()))
    select_str = column_dict_to_str(column_dict)

    # The "FROM" part of the SQL search
    tables = ['forced']
    if 'pdr2' in rerun or 's18a' in rerun:
        if psf or shape:
            tables.append('forced2')
        if aper:
            tables.append('forced4')
    if 's17a' in rerun:
        if aper:
            tables.append('forced3')
    if 's16a' in rerun or 'pdr1' in rerun:
        # Only forced catalog is used
        pass
    else:
        # TODO: need to support other reruns
        raise NameError("Wrong rerun name")
    from_str = join_table_by_id(rerun, tables)

    # The "WHERE" part of the SQL search
    # coneSearch(): Returns true if coord is within radius arcseconds from (ra, dec).
    # The sky coordinates are in degrees. Radius is in arcseconds.
    where_str = "WHERE coneSearch(coord, {0}, {1}, {2})".format(ra, dec, rad)
    if primary:
        where_str += ' AND isprimary'
    if clean:
        where_str += sql_clean_objects(rerun)
    if where_list:
        where_str += " AND ".join(where_list)

    return ' '.join([select_str, from_str, where_str])
