#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Basic configuration of HSC archive."""

import os
import sys
import getpass
import warnings

from astropy.table import Table

__all__ = ('Field', 'Server', 'Rerun', 'DrException',
           'PDR_URL', 'IDR_URL', 'AVAILABLE_DRS')


# Update: 2019-05-09
PDR_URL = "https://hsc-release.mtk.nao.ac.jp"
IDR_URL = "https://hscdata.mtk.nao.ac.jp"

AVAILABLE_DRS = ['pdr1', 'pdr2', 'dr1', 'dr2', 'dr3']


class DrException(Exception):
    """Class for error related to data release information.
    """
    pass


class Field(object):
    """Class for individual HSC SSP field.

    Examples
    --------

    The examples below illustrate common usage of the `Field` object.

        >>> from unagi.config import Field
        >>> g09 = Field('W_GAMA09H', short='g09',
                           {'filter_available: ['g', 'r']})

    Parameters
    ----------
    name : strin, 'AVAILABLE_DRS'g
        Name of the field
    """
    def __init__(self, name, *initial_dict, **kwargs):
        self.name = str(name).strip()

        for key in kwargs:
            setattr(self, key, kwargs[key])

        for dictionary in initial_dict:
            for key in dictionary:
                setattr(self, key, dictionary[key])


class Server(object):
    """Class for configuration parameters for HSC database.

    Examples
    --------
    The examples below illustrate common usage of the `SspObject` object.

        >>> from unagi.config import Server
        >>> pdr_config = Server(dr='pdr1')

    Parameters
    ----------
    dr: string
        Data release ID
    config_file: str
        Name of the configuration file that contains the username and password.
    """
    def __init__(self, dr='dr3', config_file=None, rerun=None):
        if dr.strip()[0] == 'p':
            """Use the HSC SSP public data release at:
                http://hsc.mtk.nao.ac.jp/ssp/

            More information about the PDR/DAS query:

                https://hsc-release.mtk.nao.ac.jp/das_quarry/manual.html

            So far, only DR1 is available."""
            # Gather login information:
            self._get_credential(pdr=True, config_file=config_file)

            # PDR = Public data release
            self.database = 'PDR'

            if dr == 'pdr1':
                # Outdate warning
                warnings.warn("# PDR1 is outdated and some of the links may be broken." +
                              "Please consider using PDR2 instead!")

                self.data_release = 'pdr1'
                self.deepcoadd = False

                # Time limit for connecting to HSC Server
                self.timeout = 600

                # Useful URLs
                self.base_url = PDR_URL
                # SQL catalog log search search server
                self.cat_url = PDR_URL + "/datasearch/api/catalog_jobs/"
                # Coadd image cutout server
                self.img_url = PDR_URL + "/das_quarry/cgi-bin/quarryImage?"
                # PSF picker server
                self.psf_url = PDR_URL + "/psf/pdr1/cgi/getpsf?"
                # Direct file tree
                self.file_url = PDR_URL + "/archive/filetree/"
                # DAS search server
                self.das_url = PDR_URL + "/das_search/"
                # Ancillary information
                self.map_url = PDR_URL + "/das_search/%s/images/%s/" % (
                    self.data_release, self.data_release)
                self.txt_url = PDR_URL + "/rsrc/patch_info/"

                # Available filters
                """
                Notice:
                -------
                    Only UDEEP+DEEP fields have narrow-band coverage.
                """
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I',
                                    'HSC-Z', 'HSC-Y', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y', 'nb816', 'nb921']

                # Available reruns
                self.rerun_list = ['any', 'pdr1_udeep', 'pdr1_deep', 'pdr1_wide']
                self.rerun_default = 'pdr1_wide'
                self.wide_default = 'pdr1_wide'
                self.deep_default = 'pdr1_deep'
                self.udeep_default = 'pdr1_udeep'

                # Available fields
                """
                Notice:
                -------
                    SSP_AEGIS   : Single Tract for calibration purpose.
                """
                _PDR_UD_COSMOS = {
                    'name': 'UDEEP_COSMOS',
                    'file': 'ud_cosmos',
                    'abbr': 'cos',
                    'type': 'UDEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921'],
                    'field_map': (self.map_url +
                                  'tracts_patches_UD_COSMOS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_UD-COSMOS.txt')
                    }

                _PDR_UD_SXDS = {
                    'name': 'UDEEP_SXDS',
                    'file': 'ud_sxds',
                    'abbr': 'sxd',
                    'type': 'UDEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816'],
                    'field_map': (self.map_url +
                                  'tracts_patches_UD_SXDS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_UD-SXDS.txt')
                    }

                _PDR_D_COSMOS = {
                    'name': 'DEEP_COSMOS',
                    'file': 'd_cosmos',
                    'abbr': 'cos',
                    'type': 'DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921'],
                    'field_map': (self.map_url +
                                  'tracts_patches_D_COSMOS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_D-COSMOS.txt')
                    }

                _PDR_D_DEEP2 = {
                    'name': 'DEEP_DEEP2-3',
                    'file': 'd_deep2',
                    'abbr': 'dep',
                    'type': 'DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816'],
                    'field_map': (self.map_url +
                                  'tracts_patches_D_DEEP2-3_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_D-DEEP2-3.txt')
                    }

                _PDR_D_ELAIS = {
                    'name': 'DEEP_ELAIS-N1',
                    'file': 'd_elais',
                    'abbr': 'ela',
                    'type': 'DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921'],
                    'field_map': (self.map_url +
                                  'tracts_patches_D_ELAIS-N1_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_D-ELAIS-N1.txt')
                    }

                _PDR_D_XMM = {
                    'name': 'DEEP_XMM-LSS',
                    'file': 'd_xmm',
                    'abbr': 'xmm',
                    'type': 'DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_D_XMM-LSS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_D-XMM-LSS.txt')
                    }

                _PDR_W_XMM = {
                    'name': 'WIDE_XMM-LSS',
                    'file': 'w_xmm',
                    'abbr': 'xmm',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'W-XMMLSS_i_area.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-XMM.txt')
                    }

                _PDR_W_GAMA09 = {
                    'name': 'WIDE_GAMA09H',
                    'file': 'w_gama09',
                    'abbr': 'g09',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'W-GAMA09H_i_area.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-GAMA09H.txt')
                    }

                _PDR_W_GAMA15 = {
                    'name': 'WIDE_GAMA15H',
                    'file': 'w_gama15',
                    'abbr': 'g15',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'W-GAMA15H_i_area.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-GAMA15H.txt')
                    }

                _PDR_W_WIDE12 = {
                    'name': 'WIDE_WIDE12H',
                    'file': 'w_wide12',
                    'abbr': 'w12',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'W-WIDE12H_i_area.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-GAMA15H.txt')
                    }

                _PDR_W_HECTOMAP = {
                    'name': 'WIDE_HECTOMAP',
                    'file': 'w_hectomap',
                    'abbr': 'hec',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'W-HECTOMAP_i_area.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-HECTOMAP.txt')
                    }

                _PDR_W_VVDS = {
                    'name': 'WIDE_VVDS',
                    'file': 'w_vvds',
                    'abbr': 'vvd',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'W-VVDS_i_area.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-VVDS.txt')
                    }

                _PDR_AEGIS = {
                    'name': 'SSP_AEGIS',
                    'file': 'aegis',
                    'abbr': 'aeg',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'AEGIS_i_area.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-AEGIS.txt')
                    }

                self.fields = [_PDR_UD_COSMOS, _PDR_UD_SXDS,
                               _PDR_D_COSMOS, _PDR_D_DEEP2, _PDR_D_ELAIS,
                               _PDR_D_XMM,
                               _PDR_W_XMM, _PDR_W_GAMA09, _PDR_W_GAMA15,
                               _PDR_W_WIDE12, _PDR_W_HECTOMAP,
                               _PDR_W_VVDS, _PDR_AEGIS]
            elif dr == 'pdr2':
                self.data_release = 'pdr2'
                self.deepcoadd = False

                # Time limit for connecting to HSC Server
                self.timeout = 600

                # Useful URLs
                self.base_url = PDR_URL
                # SQL catalog log search search server
                self.cat_url = PDR_URL + "/datasearch/api/catalog_jobs/"
                # Coadd image cutout server
                self.img_url = PDR_URL + "/das_cutout/pdr2/cgi-bin/cutout?"
                # Coadd patch image url
                self.patch_url = PDR_URL + "/archive/filetree/pdr2_wide/deepCoadd-results/"
                # PSF picker server
                self.psf_url = PDR_URL + "/psf/pdr2/cgi/getpsf?"
                # Direct file tree
                self.file_url = PDR_URL + "/archive/filetree/"
                # DAS search server
                self.das_url = PDR_URL + "/das_search/"
                # Ancillary information
                self.map_url = PDR_URL + "/rsrc/pdr2/koike/survey-area/fig/"
                self.txt_url = PDR_URL + "/rsrc/pdr2/koike/survey-area/info/"

                # Available filters
                """
                Notice:
                -------
                    Only UDEEP+DEEP fields have narrow-band coverage.
                """
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I',
                                    'HSC-Z', 'HSC-Y', 'NB0387', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y', 'nb387', 'nb816', 'nb921']

                # Available reruns
                self.rerun_list = ['any', 'pdr2_dud', 'pdr2_wide',
                                   'pdr2_cosmos_wide_depth_best',
                                   'pdr2_cosmos_wide_depth_median',
                                   'pdr2_cosmos_wide_depth_worst']
                self.rerun_default = 'pdr2_wide'
                self.wide_default = 'pdr2_wide'
                self.deep_default = 'pdr2_dud'
                self.udeep_default = 'pdr2_dud'

                # Available fields
                """
                Notice:
                -------
                    From PDR2, the DEEP and UDEEP fields are released together.

                    DUD_XMM-LSS : XMM-LSS + SXDS
                    WIDE_W02 : The old XMM-LSS wide field
                    WIDE_W03 : The old GAMA09H wide field
                    WIDE_W04 : The old WIDE12+GAMA15H wide field
                    WIDE_W05 : The old VVDS wide field
                    WIDE_W06 : The old HECTOMAP wide field
                    WIDE_W07 : The old AEGIS field
                """
                _PDR_DUD_COSMOS = {
                    'name': 'DUD_COSMOS',
                    'file': 'dud_cosmos',
                    'abbr': 'cos',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_COSMOS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-COSMOS.txt')
                    }

                _PDR_DUD_DEEP2 = {
                    'name': 'DUD_DEEP2-3',
                    'file': 'dud_deep2',
                    'abbr': 'dep',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_DEEP2-3_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-DEEP2-3.txt')
                    }

                _PDR_DUD_ELAIS = {
                    'name': 'DUD_ELAIS-N1',
                    'file': 'dud_elais',
                    'abbr': 'ela',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_ELAIS-N1_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-ELAIS-N1.txt')
                    }

                _PDR_DUD_XMM = {
                    'name': 'DUD_XMM-LSS',
                    'file': 'dud_xmm',
                    'abbr': 'xmm',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_XMM-LSS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-XMM-LSS.txt')
                    }

                _PDR_W_WIDE01 = {
                    'name': 'WIDE_WIDE01',
                    'file': 'w_wide01',
                    'abbr': 'w01',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_w01_HSC-R.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w01.txt')
                    }

                _PDR_W_WIDE02 = {
                    'name': 'WIDE_WIDE02',
                    'file': 'w_wide02',
                    'abbr': 'w02',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_w02_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w02.txt')
                    }

                _PDR_W_WIDE03 = {
                    'name': 'WIDE_WIDE03',
                    'file': 'w_wide03',
                    'abbr': 'w03',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_w03_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w03.txt')
                    }

                _PDR_W_WIDE04 = {
                    'name': 'WIDE_WIDE04',
                    'file': 'w_wide04',
                    'abbr': 'w04',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_w04_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w04.txt')
                    }

                _PDR_W_WIDE05 = {
                    'name': 'WIDE_WIDE05',
                    'file': 'w_wide05',
                    'abbr': 'w05',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_w05_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w05.txt')
                    }

                _PDR_W_WIDE06 = {
                    'name': 'WIDE_WIDE06',
                    'file': 'w_wide06',
                    'abbr': 'w06',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_w06_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w06.txt')
                    }

                _PDR_W_WIDE07 = {
                    'name': 'WIDE_WIDE07',
                    'file': 'w_wide07',
                    'abbr': 'w07',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_AEGIS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-AEGIS.txt')
                    }

                self.fields = [_PDR_DUD_COSMOS, _PDR_DUD_DEEP2,
                               _PDR_DUD_ELAIS, _PDR_DUD_XMM,
                               _PDR_W_WIDE01, _PDR_W_WIDE02, _PDR_W_WIDE03,
                               _PDR_W_WIDE04, _PDR_W_WIDE05, _PDR_W_WIDE06,
                               _PDR_W_WIDE07]
            else:
                raise DrException("!! Wrong information about data release !!")

            self.field_table = Table(rows=self.fields)

            self.field_name = list(self.field_table['name'].data.astype('str'))
        else:
            """Use the HSC SSP internal data release at:

            More information about the IDR/DAS query:

                https://hscdata.mtk.nao.ac.jp/das_quarry/dr2/manual.html
                https://hscdata.mtk.nao.ac.jp/das_quarry/dr3/manual.html

            So far, DR2 and DR3 are available. DR1 has been removed from NAOJ to save space.

            Information about the S17A data release:

                https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr2/s17a/doc/products_status.html

            Information about the S18A data release:

                https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr2/s18a/doc/products_status.html

            Information about the S19A data release:

                https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr3/s19a/doc/products_status.html

            Information about the S20A data release:

                https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr3/s20a/doc/products_status.html

            """
            # Gather login information:
            self._get_credential(pdr=False)

            # IDR = Internal data release
            self.database = 'IDR'

            # Time limit for connecting to HSC Server
            self.timeout = 300

            if dr == 'dr1':
                raise ValueError("# DR1 has become unavailable!")
                # TODO: Remove later
                self.data_release = 'dr1'
                self.deepcoadd = False

                # Useful URLs
                self.base_url = IDR_URL
                # SQL catalog log search search server
                self.cat_url = IDR_URL + "/datasearch/api/catalog_jobs/"
                # Coadd image cutout server
                self.img_url = IDR_URL + "/das_quarry/dr1/cgi-bin/quarryImage?"
                # PSF picker server
                self.psf_url = IDR_URL + "/psf/6/cgi/getpsf?"
                # DAS querry server
                self.das_url = IDR_URL + "/das_console/dr1/"
                # Direct file tree
                # Ancillary information
                self.map_url = IDR_URL + "/hsc_ssp/dr1/s16a/doc/fig/"
                self.txt_url = IDR_URL + "/hsc_ssp/dr1/s16a/doc/info/"

                # Available filters
                """
                Notice:
                -------
                    Only UDEEP+DEEP fields have narrow-band coverage.
                """
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I',
                                    'HSC-Z', 'HSC-Y', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y', 'nb816', 'nb921']

                # Available reruns
                """
                We ignore the early data release and the s15a release.
                All their data are covered by s15b and s16a release, and
                the data quality has been significantly improved
                """
                self.rerun_list = ['any',
                                   's15b_udeep', 's15b_deep', 's15b_wide',
                                   's16a_udeep', 's16a_deep', 's16a_wide',
                                   's16a_wide2', 's16a_chorus_cosmos',
                                   's16a_chorus_cosmos_nb', 's16a_chorus_sxds',
                                   's16a_clauds_deep', 's16a_clauds_udeep_deep_depth']
                self.rerun_default = 's16a_wide2'
                self.wide_default = 's16a_wide2'
                self.deep_default = 's16a_deep'
                self.udeep_default = 's16a_udeep'

                # Available fields
                """
                Notice:
                -------
                    SSP_AEGIS   : Single Tract for calibration purpose.
                    WIDE_WIDE01 : Only has g and r-band data.
                """
                _IDR_UD_COSMOS = {
                    'name': 'UDEEP_COSMOS',
                    'file': 'ud_cosmos',
                    'abbr': 'cos',
                    'type': 'UDEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816'],
                    'field_map': (self.map_url +
                                  'tracts_patches_UD_COSMOS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_UD-COSMOS.txt')
                    }

                _IDR_UD_SXDS = {
                    'name': 'UDEEP_SXDS',
                    'file': 'ud_sxds',
                    'abbr': 'sxd',
                    'type': 'UDEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816'],
                    'field_map': (self.map_url +
                                  'tracts_patches_UD_SXDS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_UD-SXDS.txt')
                    }

                _IDR_D_COSMOS = {
                    'name': 'DEEP_COSMOS',
                    'file': 'd_cosmos',
                    'abbr': 'cos',
                    'type': 'DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921'],
                    'field_map': (self.map_url +
                                  'tracts_patches_D_COSMOS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_D-COSMOS.txt')
                    }

                _IDR_D_DEEP2 = {
                    'name': 'DEEP_DEEP2-3',
                    'file': 'd_deep2',
                    'abbr': 'dep',
                    'type': 'DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816'],
                    'field_map': (self.map_url +
                                  'tracts_patches_D_DEEP2-3_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_D-DEEP2-3.txt')
                    }

                _IDR_D_ELAIS = {
                    'name': 'DEEP_ELAIS-N1',
                    'file': 'd_elais',
                    'abbr': 'ela',
                    'type': 'DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816'],
                    'field_map': (self.map_url +
                                  'tracts_patches_D_ELAIS-N1_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_D-ELAIS-N1.txt')
                    }

                _IDR_D_XMM = {
                    'name': 'DEEP_XMM-LSS',
                    'file': 'd_xmm',
                    'abbr': 'xmm',
                    'type': 'DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_D_XMM-LSS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_D-XMM-LSS.txt')
                    }

                _IDR_W_WIDE01 = {
                    'name': 'WIDE_WIDE01H',
                    'file': 'w_wide01',
                    'abbr': 'w01',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_WIDE01H_HSC-R.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-WIDE01H.txt')
                    }

                _IDR_W_XMM = {
                    'name': 'WIDE_XMM-LSS',
                    'file': 'w_xmm',
                    'abbr': 'xmm',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_XMM_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-XMM.txt')
                    }

                _IDR_W_GAMA09 = {
                    'name': 'WIDE_GAMA09H',
                    'file': 'w_gama09',
                    'abbr': 'g09',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_GAMA09H_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-GAMA09H.txt')
                    }

                _IDR_W_GAMA15 = {
                    'name': 'WIDE_GAMA15H',
                    'file': 'w_gama15',
                    'abbr': 'g15',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_GAMA15H_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-GAMA15H.txt')
                    }

                _IDR_W_WIDE12 = {
                    'name': 'WIDE_WIDE12H',
                    'file': 'w_wide12',
                    'abbr': 'w12',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_WIDE12H_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-WIDE12H.txt')
                    }

                _IDR_W_HECTOMAP = {
                    'name': 'WIDE_HECTOMAP',
                    'file': 'w_hectomap',
                    'abbr': 'hec',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_HECTOMAP_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-HECTOMAP.txt')
                    }

                _IDR_W_VVDS = {
                    'name': 'WIDE_VVDS',
                    'file': 'w_vvds',
                    'abbr': 'vvd',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_VVDS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-VVDS.txt')
                    }

                _IDR_AEGIS = {
                    'name': 'SSP_AEGIS',
                    'file': 'aegis',
                    'abbr': 'aeg',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_AEGIS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-AEGIS.txt')
                    }

                self.fields = [_IDR_UD_COSMOS, _IDR_UD_SXDS,
                               _IDR_D_COSMOS, _IDR_D_DEEP2, _IDR_D_ELAIS,
                               _IDR_D_XMM, _IDR_W_WIDE01,
                               _IDR_W_XMM, _IDR_W_GAMA09, _IDR_W_GAMA15,
                               _IDR_W_WIDE12, _IDR_W_HECTOMAP,
                               _IDR_W_VVDS, _IDR_AEGIS]

            elif dr == 'dr2':
                warnings.warn("# DR3 has become available now!")
                self.data_release = 'dr2'
                self.deepcoadd = False

                # Useful URLs
                self.base_url = IDR_URL
                # SQL catalog log search search server
                self.cat_url = IDR_URL + "/datasearch/api/catalog_jobs/"
                # Coadd image cutout server
                if 's17a' in rerun.strip():
                    self.img_url = IDR_URL + "/das_quarry/dr2/cgi-bin/quarryImage?"
                else:
                    self.img_url = IDR_URL + "/das_quarry/dr2.1/cgi-bin/cutout?"
                # Coadd patch image url
                self.patch_url = IDR_URL + "/hsc_ssp/dr2/s18a/data/s18a_wide/deepCoadd-results/"
                # PSF picker server
                self.psf_url = IDR_URL + "/psf/6/cgi/getpsf?"
                # Direct file tree
                self.file_url = IDR_URL + "/hsc_ssp/dr2/"
                # DAS search server
                # TODO: Not sure this is going to work
                self.das_url = IDR_URL + "/das_console/dr2.1/"
                # Ancillary information
                self.map_url = IDR_URL + "/hsc_ssp/dr2/s18a/doc/fig/"
                self.txt_url = IDR_URL + "/hsc_ssp/dr2/s18a/doc/info/"

                # Available filters
                """
                Notice:
                -------
                    Only UDEEP+DEEP fields have narrow-band coverage.
                """
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y',
                                    'NB0387', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y', 'nb0387', 'nb816', 'nb921']

                # Available reruns
                """
                2018-11: S17A and S18A releases are available
                """
                self.rerun_list = ['any', 's17a_dud', 's17a_wide',
                                   's17a_dud_wide_depth_best',
                                   's17a_dud_wide_depth_median',
                                   's17a_dud_wide_depth_median',
                                   's18a_wide', 's18a_dud',
                                   's18a_dud_wide_depth_best',
                                   's18a_dud_wide_depth_median',
                                   's18a_dud_wide_depth_median',
                                   's17a_chorus', 's18a_chorus'
                                   ]
                self.rerun_default = 's18a_wide'
                self.wide_default = 's18a_wide'
                self.deep_default = 's18a_dud'
                self.udeep_default = 's18a_dud'

                # Available fields
                """
                Notice:
                -------
                    From S17A and S18A, the DEEP and UDEEP fields are released together.

                    DUD_XMM-LSS : XMM-LSS + SXDS
                    WIDE_W02 : The old XMM-LSS wide field
                    WIDE_W03 : The old GAMA09H wide field
                    WIDE_W04 : The old WIDE12+GAMA15H wide field
                    WIDE_W05 : The old VVDS wide field
                    WIDE_W06 : The old HECTOMAP wide field
                    WIDE_W07 : The old AEGIS field
                """
                _IDR_DUD_COSMOS = {
                    'name': 'DUD_COSMOS',
                    'file': 'dud_cosmos',
                    'abbr': 'cos',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_COSMOS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-COSMOS.txt')
                    }

                _IDR_DUD_DEEP2 = {
                    'name': 'DUD_DEEP2-3',
                    'file': 'dud_deep2',
                    'abbr': 'dep',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_DEEP2-3_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-DEEP2-3.txt')
                    }

                _IDR_DUD_ELAIS = {
                    'name': 'DUD_ELAIS-N1',
                    'file': 'dud_elais',
                    'abbr': 'ela',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_ELAIS-N1_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-ELAIS-N1.txt')
                    }

                _IDR_DUD_XMM = {
                    'name': 'DUD_XMM-LSS',
                    'file': 'dud_xmm',
                    'abbr': 'xmm',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_XMM-LSS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-XMM-LSS.txt')
                    }

                _IDR_W_WIDE01 = {
                    'name': 'WIDE_WIDE01',
                    'file': 'w_wide01',
                    'abbr': 'w01',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W-w01_HSC-R.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w01.txt')
                    }

                _IDR_W_WIDE02 = {
                    'name': 'WIDE_WIDE02',
                    'file': 'w_wide02',
                    'abbr': 'w02',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W-w02_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w02.txt')
                    }

                _IDR_W_WIDE03 = {
                    'name': 'WIDE_WIDE03',
                    'file': 'w_wide03',
                    'abbr': 'w03',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W-w03_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w03.txt')
                    }

                _IDR_W_WIDE04 = {
                    'name': 'WIDE_WIDE04',
                    'file': 'w_wide04',
                    'abbr': 'w04',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W-w04_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w04.txt')
                    }

                _IDR_W_WIDE05 = {
                    'name': 'WIDE_WIDE05',
                    'file': 'w_wide05',
                    'abbr': 'w05',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W-w05_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w05.txt')
                    }

                _IDR_W_WIDE06 = {
                    'name': 'WIDE_WIDE06',
                    'file': 'w_wide06',
                    'abbr': 'w06',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W-w06_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-w06.txt')
                    }

                _IDR_W_WIDE07 = {
                    'name': 'WIDE_WIDE07',
                    'file': 'w_wide07',
                    'abbr': 'w07',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W-AEGIS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-AEGIS.txt')
                    }

                self.fields = [_IDR_DUD_COSMOS, _IDR_DUD_DEEP2,
                               _IDR_DUD_ELAIS, _IDR_DUD_XMM,
                               _IDR_W_WIDE01, _IDR_W_WIDE02, _IDR_W_WIDE03,
                               _IDR_W_WIDE04, _IDR_W_WIDE05, _IDR_W_WIDE06,
                               _IDR_W_WIDE07]

            elif dr == 'dr3':
                self.data_release = 'dr3'
                self.deepcoadd = True

                # Useful URLs
                self.base_url = IDR_URL
                # SQL catalog log search search server
                self.cat_url = IDR_URL + "/datasearch/api/catalog_jobs/"
                # Coadd image cutout server
                self.img_url = IDR_URL + "/das_quarry/dr3/cgi-bin/cutout?"
                # Coadd patch image url
                # TODO: make this not rerun specific
                self.patch_url = IDR_URL + "/hsc_ssp/dr3/s20a/data/s20a_wide/deepCoadd-results/"
                # PSF picker server
                self.psf_url = IDR_URL + "/psf/8/cgi/getpsf?"
                # Direct file tree
                self.file_url = IDR_URL + "/hsc_ssp/dr3/"
                # DAS search server
                # See: https://hscdata.mtk.nao.ac.jp/das_console/dr3.1/usage.html
                self.das_url = IDR_URL + "/das_console/dr3.1/"
                # Ancillary information
                # TODO: make this not rerun specific
                self.map_url = IDR_URL + "/hsc_ssp/dr3/s18a/doc/fig/"
                self.txt_url = IDR_URL + "/hsc_ssp/dr3/s18a/doc/info/"
                # Quality assesment dir: https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr3/s19a/doc/qa/
                # Stellar sequence dir: https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr3/s19a/doc/ss/

                # Available filters
                """
                Notice:
                -------
                    Only UDEEP+DEEP fields have narrow-band coverage.
                """
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y',
                                    'NB0387', 'NB0816', 'NB0921', 'NB1010']
                self.filter_list_short = [
                    'g', 'r', 'i', 'z', 'y', 'nb0387', 'nb816', 'nb921', 'nb1010']

                # Available reruns
                """
                2019-09: S19A releases are available
                2020-08: S20A releases are available
                """
                self.rerun_list = ['any', 's19a_dud', 's19a_wide',
                                   's20a_wide', 's20a_dud', 's19a_dud_wide_depth_worst',
                                   's19a_dud_wide_depth_best', 's19a_dud_wide_depth_median',
                                   's20a_dud_wide_depth_worst',
                                   's20a_dud_wide_depth_best', 's20a_dud_wide_depth_median']
                self.rerun_default = 's20a_wide'
                self.wide_default = 's20a_wide'
                self.deep_default = 's20a_dud'
                self.udeep_default = 's20a_dud'

                """
                Notice:
                -------
                    From S19A and S20A, the WIDE fields are separated into Spring, Autumn, and
                    HectoMap field, also the AEGIS region.

                    Notice that the Spring and Autumn fields are not continue in every filter.
                """
                _IDR_DUD_COSMOS = {
                    'name': 'DUD_COSMOS',
                    'file': 'dud_cosmos',
                    'abbr': 'cos',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387', 'NB1010'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_COSMOS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-COSMOS.txt')
                    }

                _IDR_DUD_DEEP2 = {
                    'name': 'DUD_DEEP2-3',
                    'file': 'dud_deep2',
                    'abbr': 'dep',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_DEEP2-3_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-DEEP2-3.txt')
                    }

                _IDR_DUD_ELAIS = {
                    'name': 'DUD_ELAIS-N1',
                    'file': 'dud_elais',
                    'abbr': 'ela',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_ELAIS-N1_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-ELAIS-N1.txt')
                    }

                _IDR_DUD_XMM = {
                    'name': 'DUD_XMM-LSS',
                    'file': 'dud_xmm',
                    'abbr': 'xmm',
                    'type': 'UDEEP_DEEP',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y', 'NB0921',
                                         'NB0816', 'NB0387', 'NB1010'],
                    'field_map': (self.map_url +
                                  'tracts_patches_DUD_XMM-LSS_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_DUD-XMM-LSS.txt')
                    }

                _IDR_W_SPRING = {
                    'name': 'WIDE_SPRING',
                    'file': 'w_spring',
                    'abbr': 'spring',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_spring_HSC-I.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-spring.txt')
                    }

                _IDR_W_AUTUMN = {
                    'name': 'WIDE_AUTUMN',
                    'file': 'w_autumn',
                    'abbr': 'autumn',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_autumn_HSC-Z.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-autumn.txt')
                    }

                _IDR_W_HECTOMAP = {
                    'name': 'WIDE_HECTOMAP',
                    'file': 'w_hectomap',
                    'abbr': 'hectomap',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_hectomap_HSC-Z.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-hectomap.txt')
                    }

                _IDR_W_AEGIS = {
                    'name': 'WIDE_AEGIS',
                    'file': 'w_aegis',
                    'abbr': 'aegis',
                    'type': 'WIDE',
                    'filter_available': ['HSC-g', 'HSC-r', 'HSC-i',
                                         'HSC-z', 'HSC-y'],
                    'field_map': (self.map_url +
                                  'tracts_patches_W_AEGIS_HSC-Z.png'),
                    'patch_info': (self.txt_url +
                                   'tracts_patches_W-AEGIS.txt')
                    }

                self.fields = [_IDR_DUD_COSMOS, _IDR_DUD_DEEP2,
                               _IDR_DUD_ELAIS, _IDR_DUD_XMM, _IDR_W_AEGIS,
                               _IDR_W_SPRING, _IDR_W_AUTUMN, _IDR_W_HECTOMAP,]
            else:
                raise DrException("!! Wrong information about data release !!")

            self.field_table = Table(rows=self.fields)

            self.field_name = self.field_table['name'].data.astype('str')

    def _get_credential(self, pdr=False, config_file=None):
        """Get the username and password for HSC SSP database.
        """
        if config_file is not None:
            # Import HSC username and password from a file
            config = Table.read(config_file, format='ascii.no_header')['col1']
            self._username = config[0]
            self._password = config[1]
        else:
            # Get the credential information from environmental variables or user input
            if not pdr:
                try:
                    self._username = os.environ['SSP_IDR_USR']
                    self._password = os.environ['SSP_IDR_PWD']
                except KeyError:
                    get_input = input
                    if sys.version_info[:2] <= (2, 7):
                        get_input = input
                    self._username = get_input("Internal Data Release Username : ")
                    self._password = getpass.getpass("Password : ")
            else:
                try:
                    self._username = os.environ['SSP_PDR_USR']
                    self._password = os.environ['SSP_PDR_PWD']
                except KeyError:
                    get_input = input
                    if sys.version_info[:2] <= (2, 7):
                        get_input = input
                    self._username = get_input("Public Data Release Username : ")
                    self._password = getpass.getpass("Password : ")


class Rerun(Server):
    """Class for rerun in HSC data release.

    Examples
    --------

    The examples below illustrate common usage of the `Rerun` object.

        >>> from unagi.config import Rerun
        >>> s16a_wide2 = Rerun('s16a_wide2')

    Parameters
    ----------

    rerun_name : string
        Name of the rerun
    """
    def __init__(self, rerun='s18a_wide', dr='dr2', config_file=None):
        super(Rerun, self).__init__(
            dr=dr, config_file=config_file, rerun=rerun)

        if str(rerun).strip() not in self.rerun_list:
            raise DrException("!! Wrong rerun !!")
        else:
            self.rerun = rerun
            delattr(self, 'rerun_list')
            delattr(self, 'rerun_default')
            delattr(self, 'udeep_default')
            delattr(self, 'deep_default')
            delattr(self, 'wide_default')
