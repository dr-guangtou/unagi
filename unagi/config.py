#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from astropy.table import Table

import os
import sys
import warnings

__all__ = ['HscField', 'HscDasConfig']


# Update: 2017-11-22
PDR_URL = "https://hsc-release.mtk.nao.ac.jp"
IDR_URL = "https://hscdata.mtk.nao.ac.jp"


class DrException(Exception):
    """Class for error related to data release information.
    """
    pass


class HscField(object):
    """Class for individual HSC SSP field.

    Examples
    --------

    The examples below illustrate common usage of the `HscField` object.

        >>> from unagi.config import HscField
        >>> g09 = HscField('W_GAMA09H', short='g09',
                           {'filter_available: ['g', 'r']})

    Parameters
    ----------
    name : string
        Name of the field
    """
    def __init__(self, name, *initial_dict, **kwargs):
        self.name = str(name).strip()

        for key in kwargs:
            setattr(self, key, kwargs[key])

        for dictionary in initial_dict:
            for key in dictionary:
                setattr(self, key, dictionary[key])


class HscDasConfig(object):
    """Class for configuration parameters for HSC database.

    Examples
    --------
    The examples below illustrate common usage of the `HscObject` object.

        >>> from unagi.config import HscDasConfig
        >>> pdr_config = HscDasConfig(pdr=True)

    Parameters
    ----------
    dr : string
        Data release.
    pdr : boolen, optional
        Use public data release or not (Default = False).
    """
    def __init__(self, dr='dr1', pdr=False):
        if pdr:
            """Use the HSC public data release at:
                http://hsc.mtk.nao.ac.jp/ssp/

            More information about the PDR/DAS query:

                https://hsc-release.mtk.nao.ac.jp/das_quarry/manual.html

            So far, only DR1 is available.
            """
            # Gather login information:
            self._get_credential(pdr=True)

            # Data release list
            self.dr_list = ['dr1']

            # PDR = Public data release
            self.database = 'PDR'

            if dr == 'dr1':

                # Useful URLs
                self.das_url = PDR_URL + "/das_quarry/cgi-bin/"
                self.psf_url = PDR_URL + "/psf/pdr1/cgi/"
                self.map_url = PDR_URL + '/das_search/pdr1/images/pdr1/'
                self.txt_url = PDR_URL + '/rsrc/patch_info/'

                # Available filters
                """
                Notice:
                -------
                    Only UDEEP+DEEP fields have narrow-band coverage.
                """
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I',
                                    'HSC-Z', 'HSC-Y', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y',
                                          'nb816', 'nb921']

                # Available reruns
                self.rerun_list = ['any', 'pdr1_udeep', 'pdr1_deep',
                                   'pdr1_wide']

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

                self.field = [_PDR_UD_COSMOS, _PDR_UD_SXDS,
                              _PDR_D_COSMOS, _PDR_D_DEEP2, _PDR_D_ELAIS,
                              _PDR_D_XMM,
                              _PDR_W_XMM, _PDR_W_GAMA09, _PDR_W_GAMA15,
                              _PDR_W_WIDE12, _PDR_W_HECTOMAP,
                              _PDR_W_VVDS, _PDR_AEGIS]

                self.field_table = Table(rows=self.fields)

                self.field_name = self.field_table['name'].data.astype('str')

            else:
                raise DrException("!! Wrong information about data release !!")
        else:
            """Use the HSC internal data release at:

            More information about the IDR/DAS query:

                https://hscdata.mtk.nao.ac.jp/das_quarry/dr1/manual.html

            So far, DR1 and DR2 are available

            Information about S15b data release:

                https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr1/s15b/doc/S15B/products_status.html

            Information about S16a data release:

                https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr1/s16a/doc/products_status.html
            """
            # Gather login information:
            self._get_credential(pdr=False)

            # Data release list
            self.dr_list = ['dr1', 'dr2']

            # PDR = Public data release
            self.database = 'IDR'

            if dr is 'dr1':
                # Useful URLs
                self.das_url = PDR_URL + "/das_quarry/dr1/cgi-bin/"
                self.psf_url = PDR_URL + "/psf/4/cgi/"

                # Available filters
                """
                Notice:
                -------
                    Only UDEEP+DEEP fields have narrow-band coverage.
                """
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I',
                                    'HSC-Z', 'HSC-Y', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y',
                                          'nb816', 'nb921']

                # Available reruns
                """
                We ignore the early data release and the s15a release.
                All their data are covered by s15b and s16a release, and
                the data quality has been significantly improved
                """
                self.rerun_list = ['any',
                                   's15b_udeep', 's15b_deep', 's15b_wide',
                                   's16a_udeep', 's16a_deep', 's16a_wide',
                                   's16a_wide2']

                # Available fields
                """
                Notice:
                -------
                    SSP_AEGIS   : Single Tract for calibration purpose.
                    WIDE_WIDE01 : Only has g and r-band data.
                """
                self.field_list = ['UDEEP_COSMOS', 'UDEEP_SXDF',
                                   'DEEP_COSMOS', 'DEEP_DEEP2', 'DEEP_XMM',
                                   'DEEP_ELAIS', 'SSP_AEGIS', 'WIDE_WIDE01',
                                   'WIDE_GAMA09', 'WIDE_GAMA15', 'WIDE_VVDS',
                                   'WIDE_HECTOMAP']

                self.field_map = {
                    'udeep_cosmos': (IDR_URL +
                                     '/hsc_ssp/dr1/s16a/doc/fig/' +
                                     'tracts_patches_UD_COSMOS_HSC-I.png'),
                    'udeep_sxds': (IDR_URL +
                                   '/hsc_ssp/dr1/s16a/doc/fig/' +
                                   'tracts_patches_UD_SXDS_HSC-I.png'),
                    'deep_cosmos': (IDR_URL +
                                    '/hsc_ssp/dr1/s16a/doc/fig/' +
                                    'tracts_patches_D_COSMOS_HSC-I.png'),
                    'deep_deep2': (IDR_URL +
                                   '/hsc_ssp/dr1/s16a/doc/fig/' +
                                   'tracts_patches_D_DEEP2-3_HSC-I.png'),
                    'deep_xmm': (IDR_URL +
                                 '/hsc_ssp/dr1/s16a/doc/fig/' +
                                 'tracts_patches_D_XMM-LSS_HSC-I.png'),
                    'deep_elais': (IDR_URL +
                                   '/hsc_ssp/dr1/s16a/doc/fig/' +
                                   'tracts_patches_D_XMM-ELAIS-N1-I.png'),
                    'ssp_aegis': (IDR_URL +
                                  '/hsc_ssp/dr1/s16a/doc/fig/' +
                                  'tracts_patches_W_AEGIS_HSC-I.png'),
                    'wide_gama09': (IDR_URL +
                                    '/hsc_ssp/dr1/s16a/doc/fig/' +
                                    'tracts_patches_W_GAMA09H_HSC-I.png'),
                    'wide_gama15': (IDR_URL +
                                    '/hsc_ssp/dr1/s16a/doc/fig/' +
                                    'tracts_patches_W_GAMA15H_HSC-I.png'),
                    'wide_hectomap': (IDR_URL +
                                      '/hsc_ssp/dr1/s16a/doc/fig/' +
                                      'tracts_patches_W_HECTOMAP_HSC-I.png'),
                    'wide_vvds': (IDR_URL +
                                  '/hsc_ssp/dr1/s16a/doc/fig/' +
                                  'tracts_patches_W_VVDS_HSC-I.png'),
                    'wide_wide12': (IDR_URL +
                                    '/hsc_ssp/dr1/s16a/doc/fig/' +
                                    'tracts_patches_W_WIDE12H_HSC-I.png'),
                    'wide_xmm': (IDR_URL +
                                 '/hsc_ssp/dr1/s16a/doc/fig/' +
                                 'tracts_patches_W_XMM_HSC-I.png'),
                    'wide_wide01': (IDR_URL +
                                    '/hsc_ssp/dr1/s16a/doc/fig/' +
                                    'tracts_patches_W_WIDE01H_HSC-R.png')
                    }

            elif dr is 'dr2':
                # Useful URLs
                self.das_url = PDR_URL + "/das_quarry/dr2/cgi-bin/"
                self.psf_url = PDR_URL + "/psf/5/cgi/"

                # DR2 is being prepared right now! Issue a warning
                warnings.warn("DR2 is being prepared right now !")

                # Available filters
                """
                Notice:
                -------
                    Only UDEEP+DEEP fields have narrow-band coverage.
                """
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I',
                                    'HSC-Z', 'HSC-Y',
                                    'NB0387', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y',
                                          'nb0387', 'nb816', 'nb921']

                # Available reruns
                """
                We ignore the early data release and the s15a release.
                All their data are covered by s15b and s16a release, and
                the data quality has been significantly improved
                """
                self.rerun_list = ['any',
                                   's17a_dud', 's17a_wide']

                # Available fields
                """
                Notice:
                -------
                    From S17A, the DEEP and UDEEP fields are released together.
                """
                self.field_list = ['DUD_COSMOS', 'UDEEP_SXDF',
                                   'DEEP_COSMOS', 'DEEP_DEEP2', 'DEEP_XMM',
                                   'DEEP_ELAIS', 'SSP_AEGIS', 'WIDE_WIDE01',
                                   'WIDE_GAMA09', 'WIDE_GAMA15', 'WIDE_VVDS',
                                   'WIDE_HECTOMAP']

                self.field_map = {
                    'dud_cosmos': (IDR_URL +
                                   '/hsc_ssp/dr2/s17a/doc/fig/' +
                                   'tracts_patches_DUD_COSMOS_HSC-I.png'),
                    'dud_deep2': (IDR_URL +
                                  '/hsc_ssp/dr2/s17a/doc/fig/' +
                                  'tracts_patches_DUD_DEEP2-3_HSC-I.png'),
                    }
            else:
                raise DrException("!! Wrong information about data release !!")

    def _get_credential(self, pdr=False):
        """Get the username and password for HSC database.
        """
        if not pdr:
            try:
                self._username = os.environ['HSC_IDR_USR']
                self._password = os.environ['HSC_IDR_PWD']
            except KeyError:
                import getpass
                get_input = input
                if sys.version_info[:2] <= (2, 7):
                    get_input = raw_input
                self._username = get_input("Internal Data Release Username : ")
                self._password = getpass.getpass("Password : ")
        else:
            try:
                self._username = os.environ['HSC_PDR_USR']
                self._password = os.environ['HSC_PDR_PWD']
            except KeyError:
                import getpass
                get_input = input
                if sys.version_info[:2] <= (2, 7):
                    get_input = raw_input
                self._username = get_input("Public Data Release Username : ")
                self._password = getpass.getpass("Password : ")
