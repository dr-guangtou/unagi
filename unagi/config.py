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

            # PDR = Public data release
            self.database = 'PDR'

            if dr == 'dr1':
                self.data_release = 'dr1'

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

                self.field = [_PDR_UD_COSMOS, _PDR_UD_SXDS,
                              _PDR_D_COSMOS, _PDR_D_DEEP2, _PDR_D_ELAIS,
                              _PDR_D_XMM,
                              _PDR_W_XMM, _PDR_W_GAMA09, _PDR_W_GAMA15,
                              _PDR_W_WIDE12, _PDR_W_HECTOMAP,
                              _PDR_W_VVDS, _PDR_AEGIS]

            else:
                raise DrException("!! Wrong information about data release !!")

            self.field_table = Table(rows=self.fields)

            self.field_name = self.field_table['name'].data.astype('str')

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

            # PDR = Public data release
            self.database = 'IDR'

            if dr is 'dr1':
                self.data_release = 'dr1'

                # Useful URLs
                self.das_url = IDR_URL + "/das_quarry/dr1/cgi-bin/"
                self.psf_url = IDR_URL + "/psf/4/cgi/"
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

                self.field = [_IDR_UD_COSMOS, _IDR_UD_SXDS,
                              _IDR_D_COSMOS, _IDR_D_DEEP2, _IDR_D_ELAIS,
                              _IDR_D_XMM, _IDR_W_WIDE01,
                              _IDR_W_XMM, _IDR_W_GAMA09, _IDR_W_GAMA15,
                              _IDR_W_WIDE12, _IDR_W_HECTOMAP,
                              _IDR_W_VVDS, _IDR_AEGIS]

                self.field_table = Table(rows=self.fields)

                self.field_name = self.field_table['name'].data.astype('str')

            elif dr is 'dr2':
                # Useful URLs
                self.das_url = PDR_URL + "/das_quarry/dr2/cgi-bin/"
                self.psf_url = PDR_URL + "/psf/5/cgi/"
                self.map_url = IDR_URL + "/hsc_ssp/dr1/s17a/doc/fig/"
                self.txt_url = IDR_URL + "/hsc_ssp/dr1/s17a/doc/info/"

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
                self.rerun_default = 's17a_wide'
                self.wide_default = 's17a_wide'
                self.deep_default = 's17a_dud'
                self.udeep_default = 's17a_dud'

                # Available fields
                """
                Notice:
                -------
                    From S17A, the DEEP and UDEEP fields are released together.

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
                                             'NB0816'],
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

                self.field = [_IDR_DUD_COSMOS, _IDR_DUD_DEEP2,
                              _IDR_DUD_ELAIS, _IDR_DUD_XMM,
                              _IDR_W_WIDE01, _IDR_W_WIDE02, _IDR_W_WIDE03,
                              _IDR_W_WIDE04, _IDR_W_WIDE05, _IDR_W_WIDE06,
                              _IDR_W_WIDE07]

            else:
                raise DrException("!! Wrong information about data release !!")

            self.field_table = Table(rows=self.fields)

            self.field_name = self.field_table['name'].data.astype('str')

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
