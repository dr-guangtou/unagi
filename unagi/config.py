#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
import sys

# Update: 2017-11-22
HSC_FILTERS = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y',
               'NB0387', 'NB0816', 'NB0921']

HSC_FSHORT = ['g', 'r', 'i', 'z', 'y', 'nb0387', 'nb0816', 'nb0921']

HSC_IDR_RERUN = ['s15b_udeep', 's15b_deep', 's15b_wide',
                 's16a_udeep', 's16a_deep', 's16a_wide', 's16a_wide2',
                 's17a_dud', 's17a_wide', 'any']

HSC_DR1_URL = "https://hscdata.mtk.nao.ac.jp/das_quarry/dr1/cgi-bin/"
HSC_DR2_URL = "https://hscdata.mtk.nao.ac.jp/das_quarry/dr2/cgi-bin/"

HSC_DR1_PSF = "https://hscdata.mtk.nao.ac.jp/psf/4/cgi"
HSC_DR2_PSF = "https://hscdata.mtk.nao.ac.jp/psf/5/cgi"

PDR_URL = "https://hsc-release.mtk.nao.ac.jp"


class DrException(Exception):
    """Class for error related to data release information.
    """
    pass


class HscDasConfig(object):
    """Class for configuration parameters for HSC database.

    Examples
    --------
    The examples below illustrate common usage of the `HscObject` object.

        >>> from unagi.config import HscDasConfig
        >>> object_1 = HscObject(30.0, -2.0)

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

            if dr == 'dr1':
                # PDR = Public data release
                self.database = 'PDR'

                # Useful URLs
                self.das_url = PDR_URL + "/das_quarry/cgi-bin/"
                self.psf_url = PDR_URL + "/psf/pdr1/cgi/"

                # Available filters
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I',
                                    'HSC-Z', 'HSC-Y', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y',
                                          'nb816', 'nb921']

                # Available reruns
                self.rerun_list = ['any', 'pdr1_udeep', 'pdr1_deep',
                                   'pdr1_wide']

                # Available fields
                self.field_list = ['UDEEP_COSMOS', 'UDEEP_SXDF',
                                   'DEEP_COSMOS', 'DEEP_DEEP2', 'DEEP_XMM',
                                   'DEEP_ELAIS', 'WIDE_AEGIS',
                                   'WIDE_GAMA09', 'WIDE_GAMA15', 'WIDE_VVDS',
                                   'WIDE_HECTOMAP']

                self.field_map = {
                    'udeep_cosmos': (PDR_URL +
                                     '/das_search/pdr1/images/pdr1/' +
                                     'tracts_patches_UD_COSMOS_HSC-I.png'),
                    'udeep_sxds': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                   'tracts_patches_UD_SXDS_HSC-I.png'),
                    'deep_cosmos': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                    'tracts_patches_D_COSMOS_HSC-I.png'),
                    'deep_deep2': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                   'tracts_patches_D_DEEP2-3_HSC-I.png'),
                    'deep_xmm': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                 'tracts_patches_D_XMM-LSS_HSC-I.png'),
                    'deep_elais': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                   'tracts_patches_D_XMM-ELAIS-N1-I.png'),
                    'wide_aegis': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                   'AEGIS_i_area.png'),
                    'wide_gama09': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                    'W-GAMA09H_i_area.png'),
                    'wide_gama15': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                    'W-GAMA15H_i_area.png'),
                    'wide_hectomap': (PDR_URL +
                                      '/das_search/pdr1/images/pdr1/' +
                                      'W-HECTOMAP_i_area.png'),
                    'wide_vvds': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                  'W-VVDS_i_area.png'),
                    'wide_wide12': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                    'W-WIDE12H_i_area.png'),
                    'wide_xmm': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                 'W-XMMLSS_i_area.png')
                    }
            else:
                raise DrException("!! Wrong information about data release !!")
        else:
            """Use the HSC internal data release at:

            More information about the IDR/DAS query:

                https://hscdata.mtk.nao.ac.jp/das_quarry/dr1/manual.html

            So far, DR1 and DR2 are available
            """
            # Gather login information:
            self._get_credential(pdr=False)

            # Data release list
            self.dr_list = ['dr1']

            if dr is 'dr1':
                # PDR = Public data release
                self.database = 'PDR'

                # Useful URLs
                self.das_url = PDR_URL + "/das_quarry/cgi-bin/"
                self.psf_url = PDR_URL + "/psf/pdr1/cgi/"

                # Available filters
                self.filter_list = ['HSC-G', 'HSC-R', 'HSC-I',
                                    'HSC-Z', 'HSC-Y', 'NB0816', 'NB0921']
                self.filter_list_short = ['g', 'r', 'i', 'z', 'y',
                                          'nb816', 'nb921']

                # Available reruns
                self.rerun_list = ['any', 'pdr1_udeep', 'pdr1_deep',
                                   'pdr1_wide']

                # Available fields
                self.field_list = ['UDEEP_COSMOS', 'UDEEP_SXDF',
                                   'DEEP_COSMOS', 'DEEP_DEEP2', 'DEEP_XMM',
                                   'DEEP_ELAIS', 'WIDE_AEGIS',
                                   'WIDE_GAMA09', 'WIDE_GAMA15', 'WIDE_VVDS',
                                   'WIDE_HECTOMAP']

                self.field_map = {
                    'udeep_cosmos': (PDR_URL +
                                     '/das_search/pdr1/images/pdr1/' +
                                     'tracts_patches_UD_COSMOS_HSC-I.png'),
                    'udeep_sxds': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                   'tracts_patches_UD_SXDS_HSC-I.png'),
                    'deep_cosmos': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                    'tracts_patches_D_COSMOS_HSC-I.png'),
                    'deep_deep2': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                   'tracts_patches_D_DEEP2-3_HSC-I.png'),
                    'deep_xmm': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                 'tracts_patches_D_XMM-LSS_HSC-I.png'),
                    'deep_elais': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                   'tracts_patches_D_XMM-ELAIS-N1-I.png'),
                    'wide_aegis': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                   'AEGIS_i_area.png'),
                    'wide_gama09': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                    'W-GAMA09H_i_area.png'),
                    'wide_gama15': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                    'W-GAMA15H_i_area.png'),
                    'wide_hectomap': (PDR_URL +
                                      '/das_search/pdr1/images/pdr1/' +
                                      'W-HECTOMAP_i_area.png'),
                    'wide_vvds': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                  'W-VVDS_i_area.png'),
                    'wide_wide12': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                    'W-WIDE12H_i_area.png'),
                    'wide_xmm': (PDR_URL + '/das_search/pdr1/images/pdr1/' +
                                 'W-XMMLSS_i_area.png')
                    }

            elif dr is 'dr2':
                pass
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
