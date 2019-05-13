#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Core functions"""

__all__ = ['Hsc']

import urllib

import astropy.units as u
import astropy.coordinates as coord

from . import config

class Hsc():
    """
    HSC Server Class.
    """
    PIXEL_SIZE = 0.168   # arcsec / pixel
    MAX_CUTOUT = 2116 * u.arcsec
    DATABASE = ['pdr1', 'pdr2', 'dr1', 'dr2']
    FILTER_SHORT = ['g', 'r', 'i', 'z', 'y', 'nb0387', 'nb816', 'nb921']
    FILTER_LIST = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y',
                   'NB0387', 'NB0816', 'NB0921']

    def __init__(self, dr='dr2', rerun='s18a_wide', pdr=False, config_file=None):
        """
        Initialize a HSC rerun object.

        Parameters
        ----------
        dr : str
            HSC database name. Default: 'dr2'
        rerun : str
            HSC rerun dataset name. Default: 's18a_wide'
        pdr: bool
            Using public data release. Default: False
        config_file: str
            Name of the configuration file. Default: None
        """
        # Initiate the Rerun object
        assert dr in self.DATABASE
        self.rerun = config.Rerun(dr=dr, rerun=rerun, pdr=pdr, config_file=config_file)

        # Whether login to the server
        self.is_login = False
        self.opener = None
        # Try to login to the HSC archive
        self.login()

    def login(self, username=None, password=None):
        """
        Login to HSC server.

        Parameters
        ----------
        username : str
        password : str
        """
        # Get the user name and password
        if username is None:
            username = self.rerun._username
        if password is None:
            password = self.rerun._password

        # Create a password manager
        password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()

        # Add the username and password.
        password_mgr.add_password(None, self.rerun.base_url, username, password)
        handler = urllib.request.HTTPBasicAuthHandler(password_mgr)

        # create "opener" (OpenerDirector instance)
        opener = urllib.request.build_opener(handler)

        try:
            # use the opener to fetch a URL
            opener.open(self.rerun.base_url)
            # Install the opener.
            urllib.request.install_opener(opener)
            self.is_login = True
            self.opener = opener
        except urllib.error.HTTPError as e:
            print("! Can not login to HSC archive: %s" % str(e))
            self.opener = None

    def _form_cutout_url(self, coordinate, ):
        """
        Form the URL to download HSC cutout images.

        Please see details here:
            https://hscdata.mtk.nao.ac.jp/das_quarry/dr2.1/manual.html

        Example:
            1. Using the (RA, Dec) of the two diagonal corners.
            https://hscdata.mtk.nao.ac.jp/das_quarry/dr2.1/cgi-bin/cutout?ra1=135.0&dec1=0.0&ra2=135.1&dec2=0.1&type=coadd&image=on&filter=HSC-I&tract=9560&rerun=s18a_wide

            2. Uisng the (RA, Dec) of the center and the width/height of the image.
            https://hscdata.mtk.nao.ac.jp/das_quarry/dr2.1/cgi-bin/cutout?ra=135.0&dec=0.0&sw=20arcsec&sh=20arcsec&type=coadd&image=on&filter=HSC-I&tract=&rerun=s18a_wide

        About the image type:
            1. `coadd`: will return a FITS image.
            2. `wrap`: will return a `.tar` compressed file that containts all the wrapped images.

        Limitation:
            Maximum cutout image size is 3 x 3 `Patches`. Each `Patch` is 4200 x 4200 pixels.
            So the maximum size is 35 arcmin.

        """
        pass

    def _form_image_url(self, coordinate, ):
        """
        Form the URL to directly download HSC files.
        """
        pass

    def _form_psf_url(self, coordinate, ):
        """
        Form the URL to download HSC PSF model.
        """
        pass

    def _check_filter(self, filt):
        """
        Check if the filter is available in the rerun.

        If the filter string is in the short format, convert it into the formal one.

        Parameters:
        -----------
        filt: str
            Name of the filter.
        """
        if filt not in self.rerun.filter_list:
            if filt in self.rerun.filter_list_short:
                filt = self.rerun.filter_list[self.rerun.filter_list_short.index(filt)]
            else:
                raise ValueError('Unknown filter: {}'.format(filt))

        return filt

    def _download_file(self, url, file_path=None):
        """
        Function to download file from server.
        """
        pass

    def _parse_cutout_size_center(self, width, height):
        """
        Check the size of the cutout image using width and height.

        Parameters:
        -----------
        width : float
            Image width in unit of arcsec.
        height : float
            Image height in unit of arcsec.
        """
        pass
