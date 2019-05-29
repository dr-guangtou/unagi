#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Core functions"""

import os
import ssl
import json
import urllib
import shutil
import warnings

import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.utils.data import download_file

from . import config

__all__ = ['Hsc', 'DEFAULT_CUTOUT_CENTER', 'DEFAULT_CUTOUT_CORNER',
           'IMG_HDU', 'MSK_HDU', 'VAR_HDU']

DEFAULT_CUTOUT_CENTER = {
    'ra': '', 'dec': '', 'sw': '5arcsec', 'sh': '5arcsec',
    'type': 'coadd', 'image': 'on', 'mask': 'off', 'variance': 'off',
    'filter': 'HSC-I', 'rerun': ''
}

DEFAULT_CUTOUT_CORNER = {
    'ra1': '', 'dec1': '', 'ra2': '', 'dec2': '',
    'type': 'coadd', 'image': 'on', 'mask': 'off', 'variance': 'off',
    'filter': 'HSC-I', 'rerun': ''
}

SQL_OUTPUT_FORMAT = ['csv', 'csv.gz', 'sqlite3', 'fits']

IMG_HDU = 1
MSK_HDU = 2
VAR_HDU = 3

class HscException(Exception):
    """Class for error related to data release information.
    """
    pass

class QueryError(Exception):
    """Class for error related to SQL query.
    """
    pass

class Hsc():
    """
    HSC Server Class.
    """
    # HSC Pixel size in unit of arcsec / pixel
    PIXEL_SIZE = 0.168
    # The maximum allowed cutout size in unit of arcsec (~35 arcmin)
    MAX_CUTOUT = 2116 * u.arcsec
    # Default cutout image size in unit of arcsec
    DEFAULT_IMG_SIZE = 10.0 * u.arcsec
    # Available HSC database
    # TODO: dr2-citus is not supported yet
    DATABASE = ['pdr1', 'pdr2', 'dr1', 'dr2']
    # List of HSC filters
    FILTER_LIST = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y',
                   'NB0387', 'NB0816', 'NB0921']
    FILTER_SHORT = ['g', 'r', 'i', 'z', 'y', 'nb0387', 'nb816', 'nb921']

    def __init__(self, dr='dr2', rerun='s18a_wide', config_file=None):
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
        self.dr = dr

        # Check whether we are using public or internal data release.
        if dr[0] == 'p':
            self.pdr = True
        else:
            self.pdr = False

        self.rerun = rerun
        self.archive = config.Rerun(
            dr=self.dr, rerun=self.rerun, config_file=config_file)

        # SQL client version
        # TODO: figure out how to get this from HSC archive
        self.sql_version = 20181012.1

        # Whether login to the server
        self.is_login = False
        self.opener = None
        # Try to login to the HSC archive
        if not self.is_login or self.opener is None:
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
            username = self.archive._username
        if password is None:
            password = self.archive._password

        # Create a password manager
        password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()

        # Add the username and password.
        password_mgr.add_password(None, self.archive.base_url, username, password)
        handler = urllib.request.HTTPBasicAuthHandler(password_mgr)

        # Create "opener" (OpenerDirector instance)
        self.opener = urllib.request.build_opener(handler)

        try:
            # Install the opener.
            self.opener.open(self.archive.base_url)
            urllib.request.install_opener(self.opener)
            self.is_login = True
        except urllib.error.HTTPError as e:
            print("! Can not login to HSC archive: %s" % str(e))
            self.opener = None
    
    def logout(self):
        """
        Log out of the HSC server.

        TODO: This actually does not do anything.
        """
        if not self.is_login:
            print("# Has not login to HSC archive yet!")
        else:
            print("# Log out of HSC archive now!")
            self.opener = None
            self.is_login = False

    def download_cutout(self, coord, output_file, coord_2=None, w_half=None, h_half=None,
                        filt='HSC-I', img_type='coadd', image=True, variance=False, mask=False,
                        overwrite=True):
        """
        Download coadded or warped cutout image(s).

        Parameters:
        -----------
        """
        cutout_kwargs = {'coord_2': coord_2, 'w_half': w_half, 'h_half': h_half,
                         'filt': filt, 'img_type': img_type, 'image': image,
                         'variance': variance, 'mask': mask}

        if img_type == 'coadd':
            # Download FITS file for coadd image.
            cutout = self.get_cutout_image(coord, **cutout_kwargs)
            _ = cutout.writeto(output_file, overwrite=overwrite)
            return cutout
        elif img_type == 'warp':
            # Download the tarball for warpped images.
            cutout_url = self.get_cutout_image(coord, **cutout_kwargs)
            if os.path.isfile(output_file) and not overwrite:
                raise HscException("# File {} exists!".format(output_file))
            else:
                _ = shutil.move(download_file(cutout_url, show_progress=False), output_file)
            return cutout_url
        else:
            raise HscException("# Wrong image type: coadd or warp !")

    def get_cutout_image(self, coord, coord_2=None, w_half=None, h_half=None, filt='HSC-I',
                         img_type='coadd', image=True, variance=False, mask=False, verbose=False):
        """
        Get HSC cutout image.

        Parameters:
        -----------
        """
        cutout_kwargs = {'filt': filt, 'img_type': img_type, 'image': image,
                         'variance': variance, 'mask': mask}

        cutout_url = self.form_cutout_url(
            coord, coord_2=coord_2, w_half=w_half, h_half=h_half, **cutout_kwargs)

        if img_type == 'warp':
            if verbose:
                warnings.warn("# Not a coadd cutout, will return the url")
            return cutout_url

        try:
            if verbose:
                print("# Downloading FITS image from {}".format(cutout_url))
            cutout = fits.open(cutout_url)
        except urllib.error.HTTPError as e:
            print("# Error message: {}".format(e))
            raise Exception("# Can not download cutout: {}".format(cutout_url))

        return cutout

    def get_psf_model(self, coord, filt='HSC-I', img_type='coadd', centered=True, verbose=False):
        """
        Get the PSF model at a given sky position.

        Parameters:
        -----------
        """
        # Check the filter
        filt = self._check_filter(filt)

        # Centered the PSF model or not
        center_psf = 'on' if centered else 'off'

        # Image type
        if img_type is not 'coadd' and img_type is not 'warp':
            raise HscException("# Wrong image type !")

        # Default dict for generating PSF
        psf_dict = {'rerun': self.rerun, 'filter': filt, 'img_type': img_type,
                    'centered': center_psf}

        ra_str, dec_str = self._parse_coordinate(coord)
        psf_dict['ra'] = ra_str
        psf_dict['dec'] = dec_str

        psf_url = self.archive.psf_url + '&'.join(
            key + '=' + value for key, value in psf_dict.items())

        if img_type == 'warp':
            if verbose:
                warnings.warn("# Not a coadd PSF model, will return the url")
            return psf_url

        try:
            if verbose:
                print("# Downloading FITS image from {}".format(psf_url))
            psf_model = fits.open(psf_url)
        except urllib.error.HTTPError as e:
            print("# Error message: {}".format(e))
            raise Exception("# Can not download cutout: {}".format(psf_url))

        return psf_model

    def form_cutout_url(self, coord, coord_2=None, w_half=None, h_half=None, **kwargs):
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
            2. `warp`: will return a `.tar` compressed file that containts all the warpped images.

        Limitation:
            Maximum cutout image size is 3 x 3 `Patches`. Each `Patch` is 4200 x 4200 pixels.
            So the maximum size is 35 arcmin.

        """
        # Default image cutout size
        if w_half is None:
            w_half = self.DEFAULT_IMG_SIZE

        if h_half is None:
            h_half = self.DEFAULT_IMG_SIZE

        if coord_2 is not None:
            cutout_dict = self._parse_cutout_corner(coord, coord_2, **kwargs)
        else:
            cutout_dict = self._parse_cutout_center(coord, w_half, h_half, **kwargs)

        return self.archive.img_url + '&'.join(
            key + '=' + value for key, value in cutout_dict.items())

    def _form_image_url(self, coord, ):
        """
        Form the URL to directly download HSC files.
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
        if filt not in self.FILTER_LIST:
            if filt in self.FILTER_SHORT:
                filt = self.FILTER_LIST[self.FILTER_SHORT.index(filt)]
            else:
                raise ValueError('Unknown filter: {}'.format(filt))

        return filt

    def _download_file(self, url, file_path=None):
        """
        Function to download file from server.
        """
        pass

    def _parse_coordinate(self, coord, frame='icrs'):
        """
        Convert the coordinate into string for RA & Dec.

        Parameters:
        -----------
        coord: astropy.coordinates.SkyCoord object
            Sky coordinate.
        frame: str
            Name of the coordinate frame. Default: 'icrs'
        """
        return coord.transform_to(frame).to_string('decimal').split(' ')

    def _parse_size_center(self, w_half, h_half, correct=False):
        """
        Convert the image width and height into string.

        Parameters:
        -----------
        w_half: float
            Half of the image width.
        h_half: float
            Half of the image height.
        correct: bool
            Correct the bool
        """
        # Check the image width
        if w_half * 2 > self.MAX_CUTOUT:
            if correct:
                print("# Image width is too large! Just use the maximum size!")
                w_half = self.MAX_CUTOUT / 2.0
            else:
                raise Exception("# Image width is too large")

        # Check the image height
        if h_half * 2 > self.MAX_CUTOUT:
            if correct:
                print("# Image width is too large! Just use the maximum size!")
                h_half = self.MAX_CUTOUT / 2.0
            else:
                raise Exception("# Image width is too large")

        # Convert the sizes into string
        w_str = "{:f}{}".format(w_half.value, w_half.unit)
        h_str = "{:f}{}".format(h_half.value, h_half.unit)

        return w_str, h_str

    def _parse_cutout_center(self, coord, w_half, h_half, filt='HSC-I',
                             img_type='coadd', image=True,
                             variance=False, mask=False):
        """
        Organize the parameters to generate the cutout image.

        Parameters:
        -----------
        coord: astropy.coordinates.SkyCoord object
            Sky coordinate.
        w_half: float
            Half of the image width.
        h_half: float
            Half of the image height.
        """
        # Load the default dictionary
        cutout_dict = DEFAULT_CUTOUT_CENTER

        # Parse the coordinate
        ra_str, dec_str = self._parse_coordinate(coord)
        cutout_dict['ra'] = ra_str
        cutout_dict['dec'] = dec_str

        # Parse the image width and height
        w_str, h_str = self._parse_size_center(w_half, h_half)
        cutout_dict['sw'] = w_str
        cutout_dict['sh'] = h_str

        # Image filter
        filt = self._check_filter(filt)
        cutout_dict['filter'] = filt

        # Image type
        cutout_dict['type'] = img_type

        # Content of the cutout
        cutout_dict['image'] = 'on' if image else 'off'
        cutout_dict['variance'] = 'on' if variance else 'off'
        cutout_dict['mask'] = 'on' if mask else 'off'

        # Rerun name
        cutout_dict['rerun'] = self.rerun

        return cutout_dict

    def _parse_cutout_corner(self, coord1, coord2, filt='HSC-I',
                             img_type='coadd', image=True,
                             variance=False, mask=False):
        """
        Organize the parameters to generate the cutout image using corner coordinates.

        Parameters:
        -----------
        coord: astropy.coordinates.SkyCoord object
            Sky coordinate.
        w_half: float
            Half of the image width.
        h_half: float
            Half of the image height.
        """
        # Load the default dictionary
        cutout_dict = DEFAULT_CUTOUT_CORNER

        # Check the size of the cutout image
        if ((np.abs(coord2.ra - coord1.ra) >= self.MAX_CUTOUT) or (
                np.abs(coord2.dec - coord1.dec) >= self.MAX_CUTOUT)):
            raise Exception("# Image width is too large")

        # Parse the coordinate
        ra1_str, dec1_str = self._parse_coordinate(coord1)
        cutout_dict['ra1'] = ra1_str
        cutout_dict['dec1'] = dec1_str

        ra2_str, dec2_str = self._parse_coordinate(coord2)
        cutout_dict['ra2'] = ra2_str
        cutout_dict['dec2'] = dec2_str

        # Image filter
        filt = self._check_filter(filt)
        cutout_dict['filter'] = filt

        # Image type
        cutout_dict['type'] = img_type

        # Content of the cutout
        cutout_dict['image'] = 'on' if image else 'off'
        cutout_dict['variance'] = 'on' if variance else 'off'
        cutout_dict['mask'] = 'on' if mask else 'off'

        # Rerun name
        cutout_dict['rerun'] = self.rerun

        return cutout_dict

    def _http_post(self, url, data, headers):
        """
        Request data.

        Based on the NAOJ script:
        https://hsc-gitlab.mtk.nao.ac.jp/snippets/13
        """
        req = urllib.request.Request(url, data.encode('utf-8'), headers)
        res = urllib.request.urlopen(req)
        return res

    def _http_post_json(self, url, data):
        """
        Send SQL request.

        Based on the NAOJ script:
        https://hsc-gitlab.mtk.nao.ac.jp/snippets/13
        """
        data['clientVersion'] = self.sql_version
        post_data = json.dumps(data)
        return self._http_post(url, post_data, {'Content-type': 'application/json'})

    def _credential(self):
        """
        Return a dict that contains username and password.
        """
        return {'account_name': self.archive._username,
                'password': self.archive._password}

    def submit_query(self, sql, out_format, nomail=True, skip_syntax=True):
        """
        Submit SQL job to HSC archive.
        """
        url = os.path.join(self.archive.cat_url, 'submit')

        if out_format.strip().lower() not in SQL_OUTPUT_FORMAT:
            raise NameError("Wrong output format: ['csv', 'csv.gz', 'sqlite3', 'fits']")

        catalog_job = {
            'sql'                     : sql,
            'out_format'              : out_format.strip().lower(),
            'include_metainfo_to_body': True,
            'release_version'         : self.dr,
            }

        post_data = {
            'credential': self._credential(), 
            'catalog_job': catalog_job, 
            'nomail': nomail, 
            'skip_syntax_check': skip_syntax
            }

        res = self._http_post_json(url, post_data)
        job = json.load(res)
        return job
