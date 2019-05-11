#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Core functions"""

__all__ = ['HscServer']

from . import config

class Hsc():
    """
    HSC Server Class.
    """
    PIXEL_SIZE = 0.168   # arcsec / pixel
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
        assert dr in DATABASE
        self.rerun = config.rerrun(dr=dr, rerun=rerun, pdr=pdr, config_fiel=config_file)
    
    def _login(self, username=None, password=None):
        """
        Login to HSC server.

        Parameters
        ----------
        username : str
        password : str
        """
        pass

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