#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from astropy.coordinates import SkyCoord
import astropy.units as u

__all__ = ['HscObject']


class HscObject(object):
    """Class for HSC object

    Examples
    --------
    The examples below illustrate common usage of the `HscObject` object.

        >>> from unagi.object import HscObject
        >>> object_1 = HscObject(30.0, -2.0)

    Parameters
    ----------
    ra : float
        RA of the object.
    dec : float
        Dec of the object.
    frame : string, optional
        Type of coordinate frame.  Default is 'icrs'
    """
    def __init__(self, ra, dec, frame='icrs'):
        self.ra = ra
        self.dec = dec
        self.frame = frame

    @property
    def icrs(self):
        return self._sky_coord()

    def _sky_coord(self):
        """
        Create a SkyCoord object.
        """
        coord = SkyCoord(self.ra, self.dec, unit="deg", frame=self.frame)
        if self.frame is not 'icrs':
            return coord.icrs
        else:
            return coord

    def download_coadd(self, size, filter):
        """

        """
