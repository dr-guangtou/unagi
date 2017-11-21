#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


class HSCObject(object):
    """Class for HSC object

    Args:
        ra (float)  : RA of the object.
        dec (float) : Dec of the object.
    """
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec

    def download_coadd(self, size, filter):
        """

        """
