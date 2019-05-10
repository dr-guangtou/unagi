#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from astropy.coordinates import SkyCoord

__all__ = ('SspObject',)


class SspObject(object):
    """Class for HSC SSP object

    Examples
    --------
    The examples below illustrate common usage of the `SspObject` object.

        >>> from unagi.object import SspObject
        >>> object_1 = SspObject(30.0, -2.0)

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

    @property
    def galactic(self):
        return self.icrs.galactic

    @property
    def galactic_l(self):
        return self.galactic.l.value

    @property
    def galactic_b(self):
        return self.galactic.b.value

    def distance_to(self, *args, **kwargs):
        """Get the angular distance to another object.

        Examples
        --------

            >>> object_1 = SspObject(30.0, -2.0)
            >>> object_2 = SspObject(31.0, -1.9)
            >>> ang_sep = object_1.distance_to(object_2)

        Or:
            >>> coord_2 = SkyCoord(31.0, -1.9, unit='deg', frame='icrs')
            >>> ang_sep = object_1.distance_to(coord_2)

        Or:
            >>> ang_sep = object_1.distance_to(ra=30.0, dec=-1.9)
            >>> ang_sep = object_1.distance_to(l=100.0, b=20.0)

        Parameters
        ----------
            another_object : object (either SspObject or SkyCoord)
                Another HSC object or a SkyCoord object.
            ra, dec : float (optional)
                ICRS coordinate of an object.
            l, b : float (optional)
                Galactic coordinate of another object.
        """
        if args:
            if isinstance(args[0], SspObject):
                return self.icrs.separation(args[0].icrs)
            elif isinstance(args[0], SkyCoord):
                return self.icrs.separation(args[0])
            else:
                raise TypeError("Wrong type of object ",
                                "Should be SspObject or SkyCoord")
        elif ('ra' in kwargs) and ('dec' in kwargs):
            return self.icrs.separation(SkyCoord(kwargs['ra'],
                                                 kwargs['dec'],
                                                 unit="deg",
                                                 frame=self.frame))
        elif ('l' in kwargs) and ('b' in kwargs):
            return self.icrs.separation(SkyCoord(kwargs['l'], kwargs['b'],
                                                 unit="deg",
                                                 frame="galactic").icrs)
        else:
            raise Exception("Something wrong !")

    def in_region(self, *args):
        """Whether the object is within specific region.
        Examples
        --------

            >>> object_1 = SspObject(30.0, -2.0)
            >>> object_1.in_region(Tract)

        Parameters
        ----------
        """
        raise NotImplementedError

    def _sky_coord(self):
        """
        Create a SkyCoord object.
        """
        coord = SkyCoord(self.ra, self.dec, unit="deg",
                         frame=self.frame)
        if self.frame is not 'icrs':
            return coord.icrs
        else:
            return coord
