#!/usr/bin/env python
# -*- coding: utf-8 -*-

import string
import random

import astropy.units as u

__all__ = ['same_string', 'random_string', 'r_phy_to_ang']


def _passively_decode_string(a):
    try:
        return a.decode()
    except AttributeError:
        return a


def same_string(a, b):
    """Compare string in a Python2 and Python3-safe way.

    Shamelessly "borrowed" from Halotools by Andrew Hearin

    Parameters:
    -----------
    a: str
        String A to be compared.
    b: str
        String B to be compared
    """
    a = _passively_decode_string(a)
    b = _passively_decode_string(b)
    return a == b


def random_string(length=5, chars=string.ascii_uppercase + string.digits):
    """
    Random string generator.

    Based on:
    http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python

    Parameters:
    -----------
    length: int
        Length of the string. Default: 5
    chars: string object
        Types of characters allowed in the random string. Default: ASCII_Uppercase + Digits.
    """
    return ''.join(random.choice(chars) for _ in range(length))


def r_phy_to_ang(r_phy, redshift, cosmo=None, phy_unit='kpc', ang_unit='arcsec'):
    """
    Convert physical radius into angular size.
    """
    if cosmo is None:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    return (r_phy * u.Unit(phy_unit) /
            cosmo.kpc_proper_per_arcmin(redshift)).to(u.Unit(ang_unit))
