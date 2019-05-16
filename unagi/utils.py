#!/usr/bin/env python
# -*- coding: utf-8 -*-

import string
import random
import warnings

import numpy as np

import astropy.units as u

from scipy.stats import sigmaclip
from scipy.stats import gaussian_kde

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
    # Cosmology
    if cosmo is None:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # Convert the physical size into an Astropy quantity
    if not isinstance(r_phy, u.quantity.Quantity):
        r_phy = r_phy * u.Unit(phy_unit)

    return (r_phy / cosmo.kpc_proper_per_arcmin(redshift)).to(u.Unit(ang_unit))


def stats_summary(X, sigma=5.0, n_min=10, kde=True, bw=None):
    """
    Statistical summary of an array.
    """
    summary = {'low': np.nan, 'upp': np.nan, 'mean': np.nan,
               'median': np.nan, 'std': np.nan, 'kde': None, 'sigma': sigma}

    # Only use the ones with a good flux
    flag = np.isfinite(X)
    if flag.sum() <= n_min:
        warnings.warn("# Does not have enough elements: {0}".format(flag.sum()))
        return summary

    X = X[flag]

    # Sigma clipping
    if sigma is not None and sigma > 0:
        X_clipped, X_low, X_upp = sigmaclip(X, low=sigma, high=sigma)
        summary['low'] = X_low
        summary['upp'] = X_upp
    else:
        X_clipped = X
        summary['low'] = np.nanmin(X)
        summary['upp'] = np.nanmax(X)

    if len(X_clipped) <= n_min:
        warnings.warn("# Does not have enough sky object: {0}".format(len(sigma)))
        return summary

    # Mean, median, and standard deviation
    summary['mean'] = np.mean(X_clipped)
    summary['std'] = np.std(X_clipped)
    summary['median'] = np.median(X_clipped)

    if kde:
        if bw is None:
            bw = 0.2 * summary['std']
        summary['kde'] = gaussian_kde(X_clipped, bw_method=bw)

    return summary
