
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sky background related"""

import os
import copy
import warnings

import numpy as np
from scipy.stats import sigmaclip

from astropy.table import Table

__all__ = ['SkyObjs']


class SkyObjs():
    """
    Class for HSC sky objects.
    """
    # Convert the flux from erg/s/cm^2/Hz to HSC image value
    CGS_TO_MUJY = 1.7378E30

    # List of filters
    FILTER_LIST = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y',
                   'NB0387', 'NB0816', 'NB0921']
    # Nicknames of filters
    FILTER_SHORT = ['g', 'r', 'i', 'z', 'y', 'nb0387', 'nb816', 'nb921']

    def __init__(self, skyobjs, meas=False, nobj_min=5):
        """
        Initialize an object for HSC sky object catalog.
        """
        # Whether it is a forced photometry or a measurement catalog
        if meas:
            self.ra_col = 'i_ra'
            self.dec_col = 'i_dec'
            self.type = 'meas'
            self.meas = True
        else:
            self.ra_col = 'ra'
            self.dec_col = 'dec'
            self.type = 'force'
            self.meas = False

        # If skyobjs is a file name, read in the catalog
        if isinstance(skyobjs, str):
            _, file_ext = os.path.splitext(skyobjs)
            if file_ext == '.npy':
                self.skyobjs = np.load(skyobjs)
            elif file_ext == '.fits':
                self.skyobjs = Table.read(skyobjs).as_array().data
            else:
                raise TypeError("# Wrong file type: npy or fits!")
        elif isinstance(skyobjs, Table):
            try:
                self.skyobjs = skyobjs.as_array().data
            except Exception:
                self.skyobjs = skyobjs.as_array()
        elif isinstance(skyobjs, np.ndarray) or isinstance(skyobjs, np.recarray):
            self.skyobjs = skyobjs

        # Minimum number of sky objects
        self.n_min = nobj_min

        # List of Tracts
        self.tract_list = list(np.unique(self.skyobjs['tract']))
        self.n_tract = len(self.tract_list)

    def select_tract(self, tract, patch=None, n_min=10):
        """Select sky objects on one Tract (and Patch) from the catalog """
        tract_mask = self.skyobjs['tract'] == tract
        if tract_mask.sum() == 0:
            warnings.warn("# Tract {0} is not available!".format(tract))
            return None

        if patch is not None:
            tract_mask = tract_mask & (self.skyobjs['patch'] == patch)
            if tract_mask.sum() == 0:
                warnings.warn("# Tract {0}-Patch {1} is not available!".format(tract, patch))
                return None

        # Number of sky objects available
        n_skyobj = tract_mask.sum()
        if n_skyobj <= n_min:
            warnings.warn("# Less than {0} skyobjs is available.")
            return None

        return copy.deepcopy(self.skyobjs[tract_mask])

    def flux_stats(self, flux_col, sigma=5, to_mujy=True):
        """Basic statistics of the flux."""
        u_factor = self.CGS_TO_MUJY if to_mujy else 1.0

        try:
            flux = self.skyobjs[flux_col] * u_factor
        except ValueError:
            raise Exception("# Wrong flux column name: {0}".format(flux_col))

        # Only use the ones with a good flux
        mask = np.isfinite(flux)
        if mask.sum() <= self.n_min:
            warnings.warn("# Does not have enough sky object: {0}".format(mask.sum()))
        flux_use = flux[mask]

        # Sigma clipping
        if sigma is not None and sigma > 0:
            flux_use, flux_low, flux_upp = sigmaclip(flux_use, low=sigma, high=sigma)
        if len(sigma) <= self.n_min:
            warnings.warn("# Does not have enough sky object: {0}".format(len(sigma)))

        return np.mean(flux_use), np.median(flux_use), np.std(flux_use)
