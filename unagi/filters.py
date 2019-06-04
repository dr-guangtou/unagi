#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to deal with HSC filters."""

import os

import numpy as np

from astropy.table import Table

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams

import unagi

__all__ = ['filters_to_kcorrect', 'HscFilter']

plt.rc('text', usetex=True)
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'axes.titlepad': '15.0'})
rcParams.update({'font.size': 26})

# Directory that save the transmission curves of the filter
FILTER_DIR = os.path.join(os.path.dirname(unagi.__file__), 'data', 'filters')

# Formal names and nick names of the filters
FILTER_LIST = ['HSC-G', 'HSC-R', 'HSC-R2', 'HSC-I', 'HSC-I2', 'HSC-Z', 'HSC-Y',
               'HSC-NB387', 'HSC-NB816', 'HSC-NB921']
FILTER_SHORT = ['g', 'r', 'r2', 'i', 'i2', 'z', 'y', 'nb387', 'nb816', 'nb921']

# Colors to show the filter
FILTER_COLORS = ['#2ca02c', '#ff7f0e', '#ff7f0e', '#d62728', '#d62728',
                 '#8c564b', '#7f7f7f', 'c', 'm', 'purple']

class HscFilter(object):
    """Class for organizing HSC filters.

    Please see the following page for more details about HSC filters:
    https://www.subarutelescope.org/Observing/Instruments/HSC/sensitivity.html

    This class is based on the Filter() class from Ben Johnson's sedpy
    https://github.com/bd-j/sedpy/blob/master/sedpy/observate.py
    """
    def __init__(self, band, origin=False, center=False):
        """Read in the HSC filter transmission curves.
        """
        if band.strip().upper() in FILTER_LIST:
            self.index = FILTER_LIST.index(band.strip().upper())
        elif band.strip().lower() in FILTER_SHORT:
            self.index = FILTER_SHORT.index(band.strip().lower())
        else:
            raise NameError("# Wrong filter name!")

        # Name of the filter
        self.filter = FILTER_LIST[self.index]
        self.short = FILTER_SHORT[self.index]

        # Old or new version of r- and i-band filters
        if self.short == 'r' or self.short == 'i':
            self.version = 'old'
        else:
            self.version = 'new'

        # Whether the transmission curve is the total one or the just-filter one
        if origin:
            # This is the "origin" filter response curve from the HSC webpage:
            # https://www.subarutelescope.org/Observing/Instruments/HSC/sensitivity.html
            # Type of the transmission curve: area weighted one or for center
            filter_dir = os.path.join(FILTER_DIR, 'origin')
            self.type = 'just_filter'
            self.center = False
            if center:
                # This is for the center of the camera
                self.center = True
                self.name = '{}.txt'.format(self.filter)
            else:
                # This is the area-weighted mean one
                self.center = False
                self.name = 'w{}.txt'.format(self.filter)
        else:
            # This is the total transmission curves from Kawanomoto et al. 2018
            # Downloaded from here:
            # https://hsc-release.mtk.nao.ac.jp/doc/wp-content/uploads/2019/04/hsc_responses_all_rev3.tar.gz
            # Assume airmass=1.2 and PWV=1.5
            filter_dir = os.path.join(FILTER_DIR, 'total')
            self.type = 'total'
            self.name = 'hsc_{}_v2018.dat'.format(self.short)

        # Find the file and read in the transmission curve
        self.filename = os.path.join(filter_dir, self.name)
        self.wave, self.trans = self._load_filter()
        self.npts = len(self.wave)

        # Basic properties of the filter
        self.wave_effective = None
        self.wave_pivot = None
        self.wave_mean = None
        self.wave_average = None
        self.rectangular_width = None
        self.gauss_width = None
        self.effective_width = None
        self._basic_properties()

    def _load_filter(self):
        """Load and process the transimission curve."""
        if not os.path.isfile(self.filename):
            raise IOError("# Cannot find the response curve file {}".format(self.filename))

        # Read in the .txt response curve
        wave, trans = np.genfromtxt(self.filename, usecols=(0, 1), unpack=True)

        use = np.isfinite(trans) & (trans >= 0)
        order = wave[use].argsort()
        wave = wave[use][order]
        trans = trans[use][order]

        return wave, trans

    def print(self):
        """Print out basic properties of the filter."""
        print("# Filter: {}".format(self.filename))
        print("# Filter type: {}".format(self.type))
        print("# Version: {}".format(self.version))
        print("# Effective wavelength   : {:8.3f} Angstrom".format(self.wave_effective))
        print("# Pivot wavelength       : {:8.3f} Angstrom".format(self.wave_pivot))
        print("# Mean wavelength        : {:8.3f} Angstrom".format(self.wave_mean))
        print("# Effective width        : {:8.3f} Angstrom".format(self.effective_width))
        print("# Rectangular width      : {:8.3f} Angstrom".format(self.rectangular_width))

    def plot(self):
        """Plot the transmission curves of the filter."""
        fig = plt.figure(figsize=(8, 5))
        ax1 = fig.add_subplot(111)
        ax1.grid(linestyle='--', linewidth=2, alpha=0.8)
        ax1.axhline(0.0, linewidth=2, color='k', alpha=0.8)
        # Filled the transmission curve
        ax1.fill_between(self.wave, 0.0, self.trans, edgecolor='k', alpha=0.4,
                         linewidth=2.5, facecolor=FILTER_COLORS[self.index])
        ax1.text(0.9, 0.85, r'$\rm {}$'.format(self.short), transform=ax1.transAxes)
        # Highlight the effective and pivot wavelength
        ax1.axvline(self.wave_effective, linewidth=2.0, linestyle='--',
                    label=r'$\lambda_{\rm eff}$')
        ax1.axvline(self.wave_pivot, linewidth=2.0, linestyle='-.',
                    label=r'$\lambda_{\rm piv}$')
        _ = ax1.set_xlabel(r'$\mathrm{Wavelength}\ [\AA]$')
        _ = ax1.set_ylabel(r'$\mathrm{Transmission}$')
        ax1.legend(loc='lower right', fontsize=18)

        return fig

    def _basic_properties(self):
        """Get the basic properties of the filter.

        These properties include several 'effective' wavelength definitions and
        several width definitions.

        See Fukugita et al. (1996) AJ 111, 1748 for discussion and definition
        of many of these quantities.
        """
        # Calculate some useful integrals
        i0 = np.trapz(self.trans * np.log(self.wave), np.log(self.wave))
        i1 = np.trapz(self.trans, np.log(self.wave))
        i2 = np.trapz(self.trans * self.wave, self.wave)
        i3 = np.trapz(self.trans, self.wave)

        # Effective wavelength
        self.wave_effective = np.exp(i0 / i1)
        # Pivot wavelength
        self.wave_pivot = np.sqrt(i2 / i1)
        # Mean wavelength
        self.wave_mean = self.wave_effective
        self.wave_average = i2 / i3
        # Rectangular width of the filter
        self.rectangular_width = i3 / self.trans.max()

        i4 = np.trapz(self.trans * (np.log(self.wave / self.wave_effective)) ** 2.0,
                      np.log(self.wave))

        # Gaussian width of the filter
        self.gauss_width = (i4 / i1) ** (0.5)
        # Effectove width of the filter
        self.effective_width = (2.0 * np.sqrt(2. * np.log(2.)) *
                                self.gauss_width * self.wave_effective)



def filters_to_kcorrect(filename):
    """
    Convert a filter response curve to the Kcorrect format.

    This is used by Kcorrect and iSEDFit.
    """
    curve_file = os.path.join(FILTER_DIR, filename)
    if not os.path.isfile(curve_file):
        raise IOError("# Cannot find the response curve file {}".format(curve_file))

    # Read in the .txt response curve
    wave, response = np.genfromtxt(curve_file, usecols=(0,1), unpack=True)

    # Output file name
    prefix, band = os.path.splitext(filename)[0].split('-')
    output_par = os.path.join(
        FILTER_DIR, "{0}_{1}.par".format(prefix.lower().strip(), band.lower()))

    if os.path.isfile(output_par):
        print("# Curve {0} is already available".format(output_par))
    else:
        assert len(wave) == len(response), '''
            Wavelength and response curve should have the same size'''

        par = open(output_par, 'w')
        par.write(
            "# %s\n typedef struct {\n  double lambda;\n  double pass;\n } KFILTER;\n\n")
        for w, r in zip(wave, response):
            par.write("KFILTER %10.4f %11.7f\n" % (w, r))
        par.close()

    return wave, response
