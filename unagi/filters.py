#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to deal with HSC filters."""

import os
import json

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy import constants as const

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams

import unagi

__all__ = ['filters_to_kcorrect', 'hsc_filters', 'HscFilter', 'SolarSpectrum']

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
# Directory that keeps the solar spectra
SOLAR_DIR = os.path.join(os.path.dirname(unagi.__file__), 'data', 'solar')

# Formal names and nick names of the filters
FILTER_LIST = ['HSC-G', 'HSC-R', 'HSC-R2', 'HSC-I', 'HSC-I2', 'HSC-Z', 'HSC-Y',
               'HSC-NB387', 'HSC-NB816', 'HSC-NB921']
FILTER_SHORT = ['g', 'r', 'r2', 'i', 'i2', 'z', 'y', 'nb387', 'nb816', 'nb921']

# Colors to show the filter
FILTER_COLORS = ['#2ca02c', '#ff7f0e', '#ff7f0e', '#d62728', '#d62728',
                 '#8c564b', '#7f7f7f', 'c', 'm', 'purple']

# AB reference flux in unit of erg/s/cm^2/Hz
AB_FLUX = 3.631e-20
# Speed of light in unit of AA/s
C_AA_PER_SEC = const.c.to('AA/s').value


class Filter(object):
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

        # Color to show the filter
        self.color = FILTER_COLORS[self.index]

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

        # AB Zeropoint in counts
        self.ab_zero_counts = self._counts(
            self.wave, AB_FLUX * C_AA_PER_SEC / (self.wave ** 2))

        # AB absolute magnitude of Sun in this band
        self.solar_ab_mag = self._solar_ab_mag(kind='Willmer2018')

        # Convert to Kcorrect filter format
        self.to_kfilter()

        # Save it as a JSON file
        self.to_json()

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
        print("# AB absolute magnitude of the sun : {:6.3f} mag".format(self.solar_ab_mag))

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

    def _counts(self, obj_wave, obj_flux):
        """Project source spectrum onto filter and return the detector signal.
        """
        # Interpolate filter transmission to source spectrum
        newtrans = np.interp(obj_wave, self.wave, self.trans, left=0., right=0.)

        # Integrate lambda * f_lambda * R
        if True in (newtrans > 0.):
            positive = np.where(newtrans > 0.)[0]
            ind = slice(max(positive.min() - 1, 0),
                        min(positive.max() + 2, len(obj_wave)))
            counts = np.trapz(obj_wave[ind] * newtrans[ind] *
                              obj_flux[..., ind], obj_wave[ind], axis=-1)
            return np.squeeze(counts)

        return float('NaN')

    def _solar_ab_mag(self, kind='Willmer2018'):
        """Get the absolute magnitude of Sun in this band.
        """
        # Get the solar spectrum
        solar_spectrum = SolarSpectrum(kind=kind)
        counts = self._counts(solar_spectrum.wave, solar_spectrum.flux)
        return -2.5 * np.log10(counts / self.ab_zero_counts)

    def to_kfilter(self):
        """Convert the filter transmission curve into Kcorrect format.
        """
        _ = filters_to_kcorrect(self.filename)

    def as_dict(self):
        """Convert it into a dict."""
        filter_dict = self.__dict__
        filter_dict['wave'] = list(filter_dict['wave'])
        filter_dict['trans'] = list(filter_dict['trans'])
        return filter_dict

    def to_json(self):
        """Save as a JSON file."""
        pre, _ = os.path.splitext(self.filename)
        json_out = pre + '.json'
        with open(json_out, 'w') as jf:
            json.dump(self.as_dict(), jf)


def hsc_filters(origin=False, center=False, use_saved=True):
    """Get the tabel that summarizes information about HSC filters."""
    if origin:
        if center:
            tab_id = 'origin_center'
            xml_table = os.path.join(FILTER_DIR, 'hsc_filters_origin_center.xml')
        else:
            tab_id = 'origin_weighted'
            xml_table = os.path.join(FILTER_DIR, 'hsc_filters_origin_weighted.xml')
    else:
        tab_id = 'total'
        xml_table = os.path.join(FILTER_DIR, 'hsc_filters_total.xml')

    # If already created, just read it in
    if os.path.isfile(xml_table) and use_saved:
        return Table.read(xml_table, format='votable', table_id=tab_id)
    else:
        # Summarize the filters
        # Total transmission curves
        f_table = Table(
            [Filter(f, origin=origin, center=center).as_dict() for f in FILTER_LIST])
        f_table.write(xml_table, format='votable', table_id=tab_id, overwrite=True)
        return f_table


def filters_to_kcorrect(curve_file, verbose=False):
    """
    Convert a filter response curve to the Kcorrect format.

    This is used by Kcorrect and iSEDFit.
    """
    if not os.path.isfile(curve_file):
        raise IOError("# Cannot find the response curve file {}".format(curve_file))

    # Read in the .txt response curve
    wave, response = np.genfromtxt(curve_file, usecols=(0, 1), unpack=True)

    # Output file name
    prefix, _ = os.path.splitext(curve_file)
    output_par = prefix + '.par'

    if os.path.isfile(output_par):
        if verbose:
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


class SolarSpectrum(object):
    """Class to hold solar spectrum.

    On deriving Solar absolute magnitude:
    ------------------------------------
    To derive the Sun's absolute magnitudes, the IAU 2012 definitions of the astronomical unit
    (au) (Prša et al. 2016) and parsec were used, giving a distance modulus for the
    Sun of −31.5721 mag.

    To rationalize the use of solar constants, the IAU in 2015 adopted a nominal value for the Sun's
    luminosity L⊙ = 3.828 × 10^8 W (Prša et al. 2016), which corresponds to an average
    TSI of 1361 W m−2 at 1 au and an absolute bolometric magnitude of MBol = 4.74.
    """
    def __init__(self, kind='Willmer2018'):
        """Read in the solar spectrum and process it."""
        # Conversion to d=10 pc from 1 AU
        au_to_10pc = (1.0 / (3600 * 180 / np.pi * 10)) ** 2.0

        if kind.strip() == 'Willmer2018':
            # This is based on the work by Willmer 2018:
            #   https://iopscience.iop.org/article/10.3847/1538-4365/aabfdf
            # More information can be found on this webpage:
            #   http://mips.as.arizona.edu/~cnaw/sun.html
            # About the spectrum:
            #   The solar SED used here also combines observations with model spectra.
            #   The observed spectrum is a composite calculated by Haberreiter et al. (2017)
            #   using data from over 20 space-based instruments for an arbitrary date
            #   (2008 December 19, JDN = 2454820) during the solar minimum.
            #   The absolute calibration is set using the ${ATLAS}\,3$ composite spectrum
            #   of Thuillier et al. (2004), and constraining the Total Solar Irradiance (TSI)
            #   to the value measured for each day by Dudok de Wit et al. (2017).
            #   The observed composite ends at ~2.0 μm, and to extend the SED into the infrared;
            #   the model spectra of Fontenla et al. (2011) and Kurucz (2011) are used.
            self.file = os.path.join(SOLAR_DIR, 'sun_composite.fits')
            if os.path.isfile(self.file):
                solar_fits = fits.open(self.file)[1].data
                # Wavelength
                self.wave = solar_fits['WAVE']
                self.wave_unit = 'Angstrom'
                self.flux = solar_fits['FLUX'] * au_to_10pc
                self.flux_unit = 'erg/s/cm^2/AA'
            else:
                raise IOError("# Cannot find the solar spectrum {}".format(self.file))
        elif kind.strip() == 'Kurucz1993':
            # This is the theoretical spectrum of the Sun from Kurucz
            # from: ftp://ftp.stsci.edu/cdbs/grid/k93models/standards/sun_kurucz93.fits
            # The theoretical spectrum is scaled to match the observed spectrum
            # from 1.5 - 2.5 microns, and then it is used where the observed spectrum ends.
            # The theoretical model of the Sun from Kurucz93 atlas using the following
            # parameters when the Sun is at 1 au.
            #   log_Z T_eff log_g V_{Johnson} +0.0 5777 +4.44 -26.75
            self.file = os.path.join(SOLAR_DIR, 'sun_kurucz93.fits')

            # This file should be in AA and erg/s/cm^2/AA at 1AU
            if os.path.isfile(self.file):
                solar_fits = fits.open(self.file)[1].data
                # Wavelength
                self.wave = solar_fits['WAVELENGTH']
                self.wave_unit = 'Angstrom'
                self.flux = solar_fits['FLUX'] * au_to_10pc
                self.flux_unit = 'erg/s/cm^2/AA'
            else:
                raise IOError("# Cannot find the solar spectrum {}".format(self.file))
        else:
            raise NameError("# Wrong type of Solar spectrum: [Willmer2018|Kurucz1993]")

        self.wave_min, self.wave_max = self.wave.min(), self.wave.max()
        self.solar = np.column_stack((self.wave, self.flux))
