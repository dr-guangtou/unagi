#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sky background related"""

import os
import warnings

import numpy as np

from astropy.table import Table

from scipy.stats import sigmaclip
from scipy.stats import binned_statistic_2d

from . import utils
from . import plotting

__all__ = ['SkyObjs', 'AperPhot', 'S18A_APER']


class AperPhot():
    """
    Class for aperture photometry in HSC.
    """
    PIX = 0.168 # arcsec/pixe

    def __init__(self, name, rad, rerun='s18a'):
        """Start a aperture photometry object."""
        self.aper_id = name
        self.name = "aper{0}".format(self.aper_id)
        self.r_pix = rad
        self.area_pix = np.pi * (rad ** 2.0)
        self.r_arcsec = rad * self.PIX
        self.area_arcsec = np.pi * (self.r_arcsec ** 2.0)
        self.rerun = rerun

        # Name of the columns for flux and flux error
        self.flux_col = self.flux(rerun=self.rerun)
        self.err_col = self.err(rerun=self.rerun)

    def flux(self, band=None, rerun='s18a'):
        """Aperture flux column name in S18A."""
        if rerun == 's18a':
            if band is not None:
                return "{0}_apertureflux_{1}_flux".format(band.strip(), self.aper_id)
            return "apertureflux_{0}_flux".format(self.aper_id)
        else:
            raise NotImplementedError("# Only S18A data are available.")

    def err(self, band=None, rerun='s18a'):
        """Aperture flux error column name in S18A."""
        if rerun == 's18a':
            if band is not None:
                return "{0}_{1}sigma".format(band.strip(), self.flux(rerun=rerun))
            return "{0}sigma".format(self.flux(rerun=rerun))
        else:
            raise NotImplementedError("# Only S18A data are available.")

# Aperture flux in S18A
S18A_APER_ID = ['10', '15', '20', '30', '40', '57', '84',
                '118', '168', '235']
S18A_APER_RAD = [3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]
S18A_APER = {}
for ii, rr in zip(S18A_APER_ID, S18A_APER_RAD):
    S18A_APER['aper{0}'.format(ii)] = AperPhot(ii, rr)

class SkyObjs():
    """
    Class for HSC sky objects.
    """
    # Convert the flux from erg/s/cm^2/Hz to HSC image value
    CGS_TO_IMG = 1.7378E30
    # Convert the flux from erg/s/cm^2/Hz to muJy
    CGS_TO_MUJY = 1.0E29
    # Convert the from from muJy to HSC image unit
    MUJY_TO_IMG = CGS_TO_IMG / CGS_TO_MUJY

    # List of filters
    FILTER_LIST = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']

    # Nicknames of filters
    FILTER_SHORT = ['g', 'r', 'i', 'z', 'y']

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
                self.skyobjs = skyobjs.as_array()
            except Exception:
                self.skyobjs = skyobjs.as_array()
        elif isinstance(skyobjs, np.ndarray) or isinstance(skyobjs, np.recarray):
            self.skyobjs = skyobjs

        # Minimum number of sky objects
        self.n_min = nobj_min

        # List of Tracts
        self.tract_list = list(np.unique(self.skyobjs['tract']))
        self.n_tract = len(self.tract_list)

        # List of Patches and Tracts
        self.tract_patch = np.unique(
            ["{0}_{1:03d}".format(t, p) for t, p in
             zip(self.skyobjs['tract'], self.skyobjs['patch'])])
        self.n_tract_patch = len(self.tract_patch)

    def select_tract(self, tract, patch=None, n_min=10, verbose=True) -> 'SkyObjs':
        """Select sky objects on one Tract (and Patch) from the catalog """
        tract_mask = self.skyobjs['tract'] == tract
        if tract_mask.sum() == 0:
            if verbose:
                warnings.warn("# Tract {0} is not available!".format(tract))
            return SkyObjs(self.skyobjs[self.skyobjs['tract'] < 0])

        if patch is not None:
            tract_mask = tract_mask & (self.skyobjs['patch'] == patch)
            if tract_mask.sum() == 0:
                if verbose:
                    warnings.warn(
                        "# Tract {0}-Patch {1} is not available!".format(tract, patch))
                return SkyObjs(self.skyobjs[self.skyobjs['tract'] < 0])

        # Number of sky objects available
        n_skyobj = tract_mask.sum()
        if n_skyobj <= n_min:
            if patch is None:
                if verbose:
                    warnings.warn("# Tract {0} has less than {1} skyobjs: {2}".format(
                        tract, n_min, n_skyobj))
            else:
                if verbose:
                    warnings.warn("# Tract {0}-Patch {1} has < {2} skyobjs: {3}".format(
                        tract, patch, n_min, n_skyobj))
            return SkyObjs(self.skyobjs[self.skyobjs['tract'] < 0])

        return SkyObjs(self.skyobjs[tract_mask])

    def select_box(self, ra1, ra2, dec1, dec2, n_min=5, verbose=True) -> 'SkyObjs':
        """Select sky objects in a box region."""
        # Order of the coordinates
        if ra1 >= ra2:
            ra1, ra2 = ra2, ra1
        if dec1 >= dec2:
            dec1, dec2 = dec2, dec1

        # Select sky objects in that region
        box_mask = ((self.skyobjs[self.ra_col] >= ra1) &
                    (self.skyobjs[self.ra_col] <= ra2) &
                    (self.skyobjs[self.dec_col] >= dec1) &
                    (self.skyobjs[self.dec_col] <= dec2))

        if box_mask.sum() == 0:
            if verbose:
                warnings.warn(
                    "# No sky object in this region: {0}:{1}-{2}:{3}".format(
                        ra1, ra2, dec1, dec2))
            return SkyObjs(self.skyobjs[self.skyobjs['tract'] < 0])

        if box_mask.sum() <= n_min:
            if verbose:
                warnings.warn("# Only find {0} sky object(s)".format(box_mask.sum()))

        return SkyObjs(self.skyobjs[box_mask])

    def select_circle(self, ra, dec, radius, n_min=5, verbose=True):
        """Select sky objects within a circle. Radius is in astropy.units."""
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        if str(radius).replace('.', '', 1).isdigit():
            radius = radius * u.arcsec

        c = SkyCoord(ra, dec, frame='icrs', unit='deg')
        catalog = SkyCoord(self.skyobjs[self.ra_col], 
                           self.skyobjs[self.dec_col], unit='deg', frame='icrs')
        circle_mask = (catalog.separation(c) < radius)
        
        if circle_mask.sum() == 0:
            if verbose:
                warnings.warn(
                    "# No sky object in this region: RA = {0}, DEC = {1}, r = {2} arcsec".format(
                        ra, dec, radius))
            return SkyObjs(self.skyobjs[self.skyobjs['tract'] < 0])

        if circle_mask.sum() <= n_min:
            if verbose:
                warnings.warn("# Only find {0} sky object(s)".format(circle_mask.sum()))

        return SkyObjs(self.skyobjs[circle_mask])

    def flux_stats(self, aper, band, rerun='s18a', sigma=3.5,
                   kde=False, bw=None, to_mujy=True, prefix=None):
        """Basic statistics of the flux."""
        u_factor = self.CGS_TO_MUJY if to_mujy else 1.0
        assert band in self.FILTER_SHORT, "# Wrong filter name: {}".format(band)

        flux_col = aper.flux(rerun=rerun, band=band)

        try:
            flux = self.skyobjs[flux_col] * u_factor
        except ValueError:
            raise Exception("# Wrong flux column name: {0}".format(flux_col))

        return utils.stats_summary(flux, sigma=sigma, n_min=self.n_min,
                                   kde=kde, bw=bw, prefix=prefix)

    def snr_stats(self, aper, band, rerun='s18a', sigma=3.5,
                  kde=False, bw=None, prefix=None):
        """Basic statistics of the S/N."""
        assert band in self.FILTER_SHORT, "# Wrong filter name: {}".format(band)

        flux_col = aper.flux(rerun=rerun, band=band)
        err_col = aper.err(rerun=rerun, band=band)

        try:
            snr = self.skyobjs[flux_col] / self.skyobjs[err_col]
        except ValueError:
            raise Exception("# Wrong column names: {0}/{1}".format(flux_col, err_col))

        return utils.stats_summary(snr, sigma=sigma, n_min=self.n_min,
                                   kde=kde, bw=bw, prefix=prefix)

    def mu_stats(self, aper, band, to_mujy=True, rerun='s18a', sigma=3.5,
                 kde=False, bw=None, prefix=None):
        """Basic statistics of the aperture flux density."""
        u_factor = self.CGS_TO_MUJY if to_mujy else 1.0
        assert band in self.FILTER_SHORT, "# Wrong filter name: {}".format(band)

        flux_col = aper.flux(rerun=rerun, band=band)

        try:
            mu = self.skyobjs[flux_col] * u_factor / aper.area_arcsec
        except ValueError:
            raise Exception("# Wrong flux column name: {0}".format(flux_col))

        return utils.stats_summary(mu, sigma=sigma, n_min=self.n_min,
                                   kde=kde, bw=bw, prefix=prefix)

    def sum_all_filters(self, aper, **kwargs):
        """Provide a summary of sky objects in all five bands."""
        aper_sum = {}
        for band in self.FILTER_SHORT:
            # Sky flux
            flux_pre = "{0}_{1}_flux".format(aper.name, band)
            flux_stats = self.flux_stats(aper, band, prefix=flux_pre, **kwargs)
            # S/N of sky flux
            snr_pre = "{0}_{1}_snr".format(aper.name, band)
            snr_stats = self.flux_stats(aper, band, prefix=snr_pre, **kwargs)
            # Surface flux density
            mu_pre = "{0}_{1}_mu".format(aper.name, band)
            mu_stats = self.flux_stats(aper, band, prefix=mu_pre, **kwargs)
            aper_sum = {**aper_sum, **flux_stats, **snr_stats, **mu_stats}

        return aper_sum

    def sum_aper_list(self, aper_list, **kwargs):
        """Summary of sky objects in all five bands for a list of apertures."""
        if isinstance(aper_list, list):
            return {key: value for stats in [
                self.sum_all_filters(aper, **kwargs) for aper in aper_list]
                    for key, value in stats.items()}
        else:
            raise TypeError("# Need a list of AperPhot objects!")

    def sum_all_tracts(self, aper_list, patch=False, verbose=True, **kwargs):
        """Provide summary for all the Tracts-(Patches) in the catalog."""
        result = []
        if not patch:
            for t in self.tract_list:
                if isinstance(aper_list, list):
                    t_sum = self.select_tract(t, verbose=verbose).sum_aper_list(
                        aper_list, **kwargs)
                    t_sum['tract'] = t
                    result.append(t_sum)
                elif isinstance(aper_list, AperPhot):
                    t_sum = self.select_tract(t, verbose=verbose).sum_all_filters(
                        aper_list, **kwargs)
                    t_sum['tract'] = t
                    result.append(t_sum)
        else:
            for t, p in [(int(tp.split('_')[0]), int(tp.split('_')[1]))
                         for tp in self.tract_patch]:
                if isinstance(aper_list, list):
                    t_sum = self.select_tract(t, patch=p, verbose=verbose).sum_aper_list(
                        aper_list, **kwargs)
                    t_sum['tract'] = t
                    t_sum['patch'] = p
                    result.append(t_sum)
                elif isinstance(aper_list, AperPhot):
                    t_sum = self.select_tract(t, patch=p, verbose=verbose).sum_all_filters(
                        aper_list, **kwargs)
                    t_sum['tract'] = t
                    t_sum['tract'] = p
                    result.append(t_sum)

        return result

    def get_summary(self, aper, band, prop, tract=None, patch=None,
                    rerun='s18a', kde=False, bw=0.2, sigma=3.0, to_mujy=True,
                    plot=False):
        """Show histogram of the properties of sky objects."""
        assert band in self.FILTER_SHORT, "# Wrong filter name: {}".format(band)
        u_factor = self.CGS_TO_MUJY if to_mujy else 1.0

        if tract is None:
            sky = self.skyobjs
        else:
            sky = self.select_tract(tract, patch=patch).skyobjs

        # Column names
        flux_col = aper.flux(rerun=rerun, band=band)
        err_col = aper.err(rerun=rerun, band=band)

        try:
            if prop == 'flux':
                values = sky[flux_col] * u_factor
            elif prop == 'snr':
                values = sky[flux_col] / sky[err_col]
            elif prop == 'mu':
                values = (sky[flux_col] * u_factor) / aper.area_arcsec
            else:
                raise Exception("# Wrong type of properties: flux/snr/mu")
        except ValueError:
            raise Exception("# Wrong flux column name: {0}".format(flux_col))

        clipped, summary = utils.stats_summary(
            values, sigma=sigma, n_min=self.n_min, kde=kde, bw=bw,
            return_clipped=True)

        if plot:
            if tract is None:
                region = None
            if tract is not None and patch is None:
                region = r'$\mathrm{Tract\ }%5d$' % tract
            elif tract is not None and patch is not None:
                region = r'${0}:{1}$'.format(tract, patch)

            aper_str = r"$\rm {0}$".format(aper.name[0].upper() + aper.name[1:])

            hist = plotting.plot_skyobj_hist(
                clipped, summary, band, prop, region=region, aper=aper_str, fontsize=20)

            return clipped, summary, hist

        return clipped, summary

    def plot_map(self, aper, band, prop, tract=None, patch=None, boxsize=0.19,
                 rerun='s18a', sigma=3.0, to_mujy=True, region=None, y_size=4,
                 margin=0.2, fontsize=30):
        """Show histogram of the properties of sky objects."""
        assert band in self.FILTER_SHORT, "# Wrong filter name: {}".format(band)
        u_factor = self.CGS_TO_MUJY if to_mujy else 1.0

        if tract is None:
            sky = self.skyobjs
        else:
            sky = self.select_tract(tract, patch=patch).skyobjs

        # Column names
        flux_col = aper.flux(rerun=rerun, band=band)
        err_col = aper.err(rerun=rerun, band=band)

        try:
            if prop == 'flux':
                values = sky[flux_col] * u_factor
            elif prop == 'snr':
                values = sky[flux_col] / sky[err_col]
            elif prop == 'mu':
                values = (sky[flux_col] * u_factor) / aper.area_arcsec
            else:
                raise Exception("# Wrong type of properties: flux/snr/mu")
        except ValueError:
            raise Exception("# Wrong flux column name: {0}".format(flux_col))

        flag = np.isfinite(values)
        values = values[flag]

        # RA, Dec
        ra, dec = sky['ra'][flag], sky['dec'][flag]

        # Number of bins
        x_bins = np.floor((np.max(ra) - np.min(ra)) / boxsize)
        y_bins = np.floor((np.max(dec) - np.min(dec)) / boxsize)

        _, low, upp = sigmaclip(values, low=sigma, high=sigma)
        mask = (values >= low) & (values <= upp)

        n_sky, x_edges, y_edges, _ = binned_statistic_2d(
            ra[mask], dec[mask], values[mask], 'count', bins=[x_bins, y_bins])

        mean_sky, _, _, _ = binned_statistic_2d(
            ra[mask], dec[mask], values[mask], 'mean', bins=[x_bins, y_bins])

        _, low_mean, upp_mean = sigmaclip(
            mean_sky[np.isfinite(mean_sky)].flatten(), low=sigma, high=sigma)

        v_edge = np.min(np.abs([low_mean, upp_mean]))

        if region is not None:
            region_str = r'$\rm {0}$'.format(region)
        else:
            region_str = ''

        band_str = r'$\ \ \ \rm {0}-band$'.format(band)
        aper_str = r"$\ \ \ \rm {0}$".format(aper.name[0].upper() + aper.name[1:])

        skyobj_map = plotting.map_skyobjs(
            x_edges, y_edges, n_sky, mean_sky,
            label=region_str + band_str + aper_str, n_min=10,
            vmin=-v_edge, vmax=v_edge, y_size=y_size, margin=margin, fontsize=fontsize)

        return x_edges, y_edges, n_sky, mean_sky, skyobj_map
