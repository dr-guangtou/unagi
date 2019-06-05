#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Deal with the mask plane of HSC images."""

import copy

import numpy as np

from matplotlib import colors

from astropy.table import Table

from scipy.ndimage import gaussian_filter

from . import plotting

__all__ = ['Mask', 'BitMasks',
           'MASK_CMAP', 'S18A_BITMASKS', 'PDR1_BITMASKS']

MASK_CMAP = plotting.random_cmap(512, background_color='white')

class BitMasks():
    """
    Class for defining HSC bitmasks.
    """

    def __init__(self, data_release='s18a'):
        """
        Define the bitmask plane of HSC coadd image.
        """
        if data_release == 's18a':
            self._bitmasks = S18A_BITMASKS
            self.n_bits = 32
            self.type = np.uint32
        elif data_release == 'pdr1':
            self._bitmasks = PDR1_BITMASKS
            self.n_bits = 16
            self.type = np.uint16
        elif data_release == 'pdr2':
            self._bitmasks = PDR2_BITMASKS
            self.n_bits = 32
            self.type = np.uint32
        else:
            raise NotImplementedError(
                "# Rerun {0} is not available yet".format(data_release))

    @property
    def bitmasks(self):
        """
        Structured array that describe the mask planes.
        """
        return self._bitmasks

    @bitmasks.setter
    def bitmasks(self, mask_array):
        self._bitmasks = mask_array

    # Number of bits used
    @property
    def n_mask(self):
        """
        Number of mask planes used.
        """
        return len(self.bitmasks)

    def bits2name(self, idx):
        """
        Convert the bit value to the name of the mask plane.
        """
        return self.bitmasks[self.bitmasks['bits'] == idx]['name'][0]

    def check(self, name_or_bit):
        """
        Check whether a name or a bit is in this "library".
        """
        if isinstance(name_or_bit, str):
            return name_or_bit in self.bitmasks['name']
        else:
            return name_or_bit in self.bitmasks['bits']

    def name2bits(self, name):
        """
        Convert the name of the mask plane to bits value.
        """
        return self.bitmasks[self.bitmasks['name'] == name.strip().upper()]['bits'][0]

    def get_index(self, name_or_bit):
        """
        Get the index of the mask plane.
        """
        if not self.check(name_or_bit):
            raise NameError("Mask {} in not available".format(name_or_bit))
        if isinstance(name_or_bit, str):
            flag = self.bitmasks['name'] == name_or_bit.strip().upper()
        else:
            flag = self.bitmasks['bits'] == name_or_bit

        return np.where(flag)[0][0]

    def get_color(self, name_or_bit):
        """
        Return the color string used for displaying the mask plane.
        """
        return self.bitmasks[self.get_index(name_or_bit)]['color']

    def to_table(self):
        """
        Convert the bitmask structured array into astropy table.
        """
        return Table(self.bitmasks)

    def show_table(self):
        """
        Display the bitmask table.
        """
        try:
            _ = get_ipython
            return Table(self.bitmasks).show_in_notebook()
        except NameError:
            return Table(self.bitmasks).show_in_browser()


class Mask():
    """
    Class for HSC mask plane.
    """

    def __init__(self, mask, wcs=None, data_release='s18a', verbose=False):
        """
        Initialize a HSC mask plane object.

        Parameters:
        -----------
        mask: numpy.ndarray
            2-D bitmask plane from HSC image.
        """
        # WCS information
        self.wcs = wcs

        # Table for bitmask planes
        self.data_release = data_release

        # BitMask object for this mask
        self.bitmasks = BitMasks(data_release=self.data_release)

        # Key information about the mask planes
        self.library = self.bitmasks.bitmasks
        self.names = list(self.library['name'])
        self.bits = list(self.library['bits'])
        self.n_mask = self.bitmasks.n_mask
        self.type = self.bitmasks.type

        # Decode the bitmask array
        self._masks = mask.astype(self.bitmasks.type)

    @property
    def masks(self):
        """
        Decoded mask planes in a 3-D numpy array.
        """
        return self._masks

    @masks.setter
    def masks(self, mask_array):
        self._masks = mask_array

    def get_cmap(self, bit_list):
        """
        Get the colormap to show the mask plane.
        """
        if not np.all(self.check(bit_list)):
            raise NameError("One or more mask planes are not available!")

        if not isinstance(bit_list, list):
            cmap = colors.ListedColormap(
                ['white', self.bitmasks.get_color(bit_list)])
            cmap.set_under(color='w', alpha=0.0)
            return cmap
        else:
            cmap_list = []
            for b in bit_list:
                cmap = colors.ListedColormap(
                    ['white', self.bitmasks.get_color(b)])
                cmap.set_under(color='w', alpha=0.0)
                cmap_list.append(cmap)
            return cmap_list

    def display(self, bit_list, alpha_list=None, **kwargs):
        """
        Display one or multiple layers of masks.
        """
        mask = self.extract(bit_list, show=True)
        cmap = self.get_cmap(bit_list)
        if not isinstance(bit_list, list):
            return plotting.display_single(mask, scale='linear', cmap=cmap, alpha=0.9, **kwargs)
        else:
            return plotting.overplot_all(
                mask, xsize=6, ysize=6, stretch='linear', scale='minmax',
                alpha=0.7, alpha_list=alpha_list, cmap_list=cmap, vmin=1, **kwargs)

    def check(self, bit_list):
        """
        Check whether the bit or name is available.
        """
        if not isinstance(bit_list, list):
            return self.bitmasks.check(bit_list)
        else:
            return [self.bitmasks.check(b) for b in bit_list]

    def extract(self, bit_list, show=False):
        """
        Get the 2-D array of one or multiple mask plane.
        """
        if not np.all(self.check(bit_list)):
            raise NameError("One or more mask planes are not available!")

        if not isinstance(bit_list, list):
            if show:
                return self.masks & self.library[self.bitmasks.get_index(bit_list)]['value']
            return (
                self.masks & self.library[
                    self.bitmasks.get_index(bit_list)]['value'] > 0).astype(np.uint8)
        else:
            if show:
                return [
                    self.masks & self.library[self.bitmasks.get_index(b)]['value']
                    for b in bit_list]
            return [
                (self.masks & self.library[
                    self.bitmasks.get_index(b)]['value'] > 0).astype(np.uint8)
                for b in bit_list]

    def enlarge(self, name_or_bit, sigma=2.0, threshold=0.02):
        """
        Get an enlarged version of certain mask plane.
        """
        return (gaussian_filter(
            self.extract(name_or_bit, show=True), sigma=sigma) >= threshold).astype(np.uint8)

    def combine(self, name_or_bit_list):
        """
        Combine a few mask planes together.
        """
        if not isinstance(name_or_bit_list, list):
            raise TypeError("# Need to be a list of bitmask name or index")
        if len(name_or_bit_list) == 1:
            return self.extract(name_or_bit_list[0])
        return np.bitwise_or.reduce([self.extract(nb) for nb in name_or_bit_list])

    def clean(self, name_or_bit_list):
        """
        Remove one or a list of mask plane.
        """
        if not isinstance(name_or_bit_list, list):
            name_or_bit_list = [name_or_bit_list]

        if isinstance(name_or_bit_list[0], str):
            mask_use = copy.deepcopy(self.names)
        else:
            mask_use = copy.deepcopy(self.bits)

        _ = [mask_use.remove(b) for b in name_or_bit_list]

        return self.combine(mask_use)

    def mask_used(self):
        """
        Identify the mask plane that are actually used.
        """
        return [marr.sum() > 0 for marr in self.extract(self.names)]

    def name_used(self):
        """
        Get the names of the used mask planes.
        """
        return list(self.library[self.mask_used()]['name'])

    def show_used(self):
        """
        Display all mask planes that are used.
        """
        names = self.name_used()
        masks = self.extract(names, show=True)
        cmaps = self.get_cmap(names)

        labels = [r"$\rm {0}$".format(n.replace('_', '\_')) for n in names]

        fig = plotting.display_all(
            masks, n_column=4, img_size=3., label_list=labels, cmap_list=cmaps,
            label_x=0.1, label_y=0.9, fontsize=15, fontcolor='k',
            scale='minmax', scale_bar_color='k', scale_bar_y_offset=1.2,
            scale_bar_fontsize=15)

        return fig, masks


S18A_BITMASKS = np.array(
    [(0, 'BAD', 'Bad pixel',
      'red', 2 ** 0),
     (1, 'SAT', 'Source footprint includes saturated pixels',
      'tab:purple', 2 ** 1),
     (2, 'INTRP', 'Source footprint includes interpolated pixels',
      'tab:orange', 2 ** 2),
     (3, 'CR', 'Source footprint includes suspected CR pixels',
      'tab:pink', 2 ** 3),
     (4, 'EDGE', 'Source is close to the edge of the CCD',
      'tab:olive', 2 ** 4),
     (5, 'DETECTED', 'Pixel with detection above the threshold',
      'tab:blue', 2 ** 5),
     (6, 'DETECTED_NEGATIVE', 'Pixel in footprint that is detected as a negative object',
      'gray', 2 ** 6),
     (7, 'SUSPECT', 'Source footprint includes suspect pixels',
      'orangered', 2 ** 7),
     (8, 'NO_DATA', 'No useful data',
      'black', 2 ** 8),
     (9, 'BRIGHT_OBJECT', 'Bright star mask',
      'tab:brown', 2 ** 9),
     (10, 'CROSSTALK', 'Crosstalk',
      'rosybrown', 2 ** 10),
     (11, 'NOT_DEBLENDED', 'Pixel in footprint that is too large to deblend',
      'teal', 2 ** 11),
     (12, 'UNMASKEDNAN', 'NaN pixels that are interpolated over',
      'darkgreen', 2 ** 12),
     (13, 'REJECTED', 'Rejected due to a mask other than EDGE, NO_DATA, or CLIPPED',
      'violet', 2 ** 13),
     (14, 'CLIPPED', 'Pixel that has been clipped',
      'crimson', 2 ** 14),
     (15, 'SENSOR_EDGE', 'Pixel close to the edge of CCD sensor',
      'tab:olive', 2 ** 15),
     (16, 'INEXACT_PSF', 'PSF is not correct',
      'gold', 2 ** 16)],
    dtype=[('bits', np.uint8), ('name', '<U18'),
           ('meaning', '<U80'), ('color', '<U12'), ('value', np.uint32)])

PDR1_BITMASKS = np.array(
    [(0, 'BAD', 'Bad pixel',
      'red', 2 ** 0),
     (1, 'SAT', 'Source footprint includes saturated pixels',
      'tab:purple', 2 ** 1),
     (2, 'INTRP', 'Source footprint includes interpolated pixels',
      'tab:orange', 2 ** 2),
     (3, 'CR', 'Source footprint includes suspected CR pixels',
      'tab:pink', 2 ** 3),
     (4, 'EDGE', 'Source is close to the edge of the CCD',
      'tab:olive', 2 ** 4),
     (5, 'DETECTED', 'Pixel with detection above the threshold',
      'tab:blue', 2 ** 5),
     (6, 'DETECTED_NEGATIVE', 'Pixel in footprint that is detected as a negative object',
      'gray', 2 ** 6),
     (7, 'SUSPECT', 'Source footprint includes suspect pixels',
      'orangered', 2 ** 7),
     (8, 'NO_DATA', 'No useful data',
      'black', 2 ** 8),
     (9, 'BRIGHT_OBJECT', 'Bright star mask',
      'tab:brown', 2 ** 9),
     (10, 'CROSSTALK', 'Crosstalk',
      'rosybrown', 2 ** 10),
     (11, 'NOT_DEBLENDED', 'Pixel in footprint that is too large to deblend',
      'teal', 2 ** 11),
     (12, 'UNMASKEDNAN', 'NaN pixels that are interpolated over',
      'darkgreen', 2 ** 12),
     (13, 'REJECTED', 'Rejected due to a mask other than EDGE, NO_DATA, or CLIPPED',
      'violet', 2 ** 13)],
    dtype=[('bits', np.uint8), ('name', '<U18'),
           ('meaning', '<U80'), ('color', '<U12'), ('value', np.uint32)])

PDR2_BITMASKS = np.array(
    [(0, 'BAD', 'Bad pixel',
      'red', 2 ** 0),
     (1, 'SAT', 'Source footprint includes saturated pixels',
      'tab:purple', 2 ** 1),
     (2, 'INTRP', 'Source footprint includes interpolated pixels',
      'tab:orange', 2 ** 2),
     (3, 'CR', 'Source footprint includes suspected CR pixels',
      'tab:pink', 2 ** 3),
     (4, 'EDGE', 'Source is close to the edge of the CCD',
      'tab:olive', 2 ** 4),
     (5, 'DETECTED', 'Pixel with detection above the threshold',
      'tab:blue', 2 ** 5),
     (6, 'DETECTED_NEGATIVE', 'Pixel in footprint that is detected as a negative object',
      'gray', 2 ** 6),
     (7, 'SUSPECT', 'Source footprint includes suspect pixels',
      'orangered', 2 ** 7),
     (8, 'NO_DATA', 'No useful data',
      'black', 2 ** 8),
     (9, 'BRIGHT_OBJECT', 'Bright star mask',
      'tab:brown', 2 ** 9),
     (10, 'CROSSTALK', 'Crosstalk',
      'rosybrown', 2 ** 10),
     (11, 'NOT_DEBLENDED', 'Pixel in footprint that is too large to deblend',
      'teal', 2 ** 11),
     (12, 'UNMASKEDNAN', 'NaN pixels that are interpolated over',
      'darkgreen', 2 ** 12),
     (13, 'REJECTED', 'Rejected due to a mask other than EDGE, NO_DATA, or CLIPPED',
      'violet', 2 ** 13),
     (14, 'CLIPPED', 'Pixel that has been clipped',
      'crimson', 2 ** 14),
     (15, 'SENSOR_EDGE', 'Pixel close to the edge of CCD sensor',
      'tab:olive', 2 ** 15),
     (16, 'INEXACT_PSF', 'PSF is not correct',
      'gold', 2 ** 16)],
    dtype=[('bits', np.uint8), ('name', '<U18'),
           ('meaning', '<U80'), ('color', '<U12'), ('value', np.uint32)])
