#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Get key information about the camera."""

import os

import numpy as np

import unagi

__all__ = []

# Directory that save the files related to camera
CAMERA_DIR = os.path.join(os.path.dirname(unagi.__file__), 'data', 'camera')

# List of primary mirror data available
PRIMARY_MIRROR_LIST = ['']

class Camera(object):
    """Class that keeps information about HSC camera."""

    def __init__(self, qe_file='qe_ccd_HSC.txt', dewar_file='throughput_win.txt',
                 popt2_file='throughput_popt2.txt', vignetting_file='vignetting.txt',
                 mirror_date='20190426-2'):
        """Gather basic information."""
        # Quantum efficiency of FDCCD
        self.qe_file = qe_file
        self.qe = self.get_qe()
        # Transmittance of the dewar window.
        self.dewar_file = dewar_file
        self.dewar = self.get_dewar()
        # Transmittance of the Primary Focus Unit of the HSC (POpt2).
        self.popt2_file = popt2_file
        self.popt2 = self.get_popt2()
        # Vignetting
        self.vignet_file = vignetting_file
        self.vignet = self.get_vignetting()
        # Reflectivity of the Subaru Telescope's Primary Mirror
        self.mirror_list = self.get_primary_list()
        if mirror_date.strip() not in self.mirror_list:
            raise NameError("# Wrong choice of date for primary mirror reflectivity data.")
        self.mirror_date = mirror_date
        self.mirror_file = 'Subaru_M1_R_{}.txt'.format(self.mirror_date)
        self.primary_reflect = self.get_primary_mirror()

    def get_qe(self, filename=None):
        """Quantum efficiency of FDCCD.

        Please see: https://www.subarutelescope.org/Observing/Instruments/HSC/sensitivity.html
        """
        if filename is None:
            filename = self.qe_file
        wave, qe = np.genfromtxt(os.path.join(CAMERA_DIR, filename),
                                 usecols=(0, 1), unpack=True)

        return np.column_stack((wave, qe))

    def get_dewar(self, filename=None):
        """Transmittance of the dewar window.

        Please see: https://www.subarutelescope.org/Observing/Instruments/HSC/sensitivity.html
        """
        if filename is None:
            filename = self.dewar_file
        wave, dewar = np.genfromtxt(os.path.join(CAMERA_DIR, filename),
                                    usecols=(0, 1), unpack=True)

        return np.column_stack((wave, dewar))

    def get_popt2(self, filename=None):
        """Transmittance of the Primary Focus Unit of the HSC (POpt2).

        Please see: https://www.subarutelescope.org/Observing/Instruments/HSC/sensitivity.html
        """
        if filename is None:
            filename = self.popt2_file
        wave, popt2 = np.genfromtxt(os.path.join(CAMERA_DIR, filename),
                                    usecols=(0, 1), unpack=True)

        return np.column_stack((wave, popt2))

    def get_vignetting(self, filename=None):
        """Vignetting."""
        if filename is None:
            filename = self.vignet_file
        rad, vignet = np.genfromtxt(os.path.join(CAMERA_DIR, filename),
                                    usecols=(0, 1), unpack=True)

        return np.column_stack((rad, vignet))

    def get_primary_list(self):
        """List of dates that primary mirror reflectivity is available."""
        with open(os.path.join(CAMERA_DIR, 'primary_mirror.lis')) as f:
            date_list = f.read().splitlines()
        return date_list

    def get_primary_mirror(self, filename=None):
        """Vignetting."""
        if filename is None:
            filename = self.mirror_file
        wave, pos3_1, pos3_2, pos3_3 = np.genfromtxt(
            os.path.join(CAMERA_DIR, filename), usecols=(0, 1, 2, 3), unpack=True)

        return np.column_stack(
            (wave, pos3_1 / 100., pos3_2 / 100., pos3_3 / 100.))
