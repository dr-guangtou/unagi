#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to deal with HSC filters."""

import os 

import numpy as np

from astropy.table import Table

import unagi

__all__ = ['filters_to_kcorrect']

FILTER_DIR = os.path.join(os.path.dirname(unagi.__file__), 'data', 'filters')


def filters_to_kcorrect(filename):
    """
    Convert a filter response curve to the Sedpy/Kcorrect format.
    """
    curve_file = os.path.join(FILTER_DIR, filename)
    if not os.path.isfile(curve_file):
        raise IOError("# Cannot find the response curve file {}".format(curve_file))

    # Read in the .txt response curve 
    curve = Table.read(curve_file, format='ascii')
    wave, response = np.asarray(curve[curve.colnames[0]]), np.asarray(curve[curve.colnames[1]])

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
