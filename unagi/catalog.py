#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions about using HSC catalogs."""

import numpy as np

__all__ = ['moments_to_shape']


def moments_to_shape(xx, yy, xy, axis_ratio=False, radian=False):
    """
    Convert the 2nd moments into elliptical shape: radius, ellipticity, position angle.
    """
    e1 = (xx - yy) / (xx + yy)
    e2 = (2.0 * xy / (xx + yy))
    rad = np.sqrt(xx + yy)
    ell = np.sqrt(e1 ** 2.0 + e2 ** 2.0)
    theta = (0.5 * np.arctan2(e2, e1))
    if not radian:
        theta *= (180.0 / np.pi)

    if axis_ratio:
        return rad, 1.0 - ell, theta

    return rad, ell, theta
