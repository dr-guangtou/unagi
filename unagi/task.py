#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from .target import SspObject
from .config import SspRerun

__all__ = ('SspDasTask', 'SspCoaddImage',)

HSC_SSP_PIXEL = 0.168  # arcsec/pixel


class SspDasTask(object):
    """Object for HSC SSP DAS related task.

    Examples
    --------

    The examples below illustrate common usage of the `SspDasTask` object.

        >>> from unagi.task import SspDasTask

    Parameters
    ----------
    target : unagi.target.SspObject
        Target of the task.
    rerun : unagi.config.SspDasConfig
        Object for SSP DAS rerun.
    """
    def __init__(self, target, rerun):
        # Test whether the inputs are the right type
        if not isinstance(target, SspObject):
            raise TypeError("Not the right type of target.")
        if not isinstance(rerun, SspRerun):
            raise TypeError("Not the right type of target.")


class SspCoaddImage(SspDasTask):
    """Task to download HSC SSP coadd image.

    Examples
    --------

    The examples below illustrate common usage of the `SspDasTask` object.

        >>> from unagi.task import SspCoaddImage

    Parameters
    ----------
    target : unagi.target.SspObject
        Target of the task.
    rerun : unagi.config.SspDasConfig
        Object for SSP DAS rerun.
    """
    def __init__(self, target, rerun, **kwargs):
        super(SspDasTask, self).__init__(target, rerun)
