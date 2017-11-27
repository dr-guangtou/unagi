#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


__all__ = ('same_string', )


def _passively_decode_string(a):
    try:
        return a.decode()
    except AttributeError:
        return a


def same_string(a, b):
    """Compare string in a Python2 and Python3-safe way.

    Shamelessly "borrowed" from Halotools by
    """
    a = _passively_decode_string(a)
    b = _passively_decode_string(b)
    return a == b
