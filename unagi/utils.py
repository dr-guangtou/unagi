#!/usr/bin/env python
# -*- coding: utf-8 -*-

import string
import random

__all__ = ['same_string', 'random_string']


def _passively_decode_string(a):
    try:
        return a.decode()
    except AttributeError:
        return a


def same_string(a, b):
    """Compare string in a Python2 and Python3-safe way.

    Shamelessly "borrowed" from Halotools by Andrew Hearin

    Parameters:
    -----------
    a: str
        String A to be compared.
    b: str
        String B to be compared
    """
    a = _passively_decode_string(a)
    b = _passively_decode_string(b)
    return a == b


def random_string(length=5, chars=string.ascii_uppercase + string.digits):
    """
    Random string generator.

    Based on:
    http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python

    Parameters:
    -----------
    length: int
        Length of the string. Default: 5
    chars: string object
        Types of characters allowed in the random string. Default: ASCII_Uppercase + Digits.
    """
    return ''.join(random.choice(chars) for _ in range(length))
