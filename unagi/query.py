#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SQL search related functions"""

from .hsc import Hsc

__all__ = ['HELP_BASIC', 'COLUMNS_CONTAIN']

HELP_BASIC = "SELECT * FROM help('{0}');"

COLUMNS_CONTAIN = "SELECT * FROM help('{0}.%{1}%');"
