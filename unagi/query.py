#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SQL search related functions"""

__all__ = ['HELP_BASIC', 'COLUMNS_CONTAIN', 'TABLE_SCHEMA']

HELP_BASIC = "SELECT * FROM help('{0}');"

TABLE_SCHEMA = "SELECT * FROM help('{0}.{1}');"

COLUMNS_CONTAIN = "SELECT * FROM help('{0}.%{1}%');"
