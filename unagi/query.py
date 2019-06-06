#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SQL search related functions"""

__all__ = ['HELP_BASIC', 'COLUMNS_CONTAIN', 'TABLE_SCHEMA', 'PATCH_CONTAIN']

HELP_BASIC = "SELECT * FROM help('{0}');"

TABLE_SCHEMA = "SELECT * FROM help('{0}.{1}');"

COLUMNS_CONTAIN = "SELECT * FROM help('{0}.%{1}%');"

PATCH_CONTAIN = """
    --- Find coadded patch images
    SELECT
        mosaic.tract,
        mosaic.patch,
        mosaic.filter01
    FROM
        {0}.mosaic JOIN public.skymap USING (skymap_id)
    WHERE
        patch_contains(patch_area, wcs, {1}, {2})
    ;
    """

def box_search_template(template, dr='pdr2', rerun='pdr2_wide'):
    """
    Get the SQL template for box search.
    """
    pass

def cone_search_template(template, dr='pdr2', rerun='pdr2_wide'):
    """
    Get the SQL template for box search.
    """
    pass