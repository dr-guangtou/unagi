#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

# Update: 2017-11-22
HSC_FILTERS = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y',
               'NB0387', 'NB0816', 'NB0921']

HSC_FSHORT = ['g', 'r', 'i', 'z', 'y', 'nb0387', 'nb0816', 'nb0921']

HSC_IDR_RERUN = ['s15b_udeep', 's15b_deep', 's15b_wide',
                 's16a_udeep', 's16a_deep', 's16a_wide', 's16a_wide2',
                 's17a_dud', 's17a_wide', 'any']
HSC_PDR_RERUN = ['pdr1_udeep', 'pdr1_deep', 'pdr1_wide', 'any']

HSC_DR1_URL = "https://hscdata.mtk.nao.ac.jp/das_quarry/dr1/cgi-bin/"
HSC_DR2_URL = "https://hscdata.mtk.nao.ac.jp/das_quarry/dr2/cgi-bin/"
HSC_PDR_URL = "https://hsc-release.mtk.nao.ac.jp/das_quarry/cgi-bin/"

HSC_DR1_PSF = "https://hscdata.mtk.nao.ac.jp/psf/4/cgi"
HSC_DR2_PSF = "https://hscdata.mtk.nao.ac.jp/psf/5/cgi"
HSC_PDR_PSF = "https://hsc-release.mtk.nao.ac.jp/psf/pdr1/cgi/"
