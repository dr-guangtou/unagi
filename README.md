# Unagi (ウナギ)

<img src="https://github.com/dr-guangtou/unagi/blob/master/doc/unagi_logo.png" width="20%">

----

Besides being the most delicious fish in the world and one of the [signature Japanese dishes](https://en.wikipedia.org/wiki/Unagi), `unagi` can also help you navigate through the [public](https://hsc.mtk.nao.ac.jp/ssp/) or [internal data release](https://hscdata.mtk.nao.ac.jp/hsc_ssp/) of the Hyper Suprime-Cam (HSC) Subaru Strategic Program (SSP; also known as the HSC survey) -- an ambitious multi-band deeep photometric survey using the awesome prime-focus camera on the 8.2m Subaru telescope.

Also, `unagi` does not stand for anything because forced acronym is for psychopath.

If you want to learn more about the real unagi, you can check out [this video on Youtube](https://www.youtube.com/watch?v=1sqLCUuMMfo).  And here is a vlog about a famous unagi place I have tried: [Obana 尾花](https://www.youtube.com/watch?v=N26pjkM_z4A).  It is pretty close to downtown Tokyo and few statioins away from Kavli-IPMU.

If you want to learn more about this amazing food, just try `import unagi; unagi.unagi()` in your `Jupyter` Notebook.


Recent Updates
--------------

- 06/03/2019: **Support for HSC PDR2 is added** (need more tests)
  - The schema of all tables in the `PDR2_DUD` and `PDR2_WIDE` reruns are available in `unagi/data/` folder as a `.json` file.
- 06/04/2019: **New features: filters.py and camera.py**
  - Help you get important information about HSC camera (e.g. CCD QE, primary mirror reflectivity) and filters (e.g. transmission curves, absolute AB magnitude of Sun in each filter etc.).
- 06/06/2019: **New features: add hsc_check_coverage() method** to check if a coordinate is covered by one HSC `Patch`. See [here](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_hsc_sql.ipynb) for example.

Applications
------------

- [Config the access to HSC database and get basic survey information](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_hsc_config.ipynb)

<img src="https://github.com/dr-guangtou/unagi/blob/master/doc/unagi_config.png" width="70%">

- [Access to information about HSC filters and camera](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_hsc_filters.ipynb)

<img src="https://github.com/dr-guangtou/unagi/blob/master/doc/unagi_filters_camera.png" width="70%">

- [Generate cutout coadded or warped HSC images](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_hsc_cutout.ipynb).

<img src="https://github.com/dr-guangtou/unagi/blob/master/doc/unagi_cutout.png" width="70%">

- [Generate 3-color RGB picture of a small HSC region.](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_color_image.ipynb)

<img src="https://github.com/dr-guangtou/unagi/blob/master/doc/unagi_color.png" width="70%">

- [Query and download HSC PSF model.](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_hsc_psf.ipynb)

<img src="https://github.com/dr-guangtou/unagi/blob/master/doc/unagi_psf.png" width="70%">

- [Dealing with the mask planes of HSC images](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_hsc_mask.ipynb)

<img src="https://github.com/dr-guangtou/unagi/blob/master/doc/unagi_mask.png" width="70%">

- [Basic SQL search of HSC database](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_hsc_sql.ipynb)

<img src="https://github.com/dr-guangtou/unagi/blob/master/doc/unagi_sql.png" width="70%">

- Using HSC sky object to characterize residual background


TODO List
---------

- Directly download data products from HSC pipeline (e.g., coadded `Patch` or source catalogs)

- Access to HSC weak-lensing shape catalog

- Access to HSC random catalog

- Reproduce the `CModel` result

Installation
------------

- `python setup.py install` or `python setup.py develop` will do the job.
- Right now, `unagi` only supports `Python>=3`.  If you are still using `Python 2`, you should make the switch.
- `unagi` only depends on `numpy`, `scipy`, `astropy`, and `matplotlib`. All can be installed using `pip` or `conda`.

Documents
---------

I <del>promise</del>hope that documents will be available soon...but right now, please take a look at the Jupyter Notebook [demos](https://github.com/dr-guangtou/unagi/tree/master/demo) for each functionality.


Acknowledgement
---------------

Thanks the HSC collaboration for making this amazing survey happen and make these beautiful data available.  Also thank the good people at NAOJ who work tirelessly to prepare the data release.


Reporting bugs
--------------

If you notice a bug in `unagi` (and you will~), please file an detailed issue at:

https://github.com/dr-guangtou/unagi/issues



Requesting features
-------------------

If you would like to request a new feature, do the same thing.


License
-------

Copyright 2019 Song Huang and contributors.

`unagi` is free software made available under the MIT License. For details see
the LICENSE file.
