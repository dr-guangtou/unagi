# Unagi (ウナギ)

<img src="doc/unagi_logo.png" width="20%">

----

Besides being the most delicious fish in the world and one of the [signature Japanese dishes](https://en.wikipedia.org/wiki/Unagi), `unagi` can also help you navigate through the [public](https://hsc.mtk.nao.ac.jp/ssp/) or [nternal data release](https://hscdata.mtk.nao.ac.jp/hsc_ssp/) of the Hyper Suprime-Cam (HSC) Subaru Strategic Program (SSP; also known as the HSC survey) -- an ambitious multi-band deeep photometric survey using the awesome prime-focus camera on the 8.2m Subaru telescope. 

Also, `unagi` does not stand for anything because forced acronym is for psychopath. 

If you want to learn more about the real unagi, you can check out [this video on Youtube](https://www.youtube.com/watch?v=1sqLCUuMMfo).  And here is a vlog about a famous unagi place I have tried: [Obana 尾花](https://www.youtube.com/watch?v=N26pjkM_z4A).  It is pretty close to downtown Tokyo and few statioins away from Kavli-IPMU.

Applications
------------

- [Generate cutout coadded or warped HSC images](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_hsc_cutout.ipynb).
- [Generate 3-color RGB picture of a small HSC region.](https://github.com/dr-guangtou/unagi/blob/master/demo/demo_color_image.ipynb)
- Query and download HSC PSF model.
- Download HSC files (coadded `Patch` or source catalogs).


Installation
------------

Right now, please just put `unagi` under your `PYTHONPATH` environment variable.

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
