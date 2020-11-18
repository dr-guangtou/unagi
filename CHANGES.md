### v0.1.3

- 2020-11-17: Small bug fixed; allow users to use the `coadd/bg` product (the `deepCoadd`
  files before the local sky subtraction applied) in `DR3` to generate cutout images and
  color pictures.

### v0.1.2 (2020-11-16)

- First release after updates to support HSC internal release `DR3`, which include new
  reruns for `S19A` and `S20A` releases. `DR1` becomes unavailable now.
- Fix a problem related to how `astropy.fits.open()` handle password protected files. Now
  pass the `urllib.request.urlopen(url)` object to it instead of the URL itself.
- Fix a few Python warnings including:
    - Changes `is not` to `!=`
    - Make a copy of certain `matplotlib` colormap before modifying it.
- Fix a small typo related to image cutout (by Christopher Bradshaw)
- Implementing version management by `setuptools_scm` here, and adding a travis CI script 
  to just check that all dependencies install correctly, and upload the package to `PyPi` 
  in case of a new version (by Francois Lanusse)
- Add a new demo notebook to compare the differences in coadd image between `DR2` and
  `DR3`.
