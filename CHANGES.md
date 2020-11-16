### v0.1.1 (2020-11-16)

- First release after updates to support HSC internal release `DR3`, which include new
  reruns for `S19A` and `S20A` releases. `DR1` becomes unavailable now.
- Fix a problem related to how `astropy.fits.open()` handle password protected files. Now
  pass the `urllib.request.urlopen(url)` object to it instead of the URL itself.
- Fix a few Python warnings including:
    - Changes `is not` to `!=`
    - Make a copy of certain `matplotlib` colormap before modifying it.
- Add a new demo notebook to compare the differences in coadd image between `DR2` and
  `DR3`.
