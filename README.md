# solaratlas
Data and programs to work with solar spectral atlases

See *atlasfiles* directory for description of available atlas files.

**atlastools.py** has the main procedures for listing atlas file contents, loading those atlases into memory, and plotting the atlas over a specified range

**Atlas_Plot_Example.ipynb** is a Jupyter notebook that loads a predetermined atlas and plots the spectrum over a given range.


*Simple usage example:*

`file  = 'atlasfiles/neckel.hamburg.atlas.wvscl_smooth.bintab.v3.fits.gz'`\
`atlas = atlastools.make_atlas(file, 1, loaddata=1, startwave=600*u.nm, endwave=660*u.nm)`\
`atlastools.atlas_spectrum_plot('Local Intensity   1', atlas, startwave=630.0*u.nm, endwave=630.4*u.nm, plot_unit='angstrom')`
`atlas_wavescale = atlas.sun.spectral_axis`\
`atlas_intensity = atlas.sun.flux`

