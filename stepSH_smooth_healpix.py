im

######### 0.0) Load the numpy and healpy packages and set settings:

import numpy as np
import healpy as hp

fwhm =  0.17  # Should be given in radians!

##### 0.1) Read in the HEALPix Maps Before starting the LOOP OF ALL REGIONS

print "Reading HEALPix Maps"
maps = ["im_iras12","im_iras25,"im_iras60","im_akari65","im_akari90","im_iras100","im_akari140","im_akari160","im_planck857","im_planck545"]
for i in range (0,Nmaps)
  map_in = hp.read_map(maps[i]".fits", nest = True)

print "Finished reading HEALPix Map" maps[i]

print "Masking bad pixels"

###### 0.2) We'll need to mask the bad pixels before feeding the maps to the visualization and smoothing routines

hp.pixelfunc.ma(map_in, badval=None, rtol=1e-05, atol=1e-08, copy=False)

###### 0.3) Take a quick look at the map before smoothing it:

hp.mollview(map_in, title='Histogram equalized Galactic', unit='MJy/sr', norm='hist', min=0,max=2000, xsize=2000)
hp.graticule()


####### 1.0) This is the function that actually does the smoothing
############# http://healpy.readthedocs.org/en/latest/generated/healpy.sphtfunc.smoothing.html#healpy.sphtfunc.smoothing
print "Smoothing HEALPix maps to" FWHM

for i in range(0,nmaps):
  map_in = maps[i]
  map_out = hp.sphtfunc.smoothing(map_in, fwhm = fwhm, arcmin = True)

###### 1.1) Take a look at the map after smoothing it:
hp.mollview(map_out, title='Histogram equalized Galactic, Smoothed', unit='MJy/sr', norm='hist', min=0,max=2000, xsize=2000)
hp.graticule()

###### 2.0) Now that the smoothing is finished, we'll write the smoothed maps to a new file:

hp.write_map(images[i]"_smoothed.fits", map_out)

print "Smoothed" images[i] "to FWHM" fwhm "radians"

