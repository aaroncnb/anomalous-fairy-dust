

######### 0.0) Load the numpy and healpy packages and set settings:

import numpy as np
import healpy as hp

fwhm =  0.0174532925 # Should be given in radians!

########0.1) Read in the HEALPix Maps Before starting the LOOP OF ALL REGIONS

print "Reading HEALPix Maps"
maps = ["im_iras12","im_iras25,"im_iras60","im_akari65","im_akari90","im_iras100","im_akari140","im_akari160","im_planck857","im_planck545"]
for i in range
map12 = hp.read_map("../Data/im_iras12.fits", nest = True)
map25 = hp.read_map("../Data/im_iras25.fits", nest = True)
map60 = hp.read_map("../Data/im_iras60.fits", nest = True)
map100 = hp.read_map("../Data/healpix10/im_akari65.fits", nest = True)
map90 = hp.read_map("../Data/healpix10/im_akari90.fits", nest = True)
map140 = hp.read_map("../Data/healpix10/im_akari140.fits", nest = True)
map160 = hp.read_map("../Data/healpix10/im_akari160.fits", nest = True)
map857 = hp.read_map("../Data/im_planck857.fits", nest = True)
map545 = hp.read_map("../Data/im_planck545.fits", nest = True)

print "Finished reading HEALPix Maps"

####### 1.0) This is the function that actually does the smoothing
############# http://healpy.readthedocs.org/en/latest/generated/healpy.sphtfunc.smoothing.html#healpy.sphtfunc.smoothing
print "Smoothing HEALPix maps to" FWHM

for i in range(0,nmaps):
  map_in = maps[i]
  map_out = healpy.sphtfunc.smoothing(map_in, fwhm)

###### 2.0) Now that the smoothing is finished, we'll write the smoothed maps to a new file:

hp.write_map(images[i]"_smoothed.fits", map_out)

print "Smoothed" images[i] "to FWHM" fwhm "radians"

