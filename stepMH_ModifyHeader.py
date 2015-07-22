#Import the FITS in/out package for reading and writing

from astropy.io import fits


#Define an object to hold the HDU data (it'll hold the image data as well as the FITS header in a single object)
#...give the name of the .FITS file you want to read. memmap is for accessing really big files without needing to read the whole of them into RAM

hdulist = fits.open('im_iras100.fits',memmap=True)


#Now define an object to hold on to the header, specifically (be careful about which header extension you're acessing! Some fits files have two or three.)

prihdr = hdulist[1].header


#Here's where we can update keywords in the header. This one adds a keyword "FREQ" and set the value to 2998.
#...This example is for HEALPix data processing. We want to apply a Planck-data routine to non-Planck data.
#...THe Planck files usually have a FREQ keyword that tells the photometry routine what their frequency is,
#...But the AKARI and IRAS data don't have this keyword...so we have to give it one. 2998 GHz is 100 microns.

prihdr.set('FREQ', '2998')

#Finally we can write the amended FITS file (well, the header was amended. We didn't actually touch the data!)

hdulist.writeto('im_iras100_new.fits')

#Now we "close" the HDU, so python knows we're done with it, and it can be cleared from memory.

hdulist.close()



