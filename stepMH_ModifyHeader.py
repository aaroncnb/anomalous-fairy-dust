#Import the FITS in/out package for reading and writing

from astropy.io import fits

bands = ['iras12','iras25','iras60','iras100','akari65','akari90','akari140','akari160']
filenames =['im_iras12.fits','im_iras25.fits','im_iras60.fits','im_iras100.fits','im_akari65.fits','im_akari90.fits','im_akari140.fits','im_akari160.fits']
frequencies= ['24983','11992','4997','2998','4612','3331','2141','1874']
#Define an object to hold the HDU data (it'll hold the image data as well as the FITS header in a single object)
#...give the name of the .FITS file you want to read. memmap is for accessing really big files without needing to read the whole of them into RAM

for i in range(0,len(bands)):
    hdulist = fits.open(filenames[i],memmap=True)
    #Now define an object to hold on to the header, specifically (be careful about which header extension you're acessing! Some fits files have two or three.)
    prihdr = hdulist[1].header
    #Here's where we can update keywords in the header. This one adds a keyword "FREQ" and set the value to 2998.
    #...This example is for HEALPix data processing. We want to apply a Planck-data routine to non-Planck data.
    #...THe Planck files usually have a FREQ keyword that tells the photometry routine what their frequency is,
    #...But the AKARI and IRAS data don't have this keyword...so we have to give it one. 2998 GHz is 100 microns.
    prihdr.set('FREQ', frequencies[i])
    #Finally we can write the amended FITS file (well, the header was amended. We didn't actually touch the data!)
    hdulist.writeto('freq_'+filenames[i])
    #Now we "close" the HDU, so python knows we're done with it, and it can be cleared from memory.
    hdulist.close()
#####Tested 7/23/2015######


