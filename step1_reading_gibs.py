
##*****************************************************************************
##*  An all-sky, multi-wavelength image analysis code. 
##*Designed for checking the full dust SED accross the mid infrared to far infrared wavelengths 
##* for comparison to the anomalous microwave emission.
##* Based on the example photometry IDL code Created by F. Galliano for the Herschel-Planck School 2013
##*   -Aaron C. Bell (2015)
##*****************************************************************************
###Let's import the healpix and numpy packages needed for this analysis...

from astropy import units
from astropy.coordinates import SkyCoord

#### First things first, we need to import the aperture photometry packages
from photutils import SkyCircularAperture
from photutils import aperture_photometry

import healpy as hp
import numpy as np

########0.2) Read in the FITS 

print "Reading FITS Maps"

hdulist1 = fits.open('../Data/gibson/im_irac1.fits')
hdulist2 = fits.open('../Data/gibson/im_irac2.fits')
hdulist3 = fits.open('../Data/gibson/im_irac3.fits')
hdulist4 = fits.open('../Data/gibson/im_irac4.fits')
hdulist5 = fits.open('../Data/gibson/im_mips24.fits')
hdulist6 = fits.open('../Data/gibson/im_mips70.fits')
hdulist7 = fits.open('../Data/gibson/im_mips160.fits')

#########0.3)  Save the FITS HDUs as nummpy arrays (like you would in IDL)

image1 = hdulist1[0].data
image2 = hdulist2[0].data
image3 = hdulist3[0].data
image4 = hdulist4[0].data
image5 = hdulist5[0].data
image6 = hdulist6[0].data
image7 = hdulist7[0].data

print "Finished reading FITS Maps"

##########0.4) Print the dimensions of the arrays 

print(type(image_data))
print(image_data.shape)
       

## Read in the file containing all of the data about the Planck XV 2013 AME regions.
AME = np.genfromtxt('../Data/gibs_apertures.txt', delimiter =',')


## 1.1 Circular Aperture Photemotery on the HEALPix Maps (with an annulus for background subtraction):

glon  = gibs[:,1]
glat  = gibs[:,2]
NROIs = len(gibs[:,1])
Nband = len(bands)

positions = SkyCoord(l=glon * units.deg, b=glat * units.deg,
                     frame='galactic')
apertures = SkyCircularAperture(positions, r=apSize * units.arcmin)

annulus_apertures = SkyCircularAnnulus(positions, r_in=bgSizeInner * units.arcmin, r_out=bgSizeOuter * units.arcmin)



#### IDL example code for getting the HEALPix pixel numbers in side of a circular aperture. 
####I think I need to do it this way rather than using the built-in circular apertere method above.
## get pixels in aperture
#First step in doing this is to get our coordinates into a format that "ang2vec" can understand.

#phi = glon * np.pi / 180.


#theta = (np.pi / 2.)- (glat * np.pi / 180.)

#ang2vec converts a spherical angle to a position vector, but it needs phi and theta as co-lat and co-lon

#vec0 = hp.ang2vec(theta, phi)    


#query_disc was recently translated into python from the IDL Healpix package. It isn't mentioned in the Healpy documentation, but it does exist.
#It takes phi and theta, as colatitude
#healpy.query_disc(num_sides, vec, radius_angle, nested_scheme, deg=False)

#Start the list making...

#listpix_r1 = []
#listpix_r2 = []
#listpix_r3 = []

#for i in range(0,98):
#    listpix_r1.append(hp.query_disc(nside, vec0[i], aper_inner_radius/60., True, True))
#    listpix_r2.append(hp.query_disc(nside, vec0[i], aper_outer_radius1/60., True, True))
#    listpix_r3.append(hp.query_disc(nside, vec0[i], aper_outer_radius2/60., True, True))
    
    
#find pixels in the annulus (between outerradius1 and outeradius2) 
#    outerpix_all = [listpix_r2,listpix_r3]
#    outerpix_all.sort()
#    outerpix = -1L
#    for i in range(0, nouterpix2-1):
#        temp = np.where(listpix_r2 = listpix_r3[i])
#        if (temp[0] < 0):
#               outerpix = [outerpix,outerpix2[i]]
#    outerpix = outerpix[1:*]
#    nouterpix = len(outerpix)

##Here's where the photometry actually happenss, so we'll start the for loop over all the wavebands
#for i in range(0, Nbands):
#       data = hmaps[i]
#       rawflux_table = ap(data, apertures, method ="exact")
#       bkgflux_table = ap(data, annulus_apertures, method ="exact")
#       phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
#       bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()
#       bkg_sum = bkg_mean * apertures.area()
#       final_sum = phot_table['aperture_sum_raw'] - bkg_sum
#       phot_table['residual_aperture_sum'] = final_sum
#       print(phot_table['residual_aperture_sum'])




##IDL Example Code for estimating error based on the background annulus, from Clive Dickinson's (2010) IDL Code
# estimate error based on robust sigma of the background annulus
#if noise_model== 0:
#    Npoints = (pix_area*ninnerpix)/(1.13*np.float(res_arcmin))/60*np.PI/180.)^2)
#    Npoints_outer = (pix_area*nouterpix)/(1.13*np.stddev(res_arcmin)/60*np.PI/180)^2)
#    fd_err = stddev(map[outerpix,column])*factor*ninnerpix/np.sqrt(Npoints)
#else:
#    k = np.PI/2.
#    fd_err = factor * np.sqrt(np.float(ninnerpix) + (k * np.float(ninnerpix)^2/nouterpix)) * robust_sigma(map[outerpix,column]#
#
## 4) Plotting a rectangular cut-out around the aperture for display purposes. Let's make it about 3-5 degrees wide:


