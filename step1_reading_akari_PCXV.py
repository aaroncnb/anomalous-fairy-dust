
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
from photutils import SkyCircularAperture as sca
from photutils import aperture_photometry as ap

import healpy as hp
import numpy as np





########0.2) Read in the HEALPix Maps Before starting the LOOP OF ALL REGIONS
print "Reading HEALPix Maps"
map12 = hp.read_map("ame256.fits", nest = nest)
map25 = hp.read_map("../Data/im_iras25.fits", nest = nest)
map60 = hp.read_map("../Data/im_iras60.fits", nest = nest)
map100 = hp.read_map("../Data/healpix10/im_akari65.fits", nest = nest)
map90 = hp.read_map("../Data/healpix10/im_akari90.fits", nest = nest)
map140 = hp.read_map("../Data/healpix10/im_akari140.fits", nest = nest)
map160 = hp.read_map("../Data/healpix10/im_akari160.fits", nest = nest)
map857 = hp.read_map("../Data/im_planck857.fits", nest = nest)
map545 = hp.read_map("../Data/im_planck545.fits", nest = nest)
print "Finished reading HEALPix Maps"
       

## Read in the file containing all of the data about the Planck XV 2013 AME regions.
AME = np.genfromtxt('../Data/AME.txt', delimiter =',')


## 1.1 Circular Aperture Photemotery on the HEALPix Maps (with an annulus for background subtraction):

glon  = AME[:,2]
glat  = AME[:,3]
NROIs = len(AME[:,1])
Nband = len(bands)

positions = SkyCoord(l=glon * units.deg, b=glat * units.deg,
                     frame='galactic')
apertures = SkyCircularAperture(positions, r=apSize * units.arcsec)

annulus_apertures = SkyCircularAnnulus(positions, r_in=bgSizeInner, r_out=bgSizeOuter)


##Here's where the photometry actually happenss, so we'll start the for loop over all the wavebands
for i in range(0, Nbands):
       data = hmaps[i]
       rawflux_table = ap(data, apertures, method ="exact")
       bkgflux_table = ap(data, annulus_apertures, method ="exact")
       phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
       bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()
       bkg_sum = bkg_mean * apertures.area()
       final_sum = phot_table['aperture_sum_raw'] - bkg_sum
       phot_table['residual_aperture_sum'] = final_sum
       print(phot_table['residual_aperture_sum'])

#### IDL example code for getting the HEALPix pixel numbers in side of a circular aperture. 
####I think I need to do it this way rather than using the built-in circular apertere method above.
## get pixels in aperture
    phi = lon*!pi/180.
    theta = !pi/2.-lat*!pi/180.
    ang2vec, theta, phi, vec0
    IF (ordering EQ 'NESTED') THEN BEGIN
        query_disc, nside, vec0, aper_inner_radius/60., innerpix, ninnerpix, /deg, /nest
        query_disc, nside, vec0, aper_outer_radius1/60., outerpix1, nouterpix1, /deg, /nest 
        query_disc, nside, vec0, aper_outer_radius2/60., outerpix2, nouterpix2, /deg, /nest 
    ENDIF ELSE BEGIN
        query_disc, nside, vec0, aper_inner_radius/60., innerpix, ninnerpix, /deg
        query_disc, nside, vec0, aper_outer_radius1/60., outerpix1, nouterpix1, /deg 
        query_disc, nside, vec0, aper_outer_radius2/60., outerpix2, nouterpix2, /deg 
    ENDELSE

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


