
##*****************************************************************************
##*  An all-sky, multi-wavelength image analysis code. 
##*Designed for checking the full dust SED accross the mid infrared to far infrared wavelengths 
##* for comparison to the anomalous microwave emission.
##* Based on the example photometry IDL code Created by F. Galliano for the Herschel-Planck School 2013
##*   -Aaron C. Bell (2015)
##*****************************************************************************

from photutils import SkyCircularAperture as sca
from photutils import aperture_photometry as ap
from astropy import units as u
from astropy.coordinates import SkyCoord

###Let's import the healpix and numpy packages needed for this analysis...

import healpy as hp
import numpy as np

########0.2) Read in the HEALPix Maps Before starting the LOOP OF ALL REGIONS
print "Reading HEALPix Maps"
map12 = hp.read_map("../Data/im_iras12.fits")
map25 = hp.read_map("../Data/im_iras25.fits")
map60 = hp.read_map("../Data/im_iras60.fits")
map100 = hp.read_map("../Data/healpix10/im_akari65.fits")
map90 = hp.read_map("../Data/healpix10/im_akari90.fits")
map140 = hp.read_map("../Data/healpix10/im_akari140.fits")
map160 = hp.read_map("../Data/healpix10/im_akari160.fits")
map857 = hp.read_map("../Data/im_planck857.fits")
map545 = hp.read_map("../Data/im_planck545.fits")
print "Finished reading HEALPix Maps"
       

## Read in the file containing all of the data about the Planck XV 2013 AME regions.
AME = np.genfromtxt('../Data/AME.txt', delimiter =',')


## 1.1 Circular Aperture Photemotery on the HEALPix Maps (with an annulus for background subtraction):
#### First things first, we need to import the aperture photometry package

glon  = AME[:,2]
glat  = AME[:,3]
NROIs = len(AME[:,1])
Nband = len(bands)

positions = SkyCoord(l=glon * u.deg, b=glat * u.deg,
                     frame='galactic')
apertures = SkyCircularAperture(positions, r=apSize. * u.arcsec)

annulus_apertures = SkyCircularAnnulus(positions, r_in=bgSizeInner, r_out=bgSizeOuter)


##Here's where the photometry actually happenss, so we'll start the for loop over all the wavebands
for i in range(0, Nbands)

data = hmaps[i]
rawflux_table = ap(data, apertures, method ="exact")
bkgflux_table = ap(data, annulus_apertures, method ="exact")

phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()
bkg_sum = bkg_mean * apertures.area()
final_sum = phot_table['aperture_sum_raw'] - bkg_sum
phot_table['residual_aperture_sum'] = final_sum
print(phot_table['residual_aperture_sum'])




##IDL Example Code for estimating error based on the background annulus, from Clive Dickinson's (2010) IDL Code
# estimate error based on robust sigma of the background annulus
     if noise_model == 0 :

            Npoints = (pix_area*ninnerpix) / $
                      (1.13*(float(res_arcmin)/60. *!PI/180.)^2)
            Npoints_outer = (pix_area*nouterpix) / $
                      (1.13*(float(res_arcmin)/60. *!PI/180.)^2)
            fd_err = stddev(map[outerpix,column]) * factor * ninnerpix / sqrt(Npoints)
else:
            k = !PI/2.

            fd_err = factor * sqrt(float(ninnerpix) + $
                (k * float(ninnerpix)^2/nouterpix)) * robust_sigma(map[outerpix,column])


  ## Size of the image
  sizim = SIZE(img)
  Nx = sizim[1]
  Ny = sizim[2]


  ## 2) Make a mask
  ##---------------
  ## Sometimes, pixels are NaN (Not a Number), either because they
  ## have been masked out during the data processing, or they are simply off
  ## the field. Here we build a mask for each map. This mask will be 0 if 
  ## the pixel is OK and 1 otherwise.

  mask = BYTARR(Nx,Ny)
  whbad = WHERE((img LE 0.D) AND (FINITE(img) EQ 0), Nbad)
  IF (Nbad GT 0) THEN mask[whbad] = 1


  ## 3) Statistics on the data
  ##--------------------------
  ## Signal-to-noise keeping only the NaNaNs (Not a NaN).
    SovN = (ABS(img/rms))[WHERE(mask EQ 0)]
    PRINT, "  The signal-to-noise ratio for the map is on average " $
         + STRTRIM(MEAN(SovN),2)+" between "+STRTRIM(MIN(SovN),2)+" and " $
          + STRTRIM(MAX(SovN),2)+"."

  ## In all these programs, we use the STRTRIM(var,2) formatting
  ## command, that converts a number to a character string and removes
  ## the extra spaces before and after.


  ## 4) Plotting the image
  ##----------------------
  ## Plot the image, with a reduce size.

  maxsize = 500.
  reduction = maxsize/MAX([Nx,Ny])
  Nx_plot = ROUND(reduction*Nx)
  Ny_plot = ROUND(reduction*Ny)
  WINDOW, i, XSIZE=Nx_plot, YSIZE=Ny_plot, TITLE=bands[i], $
             XPOS=(i/2 MOD 2)*maxsize, YPOS=((i+1)/2 MOD 2)*maxsize
  LOADCT, 15
  TVSCL, CONGRID(img*(1-mask),Nx_plot,Ny_plot)

  #CONTOUR, DIST(31,41), POSITION=[0.15, 0.15, 0.95, 0.75], $
  #C_COLORS=INDGEN(25)*4, NLEVELS=25
  
  ## The IDL CONGRID function interpolates an array on a new
  ## grid. It has to be handled carefully when performing photometry
  ## because of possible interpolation errors. However, for display
  ## purposes, as in this case, it works fine.


  ## 5) Savings
  ##-----------
  ## We write one IDL save file per waveband.
  filexdr = "../Save/step1_"+bands[i]+"_"+ROI+".xdr"
  SAVE, FILE=filexdr, Nx, Ny, img, mask, rms, hdr
#  SAVE, FILE=filexdr, Nx, Ny, img, mask, hdr
  PRINT, " - File "+filexdr+" has been written."
    ##5b) We need to clear the 'hdr' variable, since GNOMDRIZZ will use it as input header if it's already defined, i.e. by READFITS  
    ##  Using "UNDEFINE" as suggested by CoyoteIDL, since "DELVAR" deletes at the main level and causes problems for reasons.
     UNDEFINE,hdr
      

  PRINT
ENDFOR

## 6) Save a general file
##-----------------------
filexdr = "../Save/step1_general_.xdr"
SAVE, FILE=filexdr, bands, Nband
PRINT, " - File "+filexdr+" has been written."


ENDFOR
END
