
##*****************************************************************************
##*
##*          International Young Astronomer School on Exploiting
##*                   Herschel and Planck Data (2013)
##*
##*               STEP1: Reading and visualizing the data
##*                       (Created by F. Galliano, edited in an exceedingly crude way by A. Bell)
##*
##*****************************************************************************

     ## 0) Specify the region you want to study (Need this step since we have access to all-sky data)
     ## Manually specifiying the ROI name and the Galactic coordinates. 
     ## I would like to modify this step so that the coordinates are extracted from the ROI name, 
     ## and I want the ROI to be given as a parameter when the "step1..." procedure is called.
     ## I would also like to implement a step7, which will procedurally generate a thesis
     
from photutils import SkyCircularAperture as sca
from photutils import aperture_photometry as ap
from astropy import units as u
from astropy.coordinates import SkyCoord

###Let's import the healpix and numpy packages needed for this analysis...

import healpy as hp
import numpy as np

########0.2) Read in the HEALPix Maps Before starting the LOOP OF ALL REGION
print "Reading HEALPix Maps"
       map12 = hp.read_map("../Data/im_iras12.fits")
       map25 = hp.read_map("../Data/im_iras25.fits")
       map60 = hp.read_map("../Data/im_iras60.fits")
       map100 = hp.read_map("../Data/healpix10/im_akari65.fits")
       map90= hp.read_map("../Data/healpix10/im_akari90.fits")
       map140= hp.read_map("../Data/healpix10/im_akari140.fits")
       map160= hp.read_map("../Data/healpix10/im_akari160.fits")
       map857= hp.read_map("../Data/im_planck857.fits")
       map545= hp.read_map("../Data/im_planck545.fits")
       print "Finished reading HEALPix Maps"
       

## Read in the file containing all of the data about the Planck XV 2013 AME regions.
  AME = np.genfromtxt('../Data/AME.txt', delimiter =',')


## 1.1 Circular Aperture Photemotery on the HEALPix Maps (with an annulus for background subtraction):
 #### First things first, we need to import the aperture photometry package

glon  = AME[:,2]
glat  = AME[:,3]
NROIs = len(AME[:,1])

positions = SkyCoord(l=glon * u.deg, b=glat * u.deg,
                     frame='galactic')
apertures = SkyCircularAperture(positions, r=apSize. * u.arcsec)

annulus_apertures = SkyCircularAnnulus(positions, r_in=bgSizeInner, r_out=bgSizeOuter)

rawflux_table =ap(data, apertures)
bkgflux_table = ap(data, annulus_apertures)
phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()
bkg_sum = bkg_mean * apertures.area()
final_sum = phot_table['aperture_sum_raw'] - bkg_sum
phot_table['residual_aperture_sum'] = final_sum
print(phot_table['residual_aperture_sum'])


data = map
phot_table = ap(data, apertures, method ="exact")
print(phot_table)


##IDL Example Code for estimating error based on the background annulus, from Clive Dickinson's (2010) IDL Code
; estimate error based on robust sigma of the background annulus
        IF (keyword_set(noise_model)) EQ 0 THEN $
            noise_model = 0 $
        ELSE noise_model = 1

; new version (2-Dec-2010) that has been tested with simulations for
; Planck early paper and seems to be reasonable for many applications

        IF (noise_model EQ 0) THEN BEGIN 
            Npoints = (pix_area*ninnerpix) / $
                      (1.13*(float(res_arcmin)/60. *!PI/180.)^2)
            Npoints_outer = (pix_area*nouterpix) / $
                      (1.13*(float(res_arcmin)/60. *!PI/180.)^2)
            fd_err = stddev(map[outerpix,column]) * factor * ninnerpix / sqrt(Npoints)
        ENDIF

; (works exactly for white uncorrelated noise only!)
        IF (noise_model EQ 1) THEN BEGIN
            k = !PI/2.

            fd_err = factor * sqrt(float(ninnerpix) + $
                (k * float(ninnerpix)^2/nouterpix)) * robust_sigma(map[outerpix,column])
        ENDIF


for cerberus in range(0,NROIs-1) 

#ROI = AME[cerberus,1]

##There seems to be a problem with converting the ROI name into GLON and GLAT. It's writes them into the header as a string...
##Somehow GNOMDRIZZ can interpret the string, but later on HASTROM (in step3) cannot
#glon = STRMID(ROI, 1,6)
#glat = STRMID(ROI,7)


Nband = len(bands)

#for i in range (0, Nband-1)

  ## 1) Extract rectangular cutouts from the HEALPix arrays (loaded above) using gnomdrizz
   ## As far as smoothing goes, we can apply a smoothing to the whole map from the start using the following function:
    ##    wmap_map_I_smoothed = hp.smoothing(wmap_map_I, fwhm=60, arcmin=True)
   ## It will take a long time but we should only have to do it once, assuming the smoothing criteria doesn't change a whole lot

    CASE bands[i] OF
       'akari9':BEGIN
       img = readfits("../Data/akari9/"+ROI+"_"+bands[i]+".fits",hdr)
       img[where(img[*] eq 0.00000)] = !Values.F_NAN       
       #img = img*0.31       
       #backsub, img
       #img = FILTER_IMAGE(img, FWHM=  120, /ALL_PIXELS)
    END
       'akari18':BEGIN
       img = readfits("../Data/akari18/"+ROI+"_"+bands[i]+".fits",hdr)
       img[where(img[*] eq 0.00000)] = !Values.F_NAN
       img = img*0.48
       #backsub, img
       #img = FILTER_IMAGE(img, FWHM=  120, /ALL_PIXELS)
    END
      'iras12':BEGIN
       gnomdrizz,hparrays.map12[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'iras25':BEGIN
       gnomdrizz,hparrays.map25[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'iras60':BEGIN
       gnomdrizz,hparrays.map60[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
      'akari65':BEGIN
       gnomdrizz,hparrays.map65[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'akari90':BEGIN
        gnomdrizz,hparrays.map90[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'iras100':BEGIN
       gnomdrizz,hparrays.map100[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'akari140':BEGIN
       gnomdrizz,hparrays.map140[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'akari160':BEGIN
       gnomdrizz,hparrays.map160[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'spire250':BEGIN
       img = readfits("../Data/"+ROI+"/im_"+bands[i]+".fits",hdr)
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'planck857':BEGIN
       gnomdrizz,hparrays.map857[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END
       'planck545':BEGIN
       gnomdrizz,hparrays.map545[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    img = FILTER_IMAGE(img, FWHM= 35, /ALL_PIXELS)
    END

ENDCASE



  ##  - Extension 1: RMS in MJy/sr.
  ## In the case of IRC, just make an array with the same dimension as the image. Then set all values to 10% of the signal
  
rms = img*0.1D
#rms_total = 

 
 
 
  ## 1.2) Background Subtraction: 

#backsub,img
#skyfit,img,skyimg,YORDER=0,XORDER=0, NDEG=1
#img = img-skyim

#writefits, ROI+bands[i]+'.fits',img, hdr



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
