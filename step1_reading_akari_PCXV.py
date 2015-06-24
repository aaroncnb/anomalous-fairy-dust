
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
       map160= ("../Data/healpix10/im_akari160.fits")
       map857= hp.read_map("../Data/im_planck857.fits")
       map545= hp.read_map("../Data/im_planck545.fits")
       print "Finished reading HEALPix Maps"

AME = READ_CSV('../Data/AME.txt', COUNT = amy, HEADER = amyhead, MISSING_VALUE='')
NROIs = N_ELEMENTS(AME.field01)
FOR cerberus=0,NROIs-1 DO BEGIN

ROI  = AME.field01[cerberus]

##There seems to be a problem with converting the ROI name into GLON and GLAT. It's writes them into the header as a string...
##Somehow GNOMDRIZZ can interpret the string, but later on HASTROM (in step3) cannot
#glon = STRMID(ROI, 1,6)
#glat = STRMID(ROI,7)
glon = AME.field02[cerberus]
glat = AME.field03[cerberus]



## List of the broadbands to be used (Just throwing in the AKARI bands here)
##----------------------------------
## This array contains the string label of each of the waveband we
## want to  analyze.
## Currently testing with only the FIS and HFI bands- all of the ones available as HEAPIX.
## Adding IRC bands will need IRC HEALPIX files, 
#######or some parallel process which gets the data via the "ircscan" bash script by Ohsawa-san
#Remov
#bands = ["akari9","iras12","akari18","iras25","iras60","akari65","akari90","iras100","akari140","akari160","planck857","planck545"]



Nband = N_ELEMENTS(bands)

FOR i=0,Nband-1 DO BEGIN


  ## 1) Read the FITS file 
  ## Now that HEALPIX FITS files are available for FIS as well as HFI data, I will modify this step to extract cutouts from HPY maps
  ##----------------------
#  PRINT, "Reading the map "+bands[i]+"..."+ROI+"..."+glon+"..."+glat


#order 13 pixel spacing in degrees is 0.0071625
#0.42975 in arcminutes...
#     reso = 0.42975 ELSE (Order 13pr)
###ORDER 13 is too pixely, and IDL chokes on it :( compromising for now with ORDER 12
#order 12 pixel spacing is 0.014325 degrees
#
reso1 = 1.7
#reso2 = 0.8595
size_arcmin = 240
smoothing_size_arcmin = 60
## Read-in each band's data using the appropriate routine for each (i.e. gnomdrizz for HEALPix, simple READFITS for non HEALPix)
## The procedure READ_FITS_MAP is part of the HEALPiX IDL library.
## GNOMDRIZZ is from the DRIZZLIB package (And is a very silly-sounding routine name)

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

  ## 1.1) Background Subtraction: Let's try to do a background subtraction inside of this step.
  ## I think "SKY" is too crude, so I will instead use "BACKSUB" from the Buie library...
  ##.... It fits a linear or polynomial background to each column or each row of the image, then subtracts it
  ## The default method is to subtract a mean.

#backsub,img



writefits, ROI+bands[i]+'.fits',img, hdr

   ##This one can fit a plane or 2-d polynomial to the whole image...which is best?
   ## I think that if the zody is the dominant background, a plane is best.
   ##    Continnues to have some issue with the "LOWESS" routine
   #skyfit,img,skyimg,YORDER=0,XORDER=0, NDEG=1
   #img = img-skyimg


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
