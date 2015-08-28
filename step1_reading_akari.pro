PRO step1_reading_akari, ROI, glon, glat
;;*****************************************************************************
;;*
;;*          International Young Astronomer School on Exploiting
;;*                   Herschel and Planck Data (2013)
;;*
;;*               STEP1: Reading and visualizing the data
;;*                       (Created by F. Galliano, edited in an exceedingly crude way by A. Bell)
;;*
;;*****************************************************************************

     ;; 0) Specify the region you want to study (Need this step since we have access to all-sky data)
     ;; Manually specifiying the ROI name and the Galactic coordinates. 
     ;; I would like to modify this step so that the coordinates are extracted from the ROI name, 
     ;; and I want the ROI to be given as a parameter when the "step1..." procedure is called.
     ;; I would also like to implement a step7, which will procedurally generate a thesis


;; List of the broadbands to be used (Just throwing in the AKARI bands here)
;;----------------------------------
;; This array contains the string label of each of the waveband we
;; want to  analyze.
;; Currently testing with only the FIS and HFI bands- all of the ones available as HEAPIX.
;; Adding IRC bands will need IRC HEALPIX files, or some parallel process which gets the data via the "ircscan" bash script by Ohsawa-san


IRC_present = FILE_TEST("../Data/"+ROI+"/im_akari9.fits")
;REMOVE, WHERE(has_band EQ 0), bands
IF IRC_present EQ 1 THEN BEGIN 
bands =["akari9","iras12","akari18","iras25","iras60","akari65","akari90","iras100","akari140","akari160","planck857","planck545"]
ENDIF ELSE BEGIN
bands =["iras12","iras25","iras60","akari65","akari90","iras100","akari140","akari160","planck857","planck545"]
ENDELSE


;;bands = ["akari9", "akari18", "akari65", "akari90","akari140","akari160","planck857","planck545"]
;;bands = ["akari65","akari90","akari140","akari160","planck857","planck545"]
;;bands =["iras12","iras25","iras60","akari65","akari90","iras100","akari140","akari160","planck857","planck545"]
;bands =["akari9","iras12","akari18","iras25","iras60","akari65","akari90","iras100","akari140","akari160","planck857","planck545"]


Nband = N_ELEMENTS(bands)

FOR i=0,Nband-1 DO BEGIN


  ;; 1) Read the FITS file 
  ;; Now that HEALPIX FITS files are available for FIS as well as HFI data, I will modify this step to extract cutouts from HPY maps
  ;;----------------------
  PRINT, "Reading the map "+bands[i]


 

;order 13 pixel spacing in degrees is 0.0071625
;0.42975 in arcminutes...
;     reso = 0.42975 ELSE (Order 13pr)
;;;ORDER 13 is too pixely, and IDL chokes on it :( compromising for now with ORDER 12
;order 12 pixel spacing is 0.014325 degrees
;
reso1 = 1.7
reso2 = 0.8595
size_arcmin = 30
;; Read-in each band's data using the appropriate routine for each (i.e. gnomdrizz for HEALPix, simple READFITS for non HEALPix)
;; The procedure READ_FITS_MAP is part of the HEALPiX IDL library.
;; GNOMDRIZZ is from the DRIZZLIB package (And is a very silly-sounding routine name)

    CASE bands[i] OF
       'akari9':BEGIN
       img = readfits("../Data/"+ROI+"/im_"+bands[i]+".fits",hdr)
    END
       'akari18':BEGIN
       img = readfits("../Data/"+ROI+"/im_"+bands[i]+".fits",hdr)
    END
      'iras12':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=1.7,pxsize=ROUND(47),mapfits=img,hfits=hdr, /nested
    END
       'iras25':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    END
       'iras60':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    END
      'akari65':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso2,pxsize=ROUND(size_arcmin/reso2),mapfits=img,hfits=hdr, /nested
    END
       'akari90':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso2,pxsize=ROUND(size_arcmin/reso2),mapfits=img,hfits=hdr, /nested
    END
       'iras100':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    END
       'akari140':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso2,pxsize=ROUND(size_arcmin/reso2),mapfits=img,hfits=hdr, /nested
    END
       'akari160':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso2,pxsize=ROUND(size_arcmin/reso2),mapfits=img,hfits=hdr, /nested
    END
       'planck857':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    END
       'planck545':BEGIN
       read_fits_map,"../Data/im_"+bands[i]+".fits", map, head, xhead
       gnomdrizz,map[*,0],rot=[glon,glat],reso_arcmin=reso1,pxsize=ROUND(size_arcmin/reso1),mapfits=img,hfits=hdr, /nested
    END

ENDCASE



  ;;  - Extension 1: RMS in MJy/sr.
  ;; In the case of IRC, just make an array with the same dimension as the image. Then set all values to 10% of the signal
  ;rms = img*0.1D

  ;; 1.1) Background Subtraction: Let's try to do a background subtraction inside of this step.
  ;; I think "SKY" is too crude, so I will instead use "BACKSUB" from the Buie library...
  ;;.... It fits a linear or polynomial background to each column or each row of the image, then subtracts it
  ;; The default method is to subtract a mean.
    ;;backsub, image, [/ROW,/COLUMN]


  ;;Make a stacked FITS file for later inspection
;mwrfits, img_SED[*,*,i], '../Figures/step1_SED_AkPlIr.fits', hdr

;backsub, img, /ROW, order=1
;backsub,img
   ;;This one can fit a plane or 2-d polynomial to the whole image...which is best?
   ;; I think that if the zody is the dominant background, a plane is best.

;skyfit,img,img


  ;; Size of the image
  sizim = SIZE(img)
  Nx = sizim[1]
  Ny = sizim[2]


  ;; 2) Make a mask
  ;;---------------
  ;; Sometimes, pixels are NaN (Not a Number), either because they
  ;; have been masked out during the data processing, or they are simply off
  ;; the field. Here we build a mask for each map. This mask will be 0 if 
  ;; the pixel is OK and 1 otherwise.
  mask = BYTARR(Nx,Ny)
whbad = WHERE(img LE 0,Nbad)
  IF (Nbad GT 0) THEN mask[whbad] = 1


  ;; 3) Statistics on the data
  ;;--------------------------
  ;; Signal-to-noise keeping only the NaNaNs (Not a NaN).
;  SovN = (ABS(img/rms))[WHERE(mask EQ 0)]
;  PRINT, "  The signal-to-noise ratio for the map is on average " $
;       + STRTRIM(MEAN(SovN),2)+" between "+STRTRIM(MIN(SovN),2)+" and " $
;       + STRTRIM(MAX(SovN),2)+"."
  ;; In all these programs, we use the STRTRIM(var,2) formatting
  ;; command, that converts a number to a character string and removes
  ;; the extra spaces before and after.


  ;; 4) Plotting the image
  ;;----------------------
  ;; Plot the image, with a reduce size.
  maxsize = 500.
  reduction = maxsize/MAX([Nx,Ny])
  Nx_plot = ROUND(reduction*Nx)
  Ny_plot = ROUND(reduction*Ny)
  WINDOW, i, XSIZE=Nx_plot, YSIZE=Ny_plot, TITLE=bands[i], $
             XPOS=(i/2 MOD 2)*maxsize, YPOS=((i+1)/2 MOD 2)*maxsize
  LOADCT, 15
  TVSCL, CONGRID(img*(1-mask),Nx_plot,Ny_plot)
 ; ;CONTOUR, DIST(31,41), POSITION=[0.15, 0.15, 0.95, 0.75], $
;  ;C_COLORS=INDGEN(25)*4, NLEVELS=25
  
  ;; The IDL CONGRID function interpolates an array on a new
  ;; grid. It has to be handled carefully when performing photometry
  ;; because of possible interpolation errors. However, for display
  ;; purposes, as in this case, it works fine.


  ;; 5) Savings
  ;;-----------
  ;; We write one IDL save file per waveband.
  filexdr = "../Save/step1_"+bands[i]+"_"+ROI+".xdr"
;  SAVE, FILE=filexdr, Nx, Ny, img, mask, rms, hdr
  SAVE, FILE=filexdr, Nx, Ny, img, mask, hdr
  PRINT, " - File "+filexdr+" has been written."
    ;;5b) We need to clear the 'hdr' variable, since GNOMDRIZZ will use it as input header if it's already defined, i.e. by READFITS  
    ;;  Using "UNDEFINE" as suggested by CoyoteIDL, since "DELVAR" deletes at the main level and causes problems for reasons.
     UNDEFINE,hdr
      

  PRINT
ENDFOR

;; 6) Save a general file
;;-----------------------
filexdr = "../Save/step1_general_"+ROI+".xdr"
SAVE, FILE=filexdr, bands, Nband
PRINT, " - File "+filexdr+" has been written."



END
