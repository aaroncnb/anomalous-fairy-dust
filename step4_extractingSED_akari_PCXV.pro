;PRO step4_extractingSED_akari_, ROI
;;*****************************************************************************
;;*
;;*          International Young Astronomer School on Exploiting
;;*                   Herschel and Planck Data (2013)
;;*
;;*              STEP4: Extracting SEDs from a sample of images
;;*                       (Created by F. Galliano)
;;*
;;*****************************************************************************

;ROIs = [ $
;"G004.24+18.09", $
;"G040.52+02.53", $
;"G107.20+05.20", $
;"G123.13-06.27", $
;"G160.26-18.62", $
;"G173.56-01.76", $
;"G209.01-19.38", $
;"G289.80-01.15", $
;"G305.27+00.15", $
;"G317.51-00.11", $
;"G353.05+16.90"]
;NROIs = N_ELEMENTS(ROIs)
;FOR z=0, NROIs-1 DO BEGIN
;ROI = ROIs[z]



AME = READ_CSV('../Data/AME.txt', COUNT = amy, HEADER = amyhead, MISSING_VALUE='')
print, "AME Data Read"
NROIs = N_ELEMENTS(AME.field01)
FOR cerberus=0,NROIs-1 DO BEGIN

;I think using "i" as the iterator here caused some sort of Hell to open up inside the sub routines. 
; Using "cerberus" instead.

ROI  = AME.field01[cerberus]
glon = AME.field02[cerberus]
glat = AME.field03[cerberus]

;; 1) Rebinning the maps
;;----------------------

;;   a .Load the regridded images.
RESTORE, "../Save/step1_general_.xdr"
RESTORE, "../Save/step3_"+ROI+".xdr"

;;   b. Little trick to use HREBIN (make the dimensions multiple of 5). 
;;   Otherwise, you can always use HCONGRID.
Nx_new = Nx_ref
Ny_new = Ny_ref-1
img_new = img_ref[0:Nx_new-1,0:Ny_new-1,*]
im0 = REFORM(img_ref[*,*,0])
HEXTRACT, im0, hdr_ref, 0, Nx_new-1, 0, Ny_new-1 ;; update the header

;;   c. Actual rebinning.
rescl = 1
Nx_SED = Nx_ref/rescl
Ny_SED = Ny_ref/rescl
img_SED = DBLARR(Nx_SED,Ny_SED,Nband)
mask_SED = BYTARR(Nx_SED,Ny_SED) 
FOR i=0,Nband-1 DO BEGIN
;  whnan = 0.
;  tempimg=img_new[*,*,i]
;  whnan = WHERE(FINITE(tempimg) EQ 0, Nnan)
;  tempimg[whnan] = 0.
;  img_new[*,*,i] = tempimg 
  HREBIN, img_new[*,*,i], hdr_ref, imtemp, hdr_SED, OUTSIZE=[Nx_SED,Ny_SED]
  img_SED[*,*,i] = imtemp
  whbad = WHERE((FINITE(img_SED[*,*,i])) EQ 0, Nbad)
  IF (Nbad GT 0) THEN mask_SED[whbad] = 1

;;Make a stacked FITS file for later inspection
  ;mwrfits, img_SED[*,*,i], "../Figures/step4_SED_"+ROI+"_AkPlIr.fits", hdr
ENDFOR

print, "Rescaled"+ROI+"..."

;;   d. Conversion [MJy/sr] --> [Lsun/Hz/pixel].
;;   This part makes use of the physical constants (!MKS and !ADIM)
;;   defined in idl_init.pro.
;;   This distance is way wrong for stuff inside this galaxy. Gotta fix that, or perhaps remove this step entirely?
;; dist = 50.D ;; [kpc]
dist = 0.14D
;;      Manually setting the pixel FOV. Seems that maybe the PFOV doesn't handle the AKARI header format?
;pixfov = PFOV(hdr_SED) ;; [arcsec/pixel]
pixfov = FXPAR(hdr_sed, 'CDELT2')*3600.
pixsr = pixfov^2*!ADIM.sec2 ;; [sr/pixel]

;img_SED = img_SED * 1.D6*!MKS.Jy * pixsr ;; [W/m2/Hz/pixel]
;Lnu_SED = img_SED * 4*!DPI*(dist*!MKS.kpc)^2 / !MKS.Lsun ;; [Lsun/Hz/pixel]
Lnu_SED = img_SED
;rms_SED = Lnu_SED * 0.1D

;;   e. Display the rebinned images.
;;Edited to display only the FIS data
;;SEVERALIMAGES, { $
;;	         im1:Lnu_SED[*,*,0],  tit1:bands[0], $
;;                 im2:Lnu_SED[*,*,1],  tit2:bands[1], $
;;                 im3:Lnu_SED[*,*,2],  tit3:bands[2], $
;;                 im4:Lnu_SED[*,*,3],  tit4:bands[3] $
;;                 im5:Lnu_SED[*,*,4],  tit5:bands[4]  $
;;                 }, /STERN, /SAMPLE
SEVERALIMAGES, { $
                 im1:Lnu_SED[*,*,0],  tit1:bands[0], $
                 im2:Lnu_SED[*,*,1],  tit2:bands[1], $
                 im3:Lnu_SED[*,*,2],  tit3:bands[2], $
                 im4:Lnu_SED[*,*,3],  tit4:bands[3], $
                 im5:Lnu_SED[*,*,4],  tit5:bands[4], $
                 im6:Lnu_SED[*,*,5],  tit6:bands[5], $
                 im7:Lnu_SED[*,*,6],  tit7:bands[6], $
                 im8:Lnu_SED[*,*,7],  tit8:bands[7], $
                 im9:Lnu_SED[*,*,8],  tit9:bands[8], $
                 im10:Lnu_SED[*,*,9],  tit10:bands[9], $
                 im11:Lnu_SED[*,*,10],  tit11:bands[10], $
                 im12:Lnu_SED[*,*,11],  tit12:bands[11] $ 
         }, /STERN, /SAMPLE


;;   c. Conversion [MJy/sr] --> [Lsun/Hz/pixel].
;;   This part makes use of the physical constants (!MKS and !ADIM)
;;   defined in idl_init.pro.
;rms_SED = rms_SED * 1.D6*!MKS.Jy * pixsr ;; [W/m2/Hz/pixel]
;dLnu_SED_rms = rms_SED * 4*!DPI*(dist*!MKS.kpc)^2 / !MKS.Lsun ;; [Lsun/Hz/pixel]

;;   d. Calibration errors (from the instrument user's manuals).
;;Temporarily commenting out calibration error, since I don't know it yet for the FIS data
;;  Calibration errors are proportional to the flux. Moreover the
;;  error is the same for all pixels (correlation coefficient of 1).
;;Manually setting image-wide errors for each band...
;;      iris_abs_error=[5.1,15.1,10.4,13.5]
;;      AKARI abs error = Estimating 10% for all bands except N60 and N160 (15%)
;;      HFI = 7,7,2,2,2,2

dLnu_SED = DBLARR(Nx_SED,Ny_SED,Nband)
FOR i=0,Nband-1 DO BEGIN
  CASE bands[i] OF
      "akari9": factor_cal = 0.10D
      "akari18": factor_cal = 0.10D
      "iras12": factor_cal = 0.051D
      "iras25": factor_cal = 0.151D
      "iras60": factor_cal = 0.104D
      "iras100": factor_cal = 0.135D
;;    "irac1": factor_cal = 0.03D
;;    "irac2": factor_cal = 0.03D
;;    "irac3": factor_cal = 0.03D
;;    "irac4": factor_cal = 0.03D
;;    "mips24": factor_cal = 0.05D
;;    "mips70": factor_cal = 0.10D
;;    "mips160": factor_cal = 0.12D
;;    "spire250": factor_cal = 0.07D
;;    "spire350": factor_cal = 0.07D
;;    "spire500": factor_cal = 0.07D
;;    "akari9": factor_cal = 0.10D 
;;    "akari18": factor_cal = 0.10D
    "akari65": factor_cal = 0.10D
    "akari90": factor_cal = 0.10D
    "akari140": factor_cal = 0.15D
    "akari160": factor_cal = 0.15D
    "planck857": factor_cal = 0.07D
    "planck545": factor_cal = 0.07D
;;    "planck353": factor_cal = 0.10D
  ENDCASE
  dLnu_SED[*,*,i] = factor_cal * Lnu_SED[*,*,i]
ENDFOR

;;   e. Total error (quadratic sum of the two components).
;;Temporarily setting the total error to exclude calibration error
;dLnu_SED = SQRT(dLnu_SED_rms^2+dLnu_SED_cal^2)
;; dLnu_SED = dLnu_SED_rms


;; 3) Savings and plottings
;;-------------------------
;; Display select SEDs.
;;Changing to FIS bands wavelengths
;;wave = [3.6D,4.5D,5.8D,8.0D,24.D,70.D,160.D,250.D,350.D,500.D] ;; [microns]
;wave = [65.D,90.D,140.D,160.D,345.D,550.D] ;; [microns]
wave = [9.0D,12.0D,18.0D,25.0D,60.0D,65.0D,90.0D,100.0D,140.0D,160.0D,349.81617D,550.07791D] ;; [microns]
;wave = [9.0,12.0,18.0,25.0,60.0,65.0,90.0,100.0,140.0,160.0,345.0,550.0]
;wave = [65.000,90.000,140.000,160.000,349.81617,550.07791] ;; [microns]
;filters= ['AKARI3','AKARI4','AKARI5','AKARI6','HFI1','HFI2']
;filters= ['IRAS1','IRAS2','IRAS3','AKARI3','AKARI4','IRAS4','AKARI5','AKARI6','HFI1','HFI2']
filters= ['AKARI1','IRAS1','AKARI2','IRAS2','IRAS3','AKARI3','AKARI4','IRAS4','AKARI5','AKARI6','HFI1','HFI2']
;instruments = STRTRIM(filters)
instruments = ['AKARI','IRAS','AKARI','IRAS','IRAS','AKARI','AKARI','IRAS','AKARI','AKARI','HFI','HFI']
error = [10.0,5.1,10.0,15.1,10.4,15.0,10.0,13.5,10.0,15.0,7.0,7.0]
error = error/100


nu = !MKS.clight/!MKS.micron/wave ;; [Hz]
@get_colors
WINDOW, 1, XSIZE=800, YSIZE=550
PLOT_OO, [0], [0], /NODAT, XRANGE=[1,1000], YRANGE=[1,100], $
         /XSTYLE, /YSTYLE, XTITLE="!6Wavelength [!7l!6m]", $
         YTITLE="!7m!6L!D!7m!6!N [!D!9n!6!N/pixel]"
LOADCT, 39
whplot = WHERE(mask_SED EQ 0,Nplot)
FOR i=0,Nplot-1 DO BEGIN
  col = i*255./(Nplot-1.)
  Lnu_plot = DBLARR(Nband)
  dLnu_plot = DBLARR(Nband)
  FOR j=0,Nband-1 DO BEGIN
    Lnu_plot[j] = (Lnu_SED[*,*,j])[whplot[i]]
    dLnu_plot[j] = (dLnu_SED[*,*,j])[whplot[i]]
  ENDFOR
;  OPLOT, wave, nu*Lnu_plot, COLOR=col, PSYM=-4
  OPLOT, wave, Lnu_plot, COLOR=col, PSYM=-4
 ; OPLOTERR, wave, nu*Lnu_plot, nu*dLnu_plot
  OPLOTERR, wave, Lnu_plot, dLnu_plot
ENDFOR
@get_colors


;; Now make a quick "average SED" and feed it into The DustEM
;; But first, gotta put the SED and errors and instrument info into a structure
;; that makes sense to DustEM, like the ouput of the "READ_XCAT" routine
;;Averaging part. I don't know if this is the best way.. Perhaps turning the 2d array into a 1d string is best, first of all..
;Lnu_SED_avg = MEAN(MEAN(Lnu_SED,DIMENSION=1,/NAN),DIMENSION=1,/NAN)


    ;;Now that I have a 1D average SED, I just need to put it into the DustEM-preferred style...
    ;;... as in DUSTEM_SET_DATA

;; Savings
;;--------
filexdr = "../Save/step4_"+ROI+".xdr"

SAVE, FILE=filexdr, Nx_SED, Ny_SED, Lnu_SED, dLnu_SED_rms, dLnu_SED_cal, dLnu_SED, $
                     hdr_SED, mask_SED, wave, nu, bands, Nband
PRINT, " - File "+filexdr+" has been written."



ENDFOR
END


