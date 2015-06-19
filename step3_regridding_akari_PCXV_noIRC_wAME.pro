;PRO step3_regridding_akari_noirc_wAME.pro
;;*****************************************************************************
;;*
;;*          International Young Astronomer School on Exploiting
;;*                   Herschel and Planck Data (2013)
;;*
;;*                 STEP3: Reprojection on a common grid
;;*                       (Created by F. Galliano)
;;*
;;*****************************************************************************


AME = READ_CSV('../Data/AME.txt', COUNT = amy, HEADER = amyhead, MISSING_VALUE='')
NROIs = N_ELEMENTS(AME.field01)
FOR cerberus=0,NROIs-1 DO BEGIN

ROI  = AME.field01[cerberus]

;; 1) Reading the data
;;--------------------
;;   a. Load back the general information (variables NBAND and BANDS).
RESTORE, "../Save/noirc_wAME/step1_general_.xdr"

;;   b. We choose a reference frame (arbitrary) to reproject all the maps.
;; Using the AME Map image as the reference header (for reasons unknown).

bands =["iras12","iras25","iras60","akari65","akari90","iras100","akari140","akari160","planck857","planck545","ame1"]

Nband = N_ELEMENTS(bands)

ref = "ame1"

;;;We'll try to isolate the problem by using the Step1 (non-convolved) header as the reference

RESTORE, "../Save/noirc_wAME/step1_"+ref+"_"+ROI+".xdr"
;RESTORE,"../Save/noirc_wAME/step1_general_.xdr"

hdr_ref = hdr
Nx_ref = Nx
Ny_ref = Ny

;stop
;; 2) Regridding of the images
;;----------------------------
;;   a. IMG_REF will be a cube containing all the regridded images, since they
;; now have the same size. MASK_REF will be the mask of the common field
;; (keeping only the pixels seen at all wavelengths).

img_ref = DBLARR(Nx_ref,Ny_ref,Nband)
mask_ref = BYTARR(Nx_ref,Ny_ref)

FOR i=0, Nband-1 DO BEGIN
  PRINT, "Regridding "+bands[i]+"..."+ROI+"..."

  ;; We use the IDLastro procedure HASTROM to reproject.

  RESTORE, "../Save/noirc_wAME/step1_"+bands[i]+"_"+ROI+".xdr"

 HASTROM, img, hdr, img_new, hdr_new, hdr_ref, $
           MISSING=!VALUES.D_NaN

  img_ref[*,*,i] = img_new
  
  ;; Update the mask.
  ;; In the end, the unmasked pixels will be those where ALL wavebands
  ;; are defined. If at least one of the waveband is masked, we will
  ;; mask the entire SED of the pixel.
;  whbad = WHERE((FINITE(img_ref[*,*,i]) EQ 0) OR (img_ref[*,*,i] LE 0),Nbad)
whbad = WHERE(img_ref[*,*,i] LE 0,Nbad)
  IF (Nbad GT 0) THEN mask_ref[whbad] = 1

  PRINT
ENDFOR

;;   b. Put NaNs where the mask is active.
whoff = WHERE(mask_ref EQ 1,Noff)
IF (Noff GT 0) THEN BEGIN
  FOR i=0,Nband-1 DO BEGIN
    imtemp = img_ref[*,*,i]
    imtemp[whoff] = !VALUES.D_NaN
    img_ref[*,*,i] = imtemp
  ENDFOR
ENDIF

;;   c. Display (with an home made routine).
;;Adjusting this part for the appropriate number of bands
;SEVERALIMAGES, { $
;                 im1:img_ref[*,*,0]*(1-mask_ref),  tit1:bands[0], $
;                 im2:img_ref[*,*,1]*(1-mask_ref),  tit2:bands[1], $
;                 im3:img_ref[*,*,2]*(1-mask_ref),  tit3:bands[2], $
;                 im4:img_ref[*,*,3]*(1-mask_ref),  tit4:bands[3], $
;                 im5:img_ref[*,*,4]*(1-mask_ref),  tit5:bands[4], $
;                 im6:img_ref[*,*,5]*(1-mask_ref),  tit6:bands[5], $
;                 im7:img_ref[*,*,6]*(1-mask_ref),  tit7:bands[6], $
;		 im8:img_ref[*,*,7]*(1-mask_ref),  tit8:bands[7], $
;                 im9:img_ref[*,*,8]*(1-mask_ref),  tit9:bands[8], $
;                 im10:img_ref[*,*,9]*(1-mask_ref),  tit10:bands[9], $
;                 im11:img_ref[*,*,10]*(1-mask_ref),  tit11:bands[10], $
;                 im12:img_ref[*,*,11]*(1-mask_ref),  tit12:bands[11] $
;          }, /STERN, /SAMPLE

;;   d. Savings.
filexdr = "../Save/noirc_wAME/step3_"+ROI+".xdr"
SAVE, FILE=filexdr, Nx_ref, Ny_ref, hdr_ref, img_ref, mask_ref
PRINT, " - File "+filexdr+" have been written."


ENDFOR
END
