;PRO step2_degrading_akari_
;;*****************************************************************************
;;*
;;*          International Young Astronomer School on Exploiting
;;*                   Herschel and Planck Data (2013)
;;*
;;*               STEP2: Degrading the spatial resolution
;;*                       (Created by F. Galliano)
;;*
;;*****************************************************************************


AME = READ_CSV('../Data/AME.txt', COUNT = amy, HEADER = amyhead, MISSING_VALUE='')
NROIs = N_ELEMENTS(AME.field01)
FOR cerberus=0,NROIs-1 DO BEGIN

ROI  = AME.field01[cerberus]



;; Load back the general information (variables NBAND and BANDS).
RESTORE, "../Save/step1_general_.xdr"
time0 = SYSTIME(1)
;; Big loop on the bands.
FOR i=0,Nband-1 DO BEGIN

  PRINT, "Degrading "+bands[i]+"..."+ROI+"..."


  ;; 1) Restore the data
  ;;--------------------
  ;; Reading the precomputed convolution kernels.
  kernel = READFITS("../Kernels/"+bands[i]+"_to_planck857.fits",hdr_kernel)

  ;; Reading the image stored during step 1.
  RESTORE, "../Save/step1_"+bands[i]+"_"+ROI+".xdr"


  ;; 2) Degrade the resolution of the image
  ;;---------------------------------------
  ;;   a. We degrade the image to the lower resolution (planck857).
  

;DO_THE_CONVOLUTION, img, hdr, kernel, hdr_kernel, img_planck857, hdr_planck857, $
;                      result_kernel_image, result_kernel_header, 0

  ;;   b. Make a mask for the degraded map.
  mask = BYTARR(Nx,Ny) 
  whbad = WHERE((img LE 0.D) AND (FINITE(img) EQ 0), Nbad)
;  whbad = WHERE(img_planck857 LE 0, Nbad)
  IF (Nbad GT 0) THEN mask[whbad] = 1
;  mask_planck857[0,*] = 1
;  mask_planck857[Nx-1,*] = 1
;  mask_planck857[*,0] = 1
;  mask_planck857[*,Ny-1] = 1


  ;;   c. Plot the image, with a reduced size.
  maxsize = 500.
  reduction = maxsize/MAX([Nx,Ny])
  Nx_plot = ROUND(reduction*Nx)
  Ny_plot = ROUND(reduction*Ny)
  WINDOW, i, XSIZE=Nx_plot, YSIZE=Ny_plot, TITLE=bands[i], $
             XPOS=(i/2 MOD 2)*maxsize, YPOS=((i+1)/2 MOD 2)*maxsize
  LOADCT, 15
  TVSCL, CONGRID(img*(1-mask),Nx_plot,Ny_plot)

  ;;   d. Save the signal image.
  filexdr = "../Save/step2_"+bands[i]+"_"+ROI+".xdr"
  SAVE, FILE=filexdr, Nx, Ny, img, hdr, mask
  PRINT, " - File "+filexdr+" has been written."


  ;; 3) Error propagation
  ;;---------------------
;  do_error = 0 ;; Warning! Setting this variable to 1 results in a long 
;               ;; computation time...
;  Nmc = 5 ;; number of Monte-Carlo iterations for error calculation.
;  IF (do_error EQ 1) THEN BEGIN;
;;
;    ;;   a. Generate the randomly perturbed images assuming normal
;    ;;   independent noise.
;    img_MC = DBLARR(Nx,Ny,Nmc)
;    FOR j=0,Nmc-1 DO img_MC[*,*,j] = img + rms*RANDOMN(seed,Nx,Ny)*(1-mask)

;    ;;   b. We degrade each perturbed images to the lower resolution (planck857).
;    img_planck857_MC = DBLARR(Nx,Ny,Nmc)
;    FOR j=0,Nmc-1 DO BEGIN
;      PRINT, " Error calculation "+STRTRIM(j+1,2)+"/"+STRTRIM(Nmc,2), $
;             " ("+STRTRIM((SYSTIME(1)-time0)/60,2)+" minutes elapsed)"
;      DO_THE_CONVOLUTION, REFORM(img_MC[*,*,j]), hdr, kernel, hdr_kernel, $
;                          img_out, hdr_out, result_kernel_image, $
;                          result_kernel_header, 0
;      img_planck857_MC[*,*,j] = img_out
;    ENDFOR

;img_planck857_MC = DBLARR(Nx,Ny,Nmc)
;    FOR j=0,Nmc-1 DO BEGIN
;      PRINT, " Error calculation "+STRTRIM(j+1,2)+"/"+STRTRIM(Nmc,2), $
;             " ("+STRTRIM((SYSTIME(1)-time0)/60,2)+" minutes elapsed)"
;      DO_THE_CONVOLUTION, REFORM(img_MC[*,*,j]), hdr, kernel, hdr_kernel, $
;                          img_out, hdr_out, result_kernel_image, $
;                          result_kernel_header, 0
;      img_planck857_MC[*,*,j] = img_out
;    ENDFOR

;
;    ;;   c. Take the standard deviation of the cube of Monte-Carlo images.
;    rms_planck857 = STDDEV(img_planck857_MC,DIM=3)
;    PRINT, "   Mean S/N = "+STRTRIM(MEAN(img_planck857/rms_planck857,/NaN),2) $
;           + " between "+STRTRIM(MIN(img_planck857/rms_planck857,/NaN),2)+" and " $
;           + STRTRIM(MAX(img_planck857/rms_planck857,/NaN),2)
;
;    ;;   d. Save the noise images.
;    filexdr = "../Save/step2_"+bands[i]+"_"+ROI+"MC.xdr"
;    SAVE, FILE=filexdr, Nx, Ny, Nmc, img_planck857_MC, hdr_planck857, mask_planck857,$
;                        rms_planck857
;    PRINT, " - File "+filexdr+" has been written."
;
;  ENDIF

;  PRINT
UNDEFINE, PSF
UNDEFINE, TABLE_NUMBER
ENDFOR

ENDFOR


END
