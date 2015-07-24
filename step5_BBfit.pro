s;PRO step5_BBfit_akari_
;;*****************************************************************************
;;*
;;*          International Young Astronomer School on Exploiting
;;*                   Herschel and Planck Data (2013)
;;*
;;*                   STEP5: Modified black body fits
;;*                       (Created by F. Galliano)
;;*
;;*****************************************************************************
;;

;; Modified black body function
;;-----------------------------
;; This function will need to be updated later, since the HAPER results are in Jy (not MJy/sr).
FUNCTION modBB, wave, temperature, beta
   
  wmic = wave*!MKS.micron ;; [m]
  Bnu = 2.D*!MKS.hplanck*!MKS.clight/wmic^3 $
      / ( EXP(!MKS.hplanck*!MKS.clight/(wmic*!MKS.kboltz*temperature)) - 1.D )
  wave0 = 100.D
  Inu = ((wave0/wave)^beta)*Bnu*10^(20.) ;; [MJy/sr]

  RETURN, Inu

END

;; Interface with the fitter
;;--------------------------
PRO fitinterface, wave, parm, Inu
  Inu = EXP(parm[0]) * MODBB(wave,EXP(parm[1]),parm[2])
END

  ;;==========================================================================

;; 1) Load the data
;;-----------------
AME = READ_CSV('../Data/AME.txt', COUNT = amy, HEADER = amyhead, MISSING_VALUE='')
print, "AME Data Read"
NROIs = N_ELEMENTS(AME.field01)

RESTORE()

 ;;Arrays for the Big Loop on the Regions
temperature_all    = DBLARR(1,ns)
beta_all           = DBLARR(1,ns)
G0_all             = DBLARR(1,ns)
chi2_all          = DBLARR(1,ns)
FIR_all            = DBLARR(1,ns)
tau250_all         = DBLARR(1,ns)
bands_FIR_all      = DBLARR(nmaps,ns)
bands_G0_all       = DBLARR(nmaps,ns)
Inu_SED_all          = DBLARR(11,ns)
pixct_all               = DBLARR(nmaps,ns)
Snu_SED_all         = DBLARR(nmaps,ns)



  ;; 2) Here's where we setup the Big Loop on the Regions

FOR cerberus=0,NROIs-1 DO BEGIN
 ROI  = AME.field01[cerberus]
 RESTORE,"../Save/noirc_wAME/step4_"+ROI+".xdr"
weights = 1./dLnu_SED^2
Inu_SED = Lnu_SED


;; We restrain the wavelength range by zeroing the weight of the
;; fluxes < 50 microns.
 
 weights[*,*,WHERE(wave LT 90.D) AND WHERE(wave GE 1000.D)] = 0.D
; weights[*,*,WHERE(wave LT 90)] = 0.D

;; 2) BB fitting
;;----------------------

       ;;Arrays for the big loop on the pixels

      ;;   a. Initial guesses of the parameters:
      ;;      parm[0] = LOG(optical depth)
      ;;      parm[1] = LOG(temperature in K)
      ;;      parm[2] = beta

      Nparm = 3
      parm = DBLARR(Nparm)
      parm[2] = 2.D
      Bnumax =  MAX(Inu_SED[x,y,*]*nu^(-parm[2]),imax)
     parm[1] = ALOG( 5.1D-3 / (wave[imax]*!MKS.micron) ) ;; relation T/lambmax
;      parm[1] = ALOG( 20.D) ;; Temperature
      parm[0] = ALOG( Inu_SED[x,y,imax]/MODBB(wave[imax],EXP(parm[1]),parm[2]) )

      ;;   b. Call the least square fitter.
      parinfo = REPLICATE({value:0.0,fixed:0,limited:[0,0],limits:[0.,0.]},Nparm)
;      parinfo[2].limited = [1,1]  
      parinfo[2].limits = [0,5]
;       parinfo[2].fixed = 1
;       parinfo[2].value = 2.D

      fit = MPCURVEFIT( wave, Inu_SED[x,y,*], weights[x,y,*], parm, $
                        FUNCTION_NAME="fitinterface", PARINFO=parinfo, $
                        /NODERIVATIVE, /QUIET )

      ;;b.b Calculate color correction factors based on results of the first fitting

      Nfine    = 500
      wfine    = RAMP(Nfine,1.,1000,/POW) ;; Create a logarithmic ramp (private function)
      nufine   = !MKS.clight/!MKS.micron/wfine
;      filters  = ['IRAS3','AKARI3','AKARI4','IRAS4','AKARI5','AKARI6','HFI1','HFI2'] ;; 
      filters  = ['AKARI4','IRAS4','AKARI5','AKARI6','HFI1','HFI2'] ;;
      ;sed_cc   = dustem_cc (wfine, nufine*MODBB(wfine,EXP(parm[1]),parm[2]) $
       ;           *EXP(parm[0]),filters, cc=cc)
      sed_cc   = dustem_cc(wfine, MODBB(wfine,EXP(parm[1]),parm[2])*EXP(parm[0]),filters, cc=cc)
      print, cc

      ;;b.c. Apply color correction factors to the data, then re-run the fitting.
      ;; The FOR loop here is offset by 4. This is because the color correction is only applied to the FIR bands. We skip the first 4 bands.
      
      FOR h = 4, 9 DO BEGIN
      Inu_SED[x,y,h] = Inu_SED[x,y,h] / cc[h-4]

      ENDFOR

      fit = MPCURVEFIT ( wave, Inu_SED[x,y,*], weights[x,y,*], parm, $
                        FUNCTION_NAME="fitinterface", PARINFO=parinfo, $
                        /NODERIVATIVE, /QUIET )

      ;;   c. Store the final results.
      PCXVG0             = (17.5D)^(4+2)


      G0_corr            = 10^(2.8709-(1.4267*beta[x,y]))  ;;Correction factor when using a free beta        

      tau250[x,y]        = EXP(parm[0])
      temperature[x,y]   = EXP(parm[1])
      beta[x,y]          = parm[2]
      chi2[x,y]          = TOTAL(weights[x,y,*]*(Inu_SED[x,y,*]-fit)^2)/(Nband-Nparm-1.)
      modbb_fine         = MODBB(wfine,temperature[x,y],beta[x,y])*tau250[x,y]
      modbb_bands        = MODBB(wave,temperature[x,y],beta[x,y])*tau250[x,y]
      FIR[x,y]           = integral(wfine,modbb_fine,5,1000, /Double)       ;; Total Far-Infrared Emission
    
  G0[x,y]            = G0_corr*((temperature[x,y]/17.5D)^(4+beta[x,y])) ;; ISRF (Relative to Solar Village)

.      IRC_ratio[x,y]     = Inu_SED[x,y,0] / Inu_SED[x,y,2]                  ;; 9:18 micron ratio
      bands_FIR[x,y,*]   = Inu_SED[x,y,*] / FIR[x,y]                        ;; Each band's intensity vs. FIR
      bands_G0[x,y,*]    = Inu_SED[x,y,*] / G0[x,y]                         ;; Each band's intensity vs. the starlight field

      ;;   d. Print the results.

      PRINT, "Pixel "+STRTRIM(x+1,2)+"/"+STRTRIM(Nx_SED,2)+"," $
                     +STRTRIM(y+1,2)+"/"+STRTRIM(Ny_SED,2)
      PRINT, "  tau250 = "        , tau250[x,y],        " [m]"
      PRINT, "  temperature = " , temperature[x,y], " K"
      PRINT, "  G0 =  "         , G0[x,y],          " ISRF/ISRF_local"
      PRINT, "  beta = "        , beta[x,y]
      PRINT, "  chi2 = "        , chi2[x,y]

    ENDIF
  ENDFOR
ENDFOR

;; 3) Analysis and savings
;;------------------------
;;   a. Inspect a couple of fits. 

SovN = REFORM( Inu_SED[*,*,WHERE(bands EQ "planck857")] $
               / Inu_SED[*,*,WHERE(bands EQ "iras12")] )
;; Plot SEDs ordered by S/N at 250 microns.
isort = SORT(SovN)
isort = isort[WHERE(mask_SED[isort] EQ 0,Nsort)]
Nplot = 8
iplot = isort[FLOOR(DINDGEN(Nplot)/(Nplot-1.)*(Nsort-1.))]
Inu = DBLARR(Nband,Nplot)
FOR i=0,Nband-1 DO Inu[i,*] = (Inu_SED[*,*,i])[iplot]
;; Fine wavelength grid to display the model.
Nfine = 200
wfine = RAMP(Nfine,1.,1000,/POW) ;; Create a logarithmic ramp (private function)
nufine = !MKS.clight/!MKS.micron/wfine


;;Calculate some average properties for each map. 

       IRC_ratio_avg               = TOTAL(Inu_SED[*,*,0], /NaN) / TOTAL(Inu_SED[*,*,2], /NaN) ;Average of each band vs. FIR
       bands_FIR                   = MEAN(MEAN(bands_FIR,DIMENSION=1,/NAN),DIMENSION=1,/NAN)

         FOR catniss = 0, 10 DO BEGIN       
         Inu_SED_each        = DBLARR(Nx_SED,Ny_SED)
         Inu_SED_each      = Inu_SED[*,*,catniss]
         Inu_SED_each_mean = MEAN(Inu_SED_each[WHERE(mask_SED EQ 0)],/NaN)       
         Inu_SED_all[catniss,cerberus] = Inu_SED_each_mean
         ENDFOR
         IRC_ratio_all[*,cerberus]         = IRC_ratio_avg
         temperature_all[cerberus]       = MEAN(temperature[WHERE(mask_SED EQ 0)],/NaN)
         G0_all[cerberus]                      = MEAN(G0[WHERE(mask_SED EQ 0)], /NaN)
         tau250_all[cerberus]                = MEAN(tau250[WHERE(mask_SED EQ 0)],/NaN) 
         beta_all[cerberus]                    = MEAN(beta[WHERE(mask_SED EQ 0)],/NaN)
         pixct_all[cerberus]                    = N_ELEMENTS(G0_all[WHERE(MASK_SED EQ 0)])
         solidangle_all[cerberus]           = pixct_all[cerberus] * ((1.7*1.7)/ 4.25D10)
         Snu_SED_all[*,cerberus]  = Inu_SED_all[*,cerberus] * solidangle_all[cerberus]


;;   b. Display the maps of the physical components.
;SEVERALIMAGES, { im1:tau250, tit1:"!6Optical Depth (250um) [m/pixel]", $
;                 im2:temperature, tit2:"!6Dust temperature [K]", $
;                 im3:beta, tit3:"!6Emissivity index !7b!6", $
;                 im4:chi2, tit4:"!6Reduced chi square", $
;                 im5:FIR, tit5:"Total FIR", $
;                 im6:G0, tit7:"G0"}, $ 
;               /RAINBOW, /SAMPLE, WIN=1


;;;;;;;;;;;;;;;;;;;;;A Little Experiment To Display the Images with Colorbar via Coyote Graphics;;;;;;;;;;;;;;;;;
FOR k=0,1 DO BEGIN
    IF (k EQ 1) THEN BEGIN

fileps = "../Figures/noirc_wAME/appA_"+STRTRIM(STRING(cerberus),1)+".eps"
  SET_PLOT, 'PS'
  DEVICE, FILE=fileps, /ENCAPS, /COLOR, BITS=24, XSIZE=6.0, YSIZE=1.50, /INCH
ENDIF


 !P.Multi =0
IF (k EQ 1) THEN BEGIN
DEVICE, /CLOSE
cgPS2PDF, fileps,  /DELETE_PS, UNIX_CONVERT_CMD='epstopdf'
SET_PLOT, 'X'
ENDIF
ENDFOR

phot_error = [0.051D,0.151D,0.104D,0.10D,0.10D,0.135D,0.10D,0.10D,0.07D,0.07D,0.10D]
Inu_SED_all_error = DBLARR(11,98)

FOR catniss = 0, 10 DO BEGIN
Inu_SED_all_error[catniss,cerberus] = phot_error[catniss] * Inu_SED_all[catniss,cerberus]
ENDFOR

Y_error = Inu_SED_all_error[*,cerberus]
X_error = [4.10D,3.5D,9.97D,5.5D,20.0D,10.85D,19.D,18.D,26.2D,17.5D, 120D]
;X_error = wave[*] ;0.0001D

;;Plot the average SED of each region
;; Plot SEDs ordered by S/N at 250 microns.
;; Fine wavelength grid to display the model.
Nfine = 200
wfine = RAMP(Nfine,1.,1000,/POW) ;; Create a logarithmic ramp (private function)
nufine = !MKS.clight/!MKS.micron/wfine
;; Actual display
FOR k=0,1 DO BEGIN
 @get_colors
  col = [!RED,!GREEN,!YELLOW,!BLUE,!RED,!GREEN,!YELLOW,!BLUE,!RED,!GREEN,!YELLOW,!BLUE]
  IF (k EQ 1) THEN BEGIN
    fileps = "../Figures/appB_"+STRTRIM(STRING(cerberus),1)+".eps"
    SET_PLOT, 'PS'
    DEVICE, FILE=fileps, /ENCAPS, /COLOR, BITS=8, XSIZE=7, YSIZE=5, /INCH
  ENDIF ELSE WINDOW, 2, XSIZE=600, YSIZE=400
  PLOT_OO, [0.], [0.], XRANGE=[8,1000], YRANGE=[1D-2,2D3], $ 
           XTITLE="!6Wavelength [!7l!6m] ", $
           YTITLE="I("+GREEK('lambda')+") [MJy/sr]", /XSTYLE, /YSTYLE, THICK=2, $
           TITLE=" "+ROI+" ", $
           YTICKFORMAT='exponent',$
           CHARSIZE=1.5

  FOR i=0, 10 DO BEGIN
    OPLOT, wfine, MODBB(wfine,temperature_all[cerberus],beta_all[cerberus]) $
                  *tau250_all[cerberus], COLOR=1, THICK=2
    OPLOT, wave[*], Inu_SED_all[*,cerberus], COLOR=col[3], PSYM=4   ;, SYMSIZE = 5
     OPLOTERROR, wave[*], Inu_SED_all[*,cerberus], X_error, Y_error , COLOR=col[3], PSYM=3; X_error
;PRINT, "Pixel "+STRTRIM(x+1,2)+"/"+STRTRIM(Nx_SED,2)+"," $
;                    +STRTRIM(y+1,2)+"/"+STRTRIM(Ny_SED,2)
     XYOUTS, 100, 0.4," "+textoidl('\tau_{250}')+" = "+STRMID(STRTRIM(STRING(tau250_all[cerberus]),1), 0,8)+"[m]", CHARSIZE = 1.5
     XYOUTS, 100, 1, "  T = "+STRMID(STRTRIM(STRING(temperature_all[cerberus]),1),0,5)+"[K]", CHARSIZE = 1.5
     XYOUTS, 100, 2.5, "  G0 = "+STRMID(STRTRIM(STRING(G0_all[cerberus]),1),0,4)+ " ", CHARSIZE = 1.5
;     XYOUTS, 1,1, "  Beta = "                , beta_all[cerberus]
;     XYOUTS, "  chi2 = "                 , MAX(chi2_all[cerberus])
  ENDFOR
  IF (k EQ 1) THEN BEGIN
  
    DEVICE, /CLOSE
     cgPS2PDF, fileps,  /DELETE_PS, UNIX_CONVERT_CMD='epstopdf'
    SET_PLOT, 'X'
    PRINT, " - File "+fileps+" has been written."
  ENDIF
ENDFOR




;;   d. Savings.
;; Save the individual results of each ROI
filexdr = "../Save/noirc_wAME/step5_betafix_"+ROI+".xdr"
SAVE, FILE=filexdr, $
tau250,             $ 
temperature,        $
beta,               $
G0,                 $
chi2,               $
Nx_SED, Ny_SED,     $
FIR,                $ 
Nband,              $
wave,               $
nu,                 $
bands_FIR,          $
bands_G0,           $
cc                 
PRINT, " - File "+filexdr+" has been written."

MWRFITS, G0, "G0_"+ROI+".fits"
MWRFITS, tau250, "tau250_"+ROI+".fits"
MWRFITS, FIR, "FIR_"+ROI+".fits"

ENDFOR 
;;Save the General file, with mean results from all ROIs. Will be used in StepP for plotting.
filexdr = "../Save/noirc_wAME/step5_betafix_general.xdr"
SAVE, FILE=filexdr,  $
Nband,               $
Inu_SED_all,         $
Inu_FIR_all,         $
IRC_ratio_all,       $
temperature_all,     $
beta_all,            $
tau250_all,           $
G0_all,              $
FIR_all,             $
bands_FIR_all,       $
bands_G0_all,      $
solidangle_all,    $
Snu_SED_all 
PRINT, " - File "+filexdr+" has been written."

filexdr = "../Save/noirc_wAME/step5_Inu_SED_all.xdr"
SAVE, FILE=filexdr  
Inu_SED_all

END
