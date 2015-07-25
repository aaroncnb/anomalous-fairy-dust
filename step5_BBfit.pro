;;PRO step5_BBfit_akari_
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
  Snu = ((wave0/wave)^beta)*Bnu*10^(20.)  ;;*[Sr per one degree aperture]/[10^6] ;; [Jy]

  RETURN, Snu

END

;; Interface with the fitter
;;--------------------------
PRO fitinterface, wave, parm, Snu
  Snu = EXP(parm[0]) * MODBB(wave,EXP(parm[1]),parm[2])
END

  ;;==========================================================================

;; 1) Load the data
;;-----------------
RESTORE, "../Data/multiepoch_photometry_akari_.sav"
print, "circular aperture photometry variables and results restored..."

 ;;Arrays for the Big Loop on the Regions
temperature_all    = DBLARR(1,ns)
beta_all           = DBLARR(1,ns)
G0_all             = DBLARR(1,ns)
chi2_all          = DBLARR(1,ns)
FIR_all            = DBLARR(1,ns)
tau250_all         = DBLARR(1,ns)
bands_FIR_all      = DBLARR(nmaps,ns)
bands_G0_all       = DBLARR(nmaps,ns)
Snu = fd_all
weights = fd_err_all
wave = [550.,345.,160.,140.,100.,90.,65.,60.,25.,12]


  ;; 2) Here's where we loop through the sources
for s=0,ns-1 do begin

   ;; We restrain the wavelength range by zeroing the weight of the
   ;; fluxes < 50 microns.
   weights[s,WHERE(wave LT 90)] = 0.D
   ;; 2) BB fitting
   ;;----------------------
   ;;This is done on a per-row basis in the case of circular aperture photometry results
      ;;   a. Initial guesses of the parameters:
      ;;      parm[0] = LOG(optical depth)
      ;;      parm[1] = LOG(temperature in K)
      ;;      parm[2] = beta

   Nparm = 3
   parm = DBLARR(Nparm)
   parm[2] = 2.D
   Bnumax =  MAX(Snu[s,*]^(-parm[2]),imax)
   ;;Note that "imax" is the subscript of the maximum, not the maximum itself...be careful!
   parm[1] = ALOG( 5.1D-3 / (wave[imax]*!MKS.micron) ) ;; relation T/lambmax
   parm[1] = ALOG( 20.D) ;; Temperature
   parm[0] = ALOG( Snu[s,imax]/MODBB(wave[imax],EXP(parm[1]),parm[2]) )

      ;;   b. Call the least square fitter.
      parinfo = REPLICATE({value:0.0,fixed:0,limited:[0,0],limits:[0.,0.]},Nparm)
   ;      parinfo[2].limited = [1,1]  
   ;      parinfo[2].limits = [0,5]
   ;       parinfo[2].fixed = 1
   ;         parinfo[2].value = 2.D

      fit = MPCURVEFIT( wave, Snu[s,*], weights, parm, $
                        FUNCTION_NAME="fitinterface", PARINFO=parinfo, $
                        /NODERIVATIVE, /QUIET )

      ;;b.b Calculate color correction factors based on results of the first fitting

      Nfine    = 500
      wfine    = RAMP(Nfine,1.,1000,/POW) ;; Create a logarithmic ramp (private function)
      nufine   = !MKS.clight/!MKS.micron/wfine
;      filters  = ['IRAS3','AKARI3','AKARI4','IRAS4','AKARI5','AKARI6','HFI1','HFI2'] ;; 
      filters  = ['AKARI4','IRAS4','AKARI5','AKARI6','HFI1','HFI2'] ;;
      sed_cc   = dustem_cc(wfine, MODBB(wfine,EXP(parm[1]),parm[2])*EXP(parm[0]),filters, cc=cc)
      print, cc

      ;;b.c. Apply color correction factors to the data, then re-run the fitting.
      ;; The FOR loop here is offset by 4. This is because the color correction is only applied to the FIR bands. We skip the first 4 bands.
      
      FOR h = 4, 9 DO BEGIN
         Snu[s,h] = Snu[s,h] / cc[h-4]
      ENDFOR

      fit = MPCURVEFIT ( wave, Snu[s,*], weights[s,*], parm, $
                        FUNCTION_NAME="fitinterface", PARINFO=parinfo, $
                        /NODERIVATIVE, /QUIET )
      
      ;;   c. Store the final results.
      PCXVG0             = (17.5D)^(4+2)
      G0_corr            = 10^(2.8709-(1.4267*beta[s]))  ;;Correction factor when using a free beta        
      tau250_all[s]        = EXP(parm[0])
      temperature_all[s]   = EXP(parm[1])
      beta_all[s]          = parm[2]
      chi2_all[s]          = TOTAL(weights[s,*]*(Snu[s,*]-fit)^2)/(Nband-Nparm-1.)
      modbb_fine         = MODBB(wfine,temperature_all[s],beta_all[s])*tau250_all[s]
      modbb_bands        = MODBB(wave,temperature_all[s],beta_all[s])*tau250_all[s]
      FIR_all[s]           = integral(wfine,modbb_fine,5,1000, /Double)       ;; Total Far-Infrared Emission
    
  G0            = G0_corr*((temperature/17.5D)^(4+beta)) ;; ISRF (Relative to Solar Village)

      bands_FIR[s,*]   = Snu[s,*] / FIR                        ;; Each band's intensity vs. FIR
      bands_G0[s,*]    = Snu[s,*] / G0                         ;; Each band's intensity vs. the starlight field

      ;;   d. Print the results.

      PRINT, "Pixel "+STRTRIM(x+1,2)+"/"+STRTRIM(Nx_SED,2)+"," $
                     +STRTRIM(y+1,2)+"/"+STRTRIM(Ny_SED,2)
      PRINT, "  tau250 = "        , tau250,        " [m]"
      PRINT, "  temperature = " , temperature, " K"
      PRINT, "  G0 =  "         , G0,          " ISRF/ISRF_local"
      PRINT, "  beta = "        , beta
      PRINT, "  chi2 = "        , chi2

  ENDFOR


;; 3) Analysis and savings
;;------------------------
;;   a. Inspect a couple of fits. 


phot_error = [0.051D,0.151D,0.104D,0.10D,0.10D,0.135D,0.10D,0.10D,0.07D,0.07D,0.10D]

;;   d. Savings.

;;Save the General file, with mean results from all sources. Will be used in plotting.
filexdr = "../Save/BBfit.xdr"
SAVE, FILE=filexdr,  $
Nband,               $
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

END
