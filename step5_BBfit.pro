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
  ;;;;The conversion factor from MJy/sr to Jy below is assuming that the circular aperture of the sources has a 1deg radius!
  ;;;; Pi^2 cancels out of the conversion, so in the end we just multiply by 180^2 [deg^-1].
  Snu = ((wave0/wave)^beta)*Bnu*(1.D26)*(0.000304617D) ; [Jy] (  0.000304617D [sr] per Pi^2 [sq.deg]. )

  RETURN, Snu

END

;; Interface with the fitter
;;--------------------------
PRO fitinterface, wave, parm, Snu
  Snu = EXP(parm[0]) * MODBB(wave,EXP(parm[1]),parm[2])
END

  ;;==========================================================================

;; 1.1) Here's where we load in the aperture photometry results from HAPER
;;-----------------
RESTORE, "../Data/multiepoch_photometry_akari_.sav"
print, "circular aperture photometry variables and results restored..."

;; 1.2) Settings and Setup:
;;  Here we create some template arrays. We'll need these at different stages throughout the code. "_all" is there to make it easy to recognize that these are arrays, holding the results of "all" of the sources.

temperature_all    = DBLARR(ns)
beta_all                 = DBLARR(ns)
G0_all                    = DBLARR(ns)
chi2_all                  = DBLARR(ns)
FIR_all                    = DBLARR(ns)
tau100_all              = DBLARR(ns)
Snu_all                   = DBLARR(nmaps,ns)
weights                  = DBLARR(nmaps,ns)

;;List of the wavelengths of the bands used:
wave        =  [12.D,25.D,60.D,65.D,90.D,100.D,140.D,160.D,345.D,550.D]

;;Total photometric errors for each of the maps used
phot_error = [0.051D,0.151D,0.104D,0.10D,0.10D,0.135D,0.10D,0.10D,0.07D,0.07D,0.10D]


;; Take the flux @ a given band from the HAPER result-array, "fd_all"
Snu_all     = fd_all

;;Use the flux-density errors as the weights for this greybody fitting routine. See Planck Coll. Int. Results XV (2013, Clive Dickinson et al.): They add the RMS uncertainty from the background annulus to the photometric errors of the maps.

FOR, serr = 0, ns-1 DO BEGIN
     weights_all[*,serr]    =  fd_err_all[*,serr] + (phot_err[*]*Snu_all[*,serr]
END FOR

;;Exclude the stochastic and small-grain dominated bands from the fitting:
     ;;  This just means that we'll fit a blackbody only to the bands that are dominated by blackbody radiation (i.e., the far-infrared bands longer than 65 um)

weights_all[WHERE(wave LT 90.D),*] = 0.D

;;  Later, we'll try and fit two blackbodies simultaneously- one fit the small grains and one for the bigger grains). The problem with doing that now, is we don't have a full dust SED model implemented. We'll need that to properly estimate the emissivity of the small grains.

;; 2) Here were setup a loop which will cycle through all of the sources in the fd_all array. There's 98 in the full Planck Coll. AME paper. "ns" gives the number of sources, determined earlier by HAPER.

for s=0,ns-1 do begin

;; 2) BB fitting
;;----------------------

      ;;   a. Initial guesses of the parameters:
      ;;      parm[0] = LOG(optical depth)
      ;;      parm[1] = LOG(temperature in K)
      ;;      parm[2] = beta

   Nparm = 3
   parm = DBLARR(Nparm)
   parm[2] = 2.D
   Bnumax =  MAX(Snu_all[WHERE(weights_all[s,*] NE 0.),s]^(-parm[2]),imax)
   ;;Note that "imax" is the subscript of the maximum, not the maximum itself...be careful!
   ;parm[1] = ALOG( 5.1D-3 / (wave[imax]*!MKS.micron) ) ;; relation T/lambmax
   parm[1] = ALOG( 20.D) ;; Temperature
   parm[0] = ALOG(Snu_all[s,imax]/MODBB(wave[imax],EXP(parm[1]),parm[2]) )

      ;;   b. Call the least square fitter.
      parinfo = REPLICATE({value:0.0,fixed:0,limited:[0,0],limits:[0.,0.]},Nparm)
          parinfo[2].limited = [1,1]  
          parinfo[2].limits = [0,5]
   ;;       parinfo[2].fixed = 1
   ;;        parinfo[2].value = 2.D

      fit = MPCURVEFIT( wave, Snu_all[*,s], weights[*,s], parm, $
                        FUNCTION_NAME="fitinterface", PARINFO=parinfo, $
                        /NODERIVATIVE, /QUIET )

      ;;b.b Calculate color correction factors based on results of the first fitting

      Nfine    = 500
      wfine    = RAMP(Nfine,1.,1000,/POW) ;; Create a logarithmic ramp (private function)
      nufine   = !MKS.clight/!MKS.micron/wfine
      filters  = ['AKARI4','IRAS4','AKARI5','AKARI6','HFI1','HFI2']  ;;
      sed_cc   = dustem_cc(wfine, MODBB(wfine,EXP(parm[1]),parm[2])*EXP(parm[0]),filters, cc=cc)
      print, cc

      ;;b.c. Apply color correction factors to the data, then re-run the fitting.
      ;; The FOR loop here is offset by 4. This is because the color correction is only applied to the FIR bands. We skip the first 4 bands.
      
FOR h = 4, 9 DO BEGIN
         Snu_all[h,s] = Snu_all[h,s] / cc[h-4]
      ENDFOR

      fit = MPCURVEFIT (wave, Snu_all[*,s], weights[*,s], parm, $
                        FUNCTION_NAME="fitinterface", PARINFO=parinfo, $
                        /NODERIVATIVE, /QUIET )
      
      ;;   c. Store the final results.
      PCXVG0             =  (17.5D)^(4+2)       
      tau100_all[s]           = EXP(parm[0])
      temperature_all[s]   = EXP(parm[1])

      beta_all[s]          = parm[2]

      chi2_all[s]          = TOTAL(weights[*,s]*(Snu_all[*,s]-fit)^2)/(nmaps-Nparm-1.)

      modbb_fine         = MODBB(wfine,temperature_all[s],beta_all[s])*tau100_all[s]

      modbb_bands        = MODBB(wave,temperature_all[s],beta_all[s])*tau100_all[s]
      
;; Total Far-Infrared Emission
      FIR_all[s]           = integral(wfine,modbb_fine,5,1000, /Double)       
      
      ;; Use G0_corr = 1 when fixing the beta value!
      ;G0_corr            = 1.
      ;; Correction factor when using a free beta 
      G0_corr            = 10^(2.8709-(1.4267*beta_all[s]))  

      ;; ISRF (Relative to Solar Village)
      ;; Use beta_all[s] if estimating G0 based on free-beta fitting
     ;; G0_all[s]            = G0_corr*((temperature_all[s]/17.5D)^(4+beta_all[s])) 
      ;; Use [2] if estimating G0 the PCXV way
       G0_all[s]            = G0_corr*((temperature_all[s]/17.5D)^(4+2)) 

      ;;   d. Print the results.

      PRINT, "  tau100 = "        , tau100_all[s],        " [m]"
      PRINT, "  temperature = " , temperature_all[s], " K"
      PRINT, "  G0 =  "         , G0_all[s],          " ISRF/ISRF_local"
      PRINT, "  beta = "        , beta_all[s]
      PRINT, "  chi2 = "        , chi2_all[s]

  ENDFOR


;; 3) Analysis and savings
;;------------------------
;;   a. Inspect a couple of fits. 

;;   d. Savings.

;;Save the General file, with mean results from all sources. Will be used in plotting.
filexdr = "../Save/BBfit.xdr"
SAVE, FILE=filexdr,  $
nmaps,               $
temperature_all,     $
beta_all,            $
tau100_all,           $
G0_all,              $
FIR_all,             $
Snu_all
PRINT, " - File "+filexdr+" has been written."

END
