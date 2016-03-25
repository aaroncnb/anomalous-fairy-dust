def healpix_phot(inputlist, maplist, radius, galactic=True, decimal=True, rinner=2.0, router=3.0):
    ##Here is the "observation data structure", which just means "a bunch of details 
    ## about the different all-sky data sources which could be used here.
    freqlist =     ['30',     '44',      '70',    '100','143','217','353','545','857','1874','2141','2998','3331','4612','4997','11992','16655','24983','33310']
    freqval =      [28.405889, 44.072241,70.421396,100.,143.,217.,353.,545.,857.,1874.,2141.,2998.,3331.,4612.,4997.,11992.,16655.,24983.,33310.]
    fwhmlist =     [33.1587,28.0852,13.0812,9.88,7.18,4.87,4.65,4.72,4.39,0.86,0.86,4.3,0.86,0.86,4,3.8,0.86,3.8,0.86] # fwhm in arcminutes
    band_names =   ["iras60","akari65","akari90","iras100","akari140","akari160","planck857", "planck545"]
    band_centers = [ 60e-6,    65e-6,    90e-6,   100e-6,   140e-6,    160e-6,    350e-6,      550e-6]
    idl = pidly.IDL()
    

    k0 = 1.0 
    k1 = rinner 
    k2 = router 
    apcor = ((1 - (0.5)**(4*k0**2))-((0.5)**(4*k1**2) - (0.5)**(4*k2**2)))**(-1)
  
    # 'galactic' overrules 'decimal' 
    if (galactic==True):
        targets = np.genfromtxt(inputlist, names='sname,glon,glat',delimiter=",")
        #idl.pro('euler', targets['glon'], targets['glat'], targets.['ra'], targets['dec'],2)
    elif (decimal==True):
        np.genfromtxt(inputlist, names= 'sname,ra,dec',delimiter=",")
        idl.pro('euler',ra,dec,glon,glat,1)
    else:
        np.genfromtxt(inputlist,names='sname,rah,ram,ras,decd,decm,decs',delimiter=",")
        ra = 15*tenv(rah,ram,ras)
        dec = tenv(decd,decm,decs)
        idl.pro('euler',ra,dec,glon,glat,1)
    

    ns = len(targets['glat'])

    fd3 = -1
    fd_err3 = -1

    fn = np.genfromtxt(maplist, delimiter=" ", dtype='str') 
    nmaps = len(fn)
    fd_all = np.zeros((98,16))
    fd_err_all = np.zeros((98,16))
    #openw,1,file_basename(inputlist+'.photo'),width=200

    #if radius == None:
    #    printf 1,'; A multiplicative aperture correction factor of ',apcor,' has been applied to the flux densities'
    #    printf 1, '; assuming that the source is a point source. Flux densities are in Jy.'
    #elif:
    #    printf 1,';No aperture correction factor has been applied. Flux densities are in Jy.'

    #printf 1, ';'
    #printf 1, ';Output format:'
    #printf 1, '; Source_Name  Map_number  GLON   GLAT   Flux (Jy) Flux_RMS (Jy) Median_Background_Flux (Jy)'
    #printf  1, ';'
    #printf  1, '; Map List:'
    # Make a header-legend that says what order the maps were processed:
    #for i = 0, nmaps-1:
    #    printf, 1, '; Map #',i, fn[i]
    #
    #printf, 1, ';'
    # Start the actual processing: Read-in the maps.
    for ct2 in range(0,nmaps):
        xtmp = pyfits.open(fn[ct2])
        freq = xtmp[1].header['FREQ'] #May need to "trim" off the leading space from the FREQ keyword value
        units = xtmp[1].header['TUNIT1'] # Double-check which header index is correct (i.e. [1] or [0]
        idx = np.where(freqlist == freq)
        cnt = len(idx)
        if (cnt > 0):
            currfreq = freqval[idx]
            if (radius == None):
                radval = fwhmlist[idx]
            else:
                radval = radius
        else:
            print 'Invalid frequency ', freq, ' in ', fn[ct2]
            exit()
        
        for ct in range(0,ns-1): 
         
    ##### Let's add an automatic GALACTIC -> CELSTIAL coord conversion
    ##### In case the target list is given in GAL, but the map
    ##### is given in CELESTIAL
            if (currfreq == 33310.) or (currfreq == 16655.):
                  
                idl.pro('EULER', glon[ct], glat[ct], ra, dec,  SELECT = 2)
  
                idl.pro('haperflux', fn[ct2], currfreq, fwhm, ra, dec, 1.*radval, \
                        rinner*radval, router*radval, units, fd, fd_err, fd_bg, nested=True, noise_mod=True)

            else:
                idl.pro('haperflux', fn[ct2], currfreq, fwhm, glon[ct], glat[ct], 1.*radval, \
                        rinner*radval, router*radval, units, fd, fd_err, fd_bg, nested = True, noise_mod=True)
        

            if (np.isfinite(fd_err) == False):
                fd = -1
                fd_err = -1
            else:
                if radius==None:
                    fd = fd*apcor
                    fd_err = fd_err*apcor

  ##########Let's change this next line so that it puts the data into a structure I can use for the blackbody fitting,
  ##########rather than a huge list of photometry results
            #printf,1,sname[ct],ct2,glon[ct],glat[ct],fd,fd_err,fd_bg,format='(A18,X,I2,X,F12.7,X,F12.7,X,E11.3,X,E11.3,X,E11.3)'
            #fd_all[ct,ct2] = fd
            #fd_err_all[ct,ct2] = fd_err
        
    

  #save, /variables, filename='multiepoch_photometry_akari_.sav'
  #save, filename='multiepoch_photometry_akari_.sav'

    #close,1
