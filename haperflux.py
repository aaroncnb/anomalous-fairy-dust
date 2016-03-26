def planckcorr(freq):
    h = 6.62606957E-34
    k = 1.3806488E-23 
    T = 2.725
    x = h*(freq*1E9)/k/T
    return (exp(x) - 1)**2/x**2/exp(x) 

def haperflux(inmap, freq, res_arcmin, lon, lat, aper_inner_radius, \
    aper_outer_radius1, aper_outer_radius2, units, fd=0, fd_err=0, fd_bg=0, \
    column=0, dopol=False, nested=False, noise_model=0, \
    centroid=0):


    #check parameters
    if len(sys.argv) > 8:
        print ''
        print 'SYNTAX:-'
        print ''
        print 'HIDL>haperflux(inmap, freq, res_arcmin, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units, fd, fd_err, fd_bg,column=column, noise_model=noise_model, /dopol, /nested)'
        print ''
        exit()
    

    #set parameters
    inmap = inmap
    thisfreq = float(freq)
    lon = float(lon)
    lat = float(lat)
    aper_inner_radius  = float(aper_inner_radius)
    aper_outer_radius1 = float(aper_outer_radius1)
    aper_outer_radius2 = float(aper_outer_radius2)

    # read in data
    s = np.size(inmap)

    if (s == 7):
        hmap,hhead = hp.read_map(inmap, hdu=1) #Check if Ring or Nested is needed

    if (s==4) or (s==5):
        map = inmap
        inmap = ''
        if (nest==True):
            ordering='RING' 
        else:
            ordering = 'NESTED'
    nside = np.sqrt(len(hmap[:,0])/12)
    if (round(nside,1)!=nside) or ((nside%2)!=0):
        print ''
        print 'Not a standard Healpix map...'
        print ''
        exit()
      

    npix = 12*nside**2
    ncolumn = len(map[0,:])

# set column number and test to make sure there is enough columns in
# the file!
    if (field == 0):
        column = 0 
    
    else:
        column=round(column,1)
        
    if (((column+1) > ncolumn) and (dopol== 0)):
        print ''
        print 'Column number requested larger than the number of columns in the file!'
        print ''
        exit()
    

# check for dopol keyword for calculating polarized intensity
    if (dopol==0):
        dopol = 0
    else:
        dopol = 1
        column = 1
        if (ncolumn < 3):
            print ''
            print 'To calculate polarized intensity (PI), requires polarization data with 3 columns or more...'
            print ''
            exit()

#-----do the centroiding here

    if (centroid==True):
        print 'Doing Re-centroiding of coordinates'
    

# get pixels in aperture
    phi = lon*np.pi/180.
    theta = np.pi/2.-lat*np.pi/180.
    vec0=hp.ang2vec(theta, phi)
        
    if (ordering == 'NESTED'):
        # According to the HP git repository- hp.query_disc is faster in RING
        innerpix, ninnerpix = hp.query_disc(nside=nside, vec=vec0, radius=aper_inner_radius/60, nest=True)
        #query_disc, nside, vec0, aper_inner_radius/60., innerpix, ninnerpix, /deg, /nest
        #query_disc, nside, vec0, aper_outer_radius1/60., outerpix1, nouterpix1, /deg, /nest 
        #query_disc, nside, vec0, aper_outer_radius2/60., outerpix2, nouterpix2, /deg, /nest 
    #else:
        #query_disc, nside, vec0, aper_inner_radius/60., innerpix, ninnerpix, /deg
        #query_disc, nside, vec0, aper_outer_radius1/60., outerpix1, nouterpix1, /deg 
        #query_disc, nside, vec0, aper_outer_radius2/60., outerpix2, nouterpix2, /deg 
    


# do not include pixels that are bad
    good0 = np.where(hmap[innerpix,column] != hp.UNSEEN)
        
    #good0 = where(map[innerpix,column] NE !healpix.bad_value, ninnerpix, $
    #            complement=bad0, ncomplement=nbad0)
    #good1 = where(map[outerpix1,column] NE !healpix.bad_value, nouterpix1, $
    #            complement=bad1, ncomplement=nbad1)
    #good2 = where(map[outerpix2,column] NE !healpix.bad_value, nouterpix2, $
    #            complement=bad2, ncomplement=nbad2)

    if (ninnerpix == 0) or (nouterpix1 == 0) or (nouterpix2 == 0):
        print ''
        print '***No good pixels inside aperture!***'
        print ''
        fd = np.nan
        fd_err = np.nan
        exit()
        #GOTO, SKIP1 - I think "GOTO" was used here because the IDL code's authors wanted to exit the code, 
        #     but keep the intermediate values    

    innerpix = innerpix[good0]
    outerpix1 = outerpix1[good1]
    outerpix2 = outerpix2[good2]

 

# find pixels in the annulus (between outerradius1 and outeradius2) 
    outerpix_all = [outerpix1,outerpix2]
    outerpix_all = outerpix_all[sort(outerpix_all)]
    outerpix = -1L
        
    for i in range(0L, nouterpix2):
        temp = outerpix1.index(outerpix2[i])
        if (temp[0] < 0):
            outerpix = [outerpix,outerpix2[i]]
    
    outerpix = outerpix[1:]
    nouterpix = len(outerpix)

# get conversion from to Jy/pix
    pix_area = 4.*np.pi / npix
    factor = 1.0

    if (units == 'K') or (units == 'K_RJ') or (units == 'KRJ'):
        factor = 2.*1381.*(thisfreq*1.0e9)**2/(2.997e8)**2 * pix_area
    
    if (units == 'mK') or (units == 'mK_RJ') or (units == 'mKRJ'):
        factor = 2.*1381.*(thisfreq*1.0e9)**2/(2.997e8)**2 * pix_area / 1.0e3
    
    if (units == 'uK') or (units == 'uK_RJ') or (units == 'uKRJ'):
        factor = 2.*1381.*(thisfreq*1.0e9)**2/(2.997e8)**2 * pix_area / 1.0e6
        
    if (units == 'K_CMB') or (units == 'KCMB'):
        factor = 2.*1381.*(thisfreq*1.0e9)**2/(2.997e8)**2 * pix_area / planckcorr(thisfreq)
        
    if (units == 'mK_CMB') or (units == 'mKCMB'):
        factor = 2.*1381.*(thisfreq*1.0e9)**2/(2.997e8)**2 * pix_area / 1.0e3 / planckcorr(thisfreq)
        
    if (units == 'uK_CMB') or (units == 'uKCMB'):
        factor = 2.*1381.*(thisfreq*1.0e9)**2/(2.997e8)**2 * pix_area / 1.0e6 / planckcorr(thisfreq)
        
    if (units == 'MJy/sr') or (units == 'MJY/SR') or (units == 'MjySr'):
        factor = pix_area * 1.0e6
        
    if (units == 'Jy/pixel') or (units == 'JY/PIXEL') or (units == 'JY/PIX') or (units == 'JyPix'):
         factor = 1.0
            
    if (units == 'average') or (units == 'avg') or (units == 'Average') or (units == 'AVG'):
         factor = 1.0 / float(ninnerpix)

# override columns if /dopol keyword is set
    if (dopol == 1):
        ncalc = 2
    
    else:
        ncalc = 1

    for i in range(1L, ncalc):

        if (dopol == 1):
            column=i

# Get pixel values in inner radius, converting to Jy/pix
        fd_jypix_inner = map[innerpix,column] * factor

# sum up integrated flux in inner
        fd_jy_inner = total(fd_jypix_inner)

# same for outer radius but take a robust estimate and scale by area
        fd_jy_outer = median(map[outerpix,column]) * factor

# subtract background
        fd_bg = fd_jy_outer
        fd = fd_jy_inner - fd_bg*float(ninnerpix)


#estimate error based on robust sigma of the background annulus
        if (noise_model == 0):
            noise_model = 0
        noise_model = 1

# new version (2-Dec-2010) that has been tested with simulations for
# Planck early paper and seems to be reasonable for many applications

        if (noise_model == 0): 
            Npoints = (pix_area*ninnerpix) / \
                      (1.13*(float(res_arcmin)/60. *np.pi/180.)**2)
            Npoints_outer = (pix_area*nouterpix) / \
                      (1.13*(float(res_arcmin)/60. *np.pi/180.)**2)
            fd_err = stddev(map[outerpix,column]) * factor * ninnerpix / sqrt(Npoints)
      

 # works exactly for white uncorrelated noise only!
        if (noise_model == 1):
            k = np.pi/2.

            fd_err = factor * sqrt(float(ninnerpix) + \
                (k * float(ninnerpix)**2/nouterpix)) * robust_sigma(map[outerpix,column])
        

# if dopol is set, then store the Q estimate the first time only
        if(dopol == 1) and (i == 1):
            fd1 = fd
            fd_err1 = fd_err
            fd_bg1 = fd_bg
        
    

# if dopol is set, combine Q and U to give PI
    if (dopol == 1): 
        fd = sqrt(fd1**2 +fd**2)
        fd_err = sqrt( 1./(fd1**2 + fd**2) * (fd1**2*fd_err1**2 + fd**2*fd_err** 2))
        fd_bg = sqrt(fd_bg1**2 + fd_bg**2)
    

#SKIP1: 
    return fd, fd_err, fd_bg


