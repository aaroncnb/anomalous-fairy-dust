import numpy as np
import healpy as hp
import sys



#http://stackoverflow.com/questions/14275986/removing-common-elements-in-two-lists
# Here is a function to remove the "common" pixels of the innter and outer rings of the background annulus
# The point is that the background-ring pixels we want, are the ones that /not/ contained by both outer rings.
# If we were calculating the rings's area, we'd subtract the innter ring area from the outer. It's essentially the same logic here
# I found an example function on stackoverflow, by "user1632861 "
# It's intended for lists, rather than numpy arrays, but I think it should work

def removeCommonElements(outerpix1, outerpix2):
    for pix in outerpix2[:]:
        if pix in outerpix1:
            outerpix2.remove(pix)
            outerpix1.remove(pix)

def removeCommonElementsNumpy(outerpix1, outerpix2):
    outpix1 = outerpix1.copy()
    outpix2 = outerpix2.copy()
    for pix in outerpix2[:]:
        if pix in outerpix1:
            np.delete(outpix1,pix)
            np.delete(outpix2,pix)
    return outpix1, outpix2

def deleteCommon(outerpix1,outerpix2):
    outerpix = np.delete(outerpix2, outerpix1)
    return outerpix


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

    if (s == 1):
        
        hmap,hhead = hp.read_map(inmap, hdu=1,h=True, nest=nested) #Check if Ring or Nested is needed
        #http://healpy.readthedocs.org/en/latest/generated/healpy.fitsfunc.read_map.html
        print np.size(hmap)
    
    if (s>1):
        hmap = inmap
        inmap =  ''
    
    if (nested==False):
        ordering='RING' 
    else:
        ordering ='NESTED'
    
    nside = np.sqrt(len(hmap)/12)
    if (round(nside,1)!=nside) or ((nside%2)!=0):
        print ''
        print 'Not a standard Healpix map...'
        print ''
        exit()
      

    npix = 12*nside**2
    ncolumn = len(hmap)

# set column number and test to make sure there is enough columns in
# the file!
    if (column == 0):
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
        
        # According to the HP git repository- hp.query_disc is faster in RING
        
    ## Get the pixels within the innermost (source) aperture
    innerpix = hp.query_disc(nside=nside, vec=vec0, radius=aper_inner_radius, nest=nested)
    #ninnerpix = len(innerpix) #since the HEALPY version of query_disc doesn't return the number of pixels
        
    ## Get the pixels within the inner-ring of the background annulus
    outerpix1 = hp.query_disc(nside=nside, vec=vec0, radius=aper_outer_radius1, nest=nested)
    #nouterpix1 = len(outerpix1)
        
    ## Get the pixels within the outer-ring of the background annulus
    outerpix2 = hp.query_disc(nside=nside, vec=vec0, radius=aper_outer_radius2, nest=nested)
    #nouterpix2 = len(outerpix2)
        
    


# Identify and remove the bad pixels
# In this scheme, all of the bad pixels should have been labeled with HP.UNSEEN in the HEALPix maps
    
    bad0 = np.where(hmap[innerpix] == hp.UNSEEN)
    innerpix_masked = np.delete(innerpix,bad0)
    ninnerpix = len(innerpix_masked)
    
    bad1 = np.where(hmap[outerpix1] == hp.UNSEEN)
    outerpix1_masked = np.delete(outerpix1,bad1)
    nouterpix1 = len(outerpix1_masked)
    
    bad2 = np.where(hmap[outerpix2] == hp.UNSEEN)
    outerpix2_masked = np.delete(outerpix2,bad2)
    nouterpix2 = len(outerpix2_masked)
    
    if (ninnerpix == 0) or (nouterpix1 == 0) or (nouterpix2 == 0):
        print ''
        print '***No good pixels inside aperture!***'
        print ''
        fd = np.nan
        fd_err = np.nan
        exit()
        
    innerpix = innerpix_masked
    outerpix1 = outerpix1_masked
    outerpix2 = outerpix2_masked

    # find pixels in the annulus (between outerradius1 and outeradius2) 
    # In other words, remove pixels of Outer Radius 2 that are enclosed within Outer Radius 1

    bgpix = np.delete(outerpix2, outerpix1)
    print "Common Elements Removed"
    
    nbgpix = len(bgpix)

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
        fd_jypix_inner = hmap[innerpix] * factor

        # sum up integrated flux in inner
        fd_jy_inner = total(fd_jypix_inner)

        # same for outer radius but take a robust estimate and scale by area
        fd_jy_outer = median(hmap[bgpix]) * factor

        # subtract background
        fd_bg = fd_jy_outer
        fd = fd_jy_inner - fd_bg*float(ninnerpix)


        #estimate error based on robust sigma of the background annulus
        # new version (2-Dec-2010) that has been tested with simulations for
        # Planck early paper and seems to be reasonable for many applications

        if (noise_model == 0): 
            Npoints = (pix_area*ninnerpix) / \
                      (1.13*(float(res_arcmin)/60. *np.pi/180.)**2)
            Npoints_outer = (pix_area*nbgpix) / \
                      (1.13*(float(res_arcmin)/60. *np.pi/180.)**2)
            fd_err = stddev(hmap[bgpix]) * factor * ninnerpix / sqrt(Npoints)
      

         
        # works exactly for white uncorrelated noise only!
        if (noise_model == 1):
            k = np.pi/2.

            fd_err = factor * sqrt(float(ninnerpix) + \
                (k * float(ninnerpix)**2/nbgpix)) * robust_sigma(hmap[bgpix])
        

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
    

    return fd, fd_err, fd_bg

