inmap = "/work1/users/aaronb/Databrary/HEALPix/AKARI_HEALPix_orig/1024_nside/planck_857_1_1024_1dres.fits"
freq  = 857
res_arcmin = 60.
lon = 160.60
lat = -12.05
aper_inner_radius = 0.0174533
aper_outer_radius1 = 0.0349066
aper_outer_radius2 = 0.0523598
units = 'MJy/sr'

fd,fd_err,fd_bg = \
    haperflux(inmap, freq, res_arcmin, lon, lat, aper_inner_radius, \
        aper_outer_radius1, aper_outer_radius2, units, fd=0, fd_err=0, fd_bg=0, \
        column=0, dopol=False, nested=False, noise_model=0, \
        centroid=0)
    
print fd, fd_err, fd_bg
