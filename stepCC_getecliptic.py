
from astropy import units asu gi
from astropy.coordinates import ICRS, Galactic, Ecliptic
from astropy.table import Table
from astropy.table import Table, Column

gcoords = np.genfromtxt('~/Dropbox/Aaron\/Bell/Databrary/AME_data/AME_gcoords.csv', delimiter =',')

glon  = gcoords[:,1]
glat  = gcoords[:,2]

c = SkyCoord(glon, glat, Galactic, unit="deg")  # defaults to ICRS frame

c.ecliptic
