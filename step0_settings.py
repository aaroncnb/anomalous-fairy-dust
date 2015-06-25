
# Set which photometric bands' data are to be used in the analysis
bands = ["iras12","iras25","iras60","akari65","akari90","iras100","akari140","akari160","planck857","planck545"]
# Set which of the bands' data are in HEALPix
# Specify the pixel scale of the input data

reso1 = 1.7

#Specify the width of the desired cutouts, in pixels

size_arcmin = 240

# Set the smoothing width (in pixels) for the Guassian PSF Smoothing
smoothing_size_arcmin = 60

# Set the source aperture radius to be used in the case of circular aperture photometry

apSize = 120 / reso1

# Set the inner and outer radii of the annulus to be used for the background subtraction

bgSizeInner = 80/reso1
bgSizeOuter = 100/reso1

