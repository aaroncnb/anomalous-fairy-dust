# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from numpy import genfromtxt

fd = genfromtxt('fd.csv',delimiter=',')
lnu = genfromtxt('lnu.csv',delimiter=',')
target_names = genfromtxt('inputlist_mep_test.list',usecols = (0), dtype = 'string', delimiter=',')
wl = [9,65,90,140,160,353,550] #IRAS data not used. AKARI 18 not used.
x = wl



def plot_phot(target):
	#Set the circular aperture results as y1
	y1 = fd[target,(2,4,6,8,10,12,14)] #Only includes AKARI 9, 65, 90, 140, 160 and Planck 857 and 545 - be careful!
	#Set the earlier pixel-averaged reults as y2
	y2 = lnu[target,(0,10,12,16,18,20,22)] #The indexes are different from y1 because we skip the IRAS and AKARI 18 data
	#Scale each data set according to the longest wavelength point (Planck 545 GHz) - for quick comparison of the SED shapes
	y1_s = y1 / fd[target,14]
	y2_s = y2 / lnu[target,22]
	
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.scatter(wl,y1_s, s=30, c='b', marker="s", label='Circ. Aper.')
	ax1.scatter(wl,y2_s, s=30, c='r', marker="o", label='Pixel-Avg.' )
	ax1.set_title(target_names[target])
	plt.ylabel('Flux / Flux @ 550 microns GHz')
	plt.xlabel('Wavelength (Microns)')
	plt.grid(True)
	plt.legend(loc='upper left');
	plt.show()
	plt.savefig('Plots/'+str(target)+'.pdf')
	plt.close()

for target in range(0,98):
	plot_phot(target)

#subprocess.Popen("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=Plots/CircAper_vs_PixAvg.pdf Plots/*.pdf")


