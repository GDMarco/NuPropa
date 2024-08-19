import os, sys, glob, math

import matplotlib.pyplot as plt
import numpy as np
# from matplotlib.pyplot import cm
# import itertools
# import matplotlib.patches as patches
# from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
# from matplotlib import font_manager
# from matplotlib.legend import Legend
# import matplotlib.ticker


# Extra the cross-section data as an array
def read_data(file_name):

	output = [ [], [], [] ]

	infile = open(file_name,'r').readlines()

	for iline in range(0,len(infile)):

		line = infile[iline]
		# skip blank lines
		if( not line.strip() ):
			continue
		
		line_data = line.split()

		# skip the information (starts with '#')
		if( line_data[0] == "#" ):
			continue

		# Otherwise extract the cross-section
		for idata in range(0,len(line_data)):
			output[idata].append( float(line_data[idata]) )

	return output

# Sums the cross-section for an array of channels
def sum_sigma(data_in,Emin_int,Emax_int):

	# Number of channels being summed
	nchans = len(data_in)
	# Structure of summed data to be output
	output = [ [],[],[] ]

	for ientry in range(Emin_int,Emax_int):

		sigma_cen = 0.0
		sigma_err2 = 0.0

		Evalue = data_in[0][0][ientry]

		for ichan in range(0,nchans):
			sigma = data_in[ichan][1][ientry]
			sigma_error = data_in[ichan][2][ientry]

			# check for numerical instabilities
			if( math.isnan(sigma) ):
				# take sigma and error from previous entry
				print('Encountered a nan / inf', 'taking cross-section from previous energy value')
				sigma = data_in[ichan][1][ientry-1]
				sigma_error = data_in[ichan][2][ientry-1]


			sigma_cen += sigma
			if( sigma != 0.0 ):
				sigma_err2 += sigma_error**2

		output[0].append( Evalue )
		output[1].append( sigma_cen )
		output[2].append( sigma_err2**0.5 )

	return output

base = 'Ecms_scan_500/'
file_structure = 'SigmaIncl_Ecms_channel'
seed = 's1'

# Look at the muon induced trident channel: three leptonic channels, then count the hadronic twice
channels_trident_CC = ['112','113','117','119','119']
sigma_trident_CC = []
for ichan in range(0,len(channels_trident_CC)):
	file = base+file_structure+channels_trident_CC[ichan]+'_'+seed+'.txt'
	sigma_trident_CC.append( read_data(file) )

# on-shell W
channels_W_CC = ['29']
sigma_W_CC = []
for ichan in range(0,len(channels_W_CC)):
	file = base+file_structure+channels_W_CC[ichan]+'_'+seed+'.txt'
	sigma_W_CC.append( read_data(file) )

# Then additionally compare all the channels: NC and CC tridents
channels_trident = ['101','102','103','104','105','106','107','108','109','112','113','117','119','119']
sigma_trident = []
for ichan in range(0,len(channels_trident)):
	file = base+file_structure+channels_trident[ichan]+'_'+seed+'.txt'
	sigma_trident.append( read_data(file) )

# Compute the ratio of the two predictions, accounting for potential zero division
def compute_ratio( a, b ):
	output = [ [], [], [] ]
	for ientry in range(0,len(a[0])):
		output[0].append( a[0][ientry] )
		if( b[1][ientry] != 0.0 ):
			ratio = a[1][ientry]/b[1][ientry]
			output[1].append( ratio )

			if( ratio != 0.0 ):
				output[2].append( ratio * ( (a[2][ientry]/a[1][ientry])**2 +  (b[2][ientry]/b[1][ientry])**2 ) )
			else:
				output[2].append(0.0)
		else:
			output[1].append(0.0)
			output[2].append(0.0)
	return output

# Compute the total cross-sections (summing over channels)
ebin_max = 480
total_W_CC = sum_sigma( sigma_W_CC, 0, ebin_max )
total_trident_CC = sum_sigma( sigma_trident_CC, 0, ebin_max )
total_trident = sum_sigma( sigma_trident, 0, ebin_max )

# Now compute some ratios
ratio_triCC_tri = compute_ratio( total_trident_CC, total_trident )
ratio_WCC_triCC = compute_ratio( total_W_CC, total_trident_CC )
ratio_WCC_tri = compute_ratio( total_W_CC, total_trident )


### Some illustrative plots of the ratios

# Energy range
x = np.array( ratio_WCC_tri[0] )

# Ratio of on-shell W / trident total
y1 = np.array( ratio_WCC_tri[1] )
y1err = np.array( ratio_WCC_tri[2] )

# Ratio of on-shell W / trident CC
y2 = np.array( ratio_WCC_triCC[1] )
y2err = np.array( ratio_WCC_triCC[2] )

# Ratio of trident CC / trident total
y3 = np.array( ratio_triCC_tri[1] )
y3err = np.array( ratio_triCC_tri[2] )

fig, (ax1) = plt.subplots()

# ax1.set_title('test')
ax1.errorbar(x, y1, yerr=y1err)

# ax1.set_title('test')
ax1.errorbar(x, y2, yerr=y2err)

# ax1.set_title('test')
ax1.errorbar(x, y3, yerr=y3err)

ax1.legend(['W / trident (all)','W / trident (CC)', 'trident (CC) / trident (all)'])

# ax1.fill_between(x,minp,maxp,alpha=0.5,label='LO')
ax1.set_ylim(0.85,1.05)
ax1.set_xlim(1e-2,2e4)
ax1.set_xlabel('Ecms [GeV]')
ax1.set_ylabel('ratio of cross-sections')
#r'$\sigma$ ratio')

plt.xscale("log")
plt.savefig('Cross_Section_ratios.pdf')
