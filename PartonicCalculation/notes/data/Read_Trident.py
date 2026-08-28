import os, sys, glob, math

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
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
		# Adjust to 3
		# for idata in range(0,len(line_data)):
		for idata in range(0,3):			
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




# New analysis based on energy scan around mW region
base = ''
file_structure = 'SigmaIncl_Ecms_channel'
seed = 's5'
# channels_mW = ['29','120']
channels_mW = ['29','113']
sigma_mW = []
for ichan in range(0,len(channels_mW)):
	file = base+file_structure+channels_mW[ichan]+'_'+seed+'.txt'
	sigma_mW.append( read_data(file) )

# Calculate the ratio between these channels
ratio_incl_over_quarks = compute_ratio( sigma_mW[0], sigma_mW[1] )

plot_mw = True
# Construct a plot
if( plot_mw ):
	# Energy range
	x = np.array( ratio_incl_over_quarks[0] )
	# Ratio of on-shell W / trident total
	factor = 1./9.
	y_ratio = []
	y_ratio_err = []
	for ientry in range(0,len(ratio_incl_over_quarks[1]) ):
		y_ratio.append( factor * ratio_incl_over_quarks[1][ientry] )
		y_ratio_err.append( factor * ratio_incl_over_quarks[2][ientry] )
	# Multipliy by 3
	print( y_ratio )
	y1 = np.array( y_ratio )
	y1err = np.array( y_ratio_err )
	fig, (ax1) = plt.subplots()
	# ax1.set_title('test')
	ax1.errorbar(x, y1, yerr=y1err)
	ax1.legend(['NWA / 2-to-3'])
	ax1.set_xscale('log')
	# ax1.fill_between(x,minp,maxp,alpha=0.5,label='LO')
	ax1.set_ylim(0.70,1.3)
	ax1.set_xlim(75,750)
	ax1.set_xlabel('E_{cms} [GeV]')
	ax1.set_ylabel('ratio of cross-sections')
	# plt.xscale("log")
	plt.grid()
	plt.savefig('Cross_Section_ratios_mw.pdf')


	# Energy range
	x = np.array(ratio_incl_over_quarks[0])

	# Ratio of on-shell W / trident total
	factor = 1. / 9.
	y_sig1 = []
	y_sig2 = []

	# FIX: Restored [0] to match the indexing inside the loop
	for ientry in range(0, len(sigma_mW[0][0])):
	    y_sig1.append(factor * sigma_mW[0][1][ientry] )
	    y_sig2.append(sigma_mW[1][1][ientry] )

	fig, ax2 = plt.subplots()

	# Plots with explicit labels
	ax2.plot(x, y_sig1, label='NWA', linestyle='--')
	ax2.plot(x, y_sig2, label='2-to-3', linestyle=':')

	# Automatically grabs labels from the plot commands
	ax2.legend()
	ax2.set_xscale('log')
	# Formatting axes
	ax2.set_ylim(0.0, 50)
	ax2.set_xlim(75, 750)
	ax2.set_xlabel(r'$E_{\mathrm{cms}}$ [GeV]')
	ax2.set_ylabel(r'$\sigma$ [pb]')

	ax2.text(110, 45, r'$\nu_{\mu} + \gamma \rightarrow \mu^{-} + ( W^{+}  \to \tau^{+} + \nu_{\tau} ) $', fontsize=12, color='black')

	# 1. Tell matplotlib to format minor ticks with normal scalar numbers (e.g., 200 instead of 10^2)
	ax2.xaxis.set_minor_formatter(ticker.ScalarFormatter())
	# 2. Prevent scientific notation formatting from taking over
	ax2.xaxis.set_major_formatter(ticker.ScalarFormatter()) 

	plt.grid()
	plt.savefig('Cross_Section_mw.pdf')

