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

	output = [ [],[],[],[],[] ]
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
		# temp_line = []
		# for idata in range(0,len(line_data)):
		for idata in range(0,len(line_data)):
			output[idata].append( float(line_data[idata]) )
		# output.append( temp_line )

	return output


file = 'SigmaIncl_Ecms_channel1_virt_s3.txt'

test = read_data(file)

print( test )

# Energy range
x = np.array( test[0] )
# Ratio of on-shell W / trident total
factor = 1./3.
y_ratio = []
y_ratio_err = []
y1 = np.array( test[4] )
y1err = np.array( test[2] )
y1err_plot = y1err/np.array( test[1] )
fig, (ax1) = plt.subplots()
# ax1.set_title('test')
ax1.errorbar(x, y1, yerr=y1err_plot)
ax1.legend(['nu_e + nu_e > nu_e nu_e'])

# ax1.fill_between(x,minp,maxp,alpha=0.5,label='LO')
ax1.set_ylim(0.0,1.15)
ax1.set_xlim(1,1e8)
ax1.set_xlabel('Ecms [GeV]')
ax1.set_ylabel('ratio of sigma(Enu) NLO EW / LO')
plt.xscale("log")
plt.grid()
plt.savefig('test.pdf')




# Calculate the ratio between these channels
# ratio_incl_over_quarks = compute_ratio( sigma_mW[0], sigma_mW[1] )

# plot_mw = True
# # Construct a plot
# if( plot_mw ):
#