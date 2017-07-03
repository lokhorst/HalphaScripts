import numpy as np
import eagle_constants_and_units as c
import cosmo_utils as csu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
from astropy import constants as const
from astropy import units as u

import get_halpha_SB

file_snapnum28_noSF = '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_noSFR.npz' 
sl = [slice(None,None,None), slice(None,None,None)] 
data = (np.load(file_snapnum28_noSF)['arr_0'])[sl]
factor = 10
data = get_halpha_SB.imreduce(data, factor, log=True, method = 'average')

#logSBarray = np.arange(0,5,0.1)
#Numarray = np.zeroes(len(logSBarray))
#indexarray = np.arange(len(logSBarray))
#for logSB,index in zip(logSBarray,indexarray):
#    logNumarray[index] = len(logSBarray[logSBarray >= logSB and logSBarray < (logSB+0.1)])    
#plt.plot()

hist, bin_edges = np.histogram(data, bins=50, range=[0.,5.],density=True)
"""
From python documentation for numpy.histogram:
density : bool, optional
If False, the result will contain the number of samples in each bin. If True, the result is the value of the probability density function at the bin, 
normalized such that the integral over the range is 1. Note that the sum of the histogram values will not be equal to 1 unless bins of unity width are chosen; 
it is not a probability mass function.
Overrides the normed keyword if given.
"""
plt.hist(data, bins=50, range=[0.,5.], density=True)