#!/usr/bin/env python

""" plots SB histograms of the EAGLE data

Usage: SB_distribution_hist.py [-h] [-v] [-p]

Options:
    -h, --help                          Show this screen.
    -v, --verbose                       Show extra information [default: False]
    -p, --plotchecks                    Plot checks to make sure working [default: False]

Examples: python SB_distribution_hist.py -v -p 

"""


import os
import numpy as np
import docopt

import eagle_constants_and_units as c
import cosmo_utils as csu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
from astropy import constants as const
from astropy import units as u

import get_halpha_SB

def testplot():
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

pixscale =  {'50Mpc': 0.237/1000.*(1.+0.0115), '100Mpc': 0.477/1000.*(1.+0.0235),'200Mpc': 0.928/1000.*(1.+0.047) , '500Mpc': 2.178/1000.*(1.+0.12)} ### Mpc / arcsec (comoving)
x_center = 47.5
y_center = 12.
x_angFOV = 3.*60.*60. # " 
y_angFOV = 2.*60.*60. # "  
x_FOV = {distance: pixscale[distance]*x_angFOV for distance in ['50Mpc','100Mpc','200Mpc','500Mpc']}  # cMpc
y_FOV = {distance: pixscale[distance]*y_angFOV for distance in ['50Mpc','100Mpc','200Mpc','500Mpc']}  # cMpc

def loaddata():
    fnames_100 = {'50Mpc':'data_50Mpc_100arcsec.npz','100Mpc':'data_100Mpc_100arcsec.npz','200Mpc':'data_200Mpc_100arcsec.npz','500Mpc':'data_500Mpc_100arcsec.npz'}  ## two level dictionary when want different resolutions?
    resolution=100. #arcsec
    # prep data locations
    data_tuple_array = []
    
    for distance,ind in zip(['50Mpc','100Mpc','200Mpc','500Mpc'],[0,1,2,3]):
        fname = fnames_100[distance]
        if os.path.isfile(fname):
            print("data tuple exists, loading %s now..."%fname)
            data_tuple_array.append(np.load(fname)['arr_0'])
        else:
            print("data not saved, loading from total res, 5Mpc slice, file now...")
            fname = 'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_total.npz'
            if os.path.isfile(fname):
                print("data exists, loading %s now..."%fname)
                sl = [slice(None,None,None), slice(None,None,None)]
                data = (np.load(fname)['arr_0'])[sl]
            else:
                print("data not saved, loading from original files now...")
                data = loaddata()
                np.savez(fname,data)
                
            data_tuple = changeres(distance,100,data)
            data_tuple_array.append(data_tuple)
            np.savez(fname,data_tuple)
            
    data_50_100_tuple=data_tuple_array[0]; data_100_100_tuple=data_tuple_array[1]; 
    data_200_100_tuple=data_tuple_array[2]; data_500_100_tuple=data_tuple_array[3];
    
    data_dict = {'50Mpc':data_50_100_tuple,'100Mpc':data_100_100_tuple,'200Mpc':data_200_100_tuple,'500Mpc':data_500_100_tuple}

    return data_dict
    
def extractFOVs():
    """
    Input the EAGLE data and crops it to the Dragonfly FOV at a distance
    """
    data_FOV = []
    for distance in ['50Mpc','100Mpc','200Mpc','500Mpc']:
        data_tuple = data_dict[distance]
        factor = data_tuple[2]; newsize = data_tuple[1]; data = data_tuple[0]; resolution = 100.
    
        xystarts = [x_center-x_FOV[distance]/2.,y_center-y_FOV[distance]/2.]
        size     = [x_FOV[distance], y_FOV[distance]]
    
        x1 = ((x_center-x_FOV[distance]/2.)/100.*(newsize/factor))
        x2 = ((x_center+x_FOV[distance]/2.)/100.*(newsize/factor))
        y1 = ((y_center-y_FOV[distance]/2.)/100.*(newsize/factor))
        y2 = ((y_center+y_FOV[distance]/2.)/100.*(newsize/factor))
        data_FOV.append(data[int(x1):int(x2),int(y1):int(y2)])
    
    data_50_100_FOV=data_FOV[0]; data_100_100_FOV=data_FOV[1]; data_200_100_FOV=data_FOV[2]; data_500_100_FOV=data_FOV[3]
    
    return data_50_100_FOV,data_100_100_FOV,data_200_100_FOV,data_500_100_FOV

def testplot():
    # plot them all together
    f, ((ax1, ax2, ax3, ax4)) = plt.subplots(1, 4, figsize = (15.5, 15.))

    for distance,data,axis in zip(['50Mpc','100Mpc','200Mpc','500Mpc'],[data_50_100_FOV,data_100_100_FOV,data_200_100_FOV,data_500_100_FOV],[ax1,ax2,ax3,ax4]):
        xystarts = [x_center-x_FOV[distance]/2.,y_center-y_FOV[distance]/2.]
        size     = [x_FOV[distance], y_FOV[distance]]
        get_halpha_SB.makemap(data,size,axis,xystarts = xystarts)
    
def plothists():
    f, ((ax1, ax2, ax3, ax4)) = plt.subplots(1, 4, figsize = (15.5, 15.))
    
    for distance,data,axis in zip(['50Mpc','100Mpc','200Mpc','500Mpc'],[data_50_100_FOV,data_100_100_FOV,data_200_100_FOV,data_500_100_FOV],[ax1,ax2,ax3,ax4]):
    
        hist, bin_edges = np.histogram(data, bins=50, range=[0.,5.],density=True)
        """
        From python documentation for numpy.histogram:
        density : bool, optional
        If False, the result will contain the number of samples in each bin. 
        If True, the result is the value of the probability density function at the bin, 
        normalized such that the integral over the range is 1. 
        Note that the sum of the histogram values will not be equal to 1 unless bins of unity width are chosen; 
        it is not a probability mass function.
        Overrides the normed keyword if given.
        """
        axis.hist(np.ravel(data), bins=50, range=[0.,5.], density=True)

#-------------------------------------- BODY OF PROGRAM STARTS HERE ---------------------------------------------#

if __name__ == "__main__":
    
    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    plotchecks      = arguments['--plotchecks']

    machine='coho'

    if verbose:
        print('Loading data...')
    data_dict=loaddata()
    if verbose:
        print('selecting out FOV data...')
    data_50_100_FOV,data_100_100_FOV,data_200_100_FOV,data_500_100_FOV=extractFOVs()
    
    if plotchecks:
        testplot()
        plt.tight_layout()
        plt.show()

    plothists()

