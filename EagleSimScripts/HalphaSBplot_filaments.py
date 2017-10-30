#!/usr/bin/env python

""" Plots a filament from EAGLE simulation (contours overlaid or noise added for mock observation)

Usage: HalphaSBplot_filaments.py [-h] [-s] [-v] [--debug] [-c] [-r RESOLUTION] [-t EXPTIME] [--mockobs]

Options:
    -h, --help                          Show this screen.
    -v, --verbose                       Show extra information [default: False]  
    -s, --save                          Save intermediate array products [default: False]
    --debug                             Show extra extra information [default: False]  

    -c, --cmos                          Use specs for new cameras in mock observation
    --mockobs                           Plot a mock observation (rather than just the EAGLE simulation). 
    -r RESOLUTION, --res RESOLUTION     The desired resolution of the filament image in arcsec. [default: 500]
    -t EXPTIME, --exptime EXPTIME       The desired exposure time to integrate the filament over in seconds. [default: 10**3*3600.]

Examples:

Longer description:
Script streamlining the process of plotting the EAGLE filaments (with surface brightness contours and binned correctly).
Takes inputs for distance to filament and desired angular resolution of filament, then bins EAGLE data array read in from 
npz files obtained from Nastasha Wijers at Leiden Observatory.  Plots filament with SB contours overlaid.
Adds noise to create mock Dragonfly observations.
"""

import numpy as np
import eagle_constants_and_units as c
import cosmo_utils as csu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
from astropy import constants as const
from astropy import units as u

import get_halpha_SB
import HalphaSBplot_addnoise

def getBackground(start,end,machine,plot=False):
    # Returns the total background flux in the wavlength interval supplied i.e. returns (flux)*(wavlength interval) 
    wavelength = []
    flux = []
    
    if machine=='chinook':
        geminiloc='/Users/lokhorst/Documents/Eagle/Gemini_skybackground.dat'
    elif machine=='coho':
        geminiloc='/Users/deblokhorst/Documents/Dragonfly/HalphaScripts/Gemini_skybackground.dat'
    
    with open(geminiloc,'r') as f:  #wavelength in nm, flux in phot/s/nm/arcsec^2/m^2
        for line in f:
            if line[0]!='#' and len(line)>5:
                tmp = line.split()
                wavelength.append(tmp[0])
                flux.append(tmp[1])
                
    wavelength = np.array(wavelength,'d')
    flux = np.array(flux,'d')
    
    start_ind = (np.abs(wavelength-start)).argmin()
    end_ind   = (np.abs(wavelength-end)).argmin()
    
    # if spacings are not even, need to add element by element
    total=0
    for index in np.arange(start_ind,end_ind):
      #  print index,index+1
      #  print total
        total = total+(flux[index]*(wavelength[index+1]-wavelength[index]))
        
    # if spacings are even, can just take the average of the flux array and times it by the total bandwidth
    np.mean(flux[start_ind:end_ind])*(wavelength[end_ind]-wavelength[start_ind])
    
   # print('start index and end index: %s and %s'%(start_ind,end_ind))
   # print(wavelength[start_ind:end_ind]-wavelength[start_ind+1:end_ind+1])
    if plot==True:
        plt.plot(wavelength[start_ind-100:end_ind+100],flux[start_ind-100:end_ind+100])
        top = max(flux[start_ind-100:end_ind+100])
        plt.plot([start,start,end,end,start],[0,top,top,0,0])
        plt.show()
        
    return total
    
    
def addnoise(data,resolution,exptime=10**3*3600.,CMOS=False):
    # Dragonfly info
    area_lens = np.pi*(14.3/2)**2 * 48. *10.                # cm^2, 48 * 14.3 cm diameter lenses
    pix_size = 2.8                                      # arcsec
    ang_size_pixel  = (pix_size * (1./206265.))**2      # rad^2, the pixel size of the CCD
    tau_l = 0.85                                        # transmittance of the Dragonfly lens
    tau_f = 1.                                          # transmittance of the Halpha filter -- assumed for now
    #B = getBackground(656.3,657.3,machine)              # *u.photon/u.second/u.arcsec**2/u.m**2  ****already multiplied by the bandwidth***
    B = 0.560633
    D = 0.04       # dark current (electrons / s) 
    
    if CMOS:
   #     print "Using new CMOS cameras..."
        QE = 0.70                                       # quantum efficiency of the CMOS detector
        R_squared = 2.**2                               # read noise (electrons)
    else:
   #     print "Using old cameras..."
        QE = 0.48                                       # quantum efficiency of the CCDs
        R_squared = 10.**2                              # read noise (electrons)
    
   # R_squared = 50.**2
    
    binpix_size = resolution # arcsec
    numpixel = round((binpix_size/pix_size)**2)
    #print "the number of pixels is %s"%numpixel
    
    
    ### total signal incident in exposure time ###
    totsignal = 10**data * exptime # ( photons / cm^2 /sr )
    ### total signal detected (accounting for system efficiency) ###
    detsignal = totsignal * QE * tau_l * tau_f * area_lens * ang_size_pixel * numpixel
    #print "the total object signal [electrons] detected ranges from: %s to %s"%(np.min(detsignal),np.max(detsignal))
    #print "an example of the object signal [electrons] is: %s"%detsignal[0]


    ### Background from stuff in space ###
    'background sky signal detected [B]=ph/s/arcsec^2/m^2, [B_sky]=ph/s (in a pixel)'
    B_sky = B * QE * tau_l * tau_f * area_lens*(1/100.)**2 * pix_size**2
    #print "the background in the bandwidth is: %s"%B
    #print "the background signal, B_sky [ph/s (in a pixel)], is: %s"%B_sky
    B_sky_inexptime = B_sky*exptime
    B_sky_total     = B_sky*exptime*numpixel    
    B_sky_array = np.zeros((data.shape[0],data.shape[1]))
    for x in range(data.shape[0]):
        for y in range(data.shape[1]):
            B_sky_array[x][y]=np.random.normal(0,np.sqrt(B_sky_total+detsignal[x][y])) 
#            B_sky_array[x][y]=np.random.poisson(B_sky_total)
    B_sky_array_total = B_sky_array
    #print "the mean total background signal, B_sky_total [electrons], is: %s"%B_sky_total
    #print "the total background noisy signal [electrons] ranges from: %s to %s"%(np.min(B_sky_array_total),np.max(B_sky_array_total))
    
    # Signal
    noiseadded_signal = detsignal + B_sky_total + B_sky_array
    
    ### ReadOutNoise ###
    numexposures = exptime/3600. # hour long exposures
    R_squared_array = np.zeros((data.shape[0],data.shape[1]))
    for x in range(data.shape[0]):
        for y in range(data.shape[1]):
            R_squared_array[x][y]=np.mean(np.random.normal(np.sqrt(R_squared),np.sqrt(np.sqrt(B_sky)),int(numpixel)))**2   
    R_squared_total = R_squared * round(numexposures)
    R_squared_total_array = R_squared_array * round(numexposures)
    #print "the R_squared value is: %s, so in %s exposures [per pixel], will have R_squared of: %s, %s"%(R_squared,numexposures,R_squared_total,R_squared_total_array[0])
    #print "the total R_squared value [electrons] multiplying by numpix read out is: %s, %s"%((R_squared_total*numpixel),(R_squared_total_array[0]*numpixel))
    
    ### DarkCurrent ###
    noise_from_detector = 0.0
    #D_total = D*exptime*numpixel
    #D_array = np.zeros((data.shape[0],data.shape[1]))
    #for x in range(data.shape[0]):
    #    for y in range(data.shape[1]):
    #        D_array[x][y]=np.random.normal(D_total,np.sqrt(D_total)) 
    #D_array_total = D_array
    #print "the total dark current [electrons] is: %s , %s"%(D_total, D_array_total[0])

    #noise_from_detector = D_array_total + R_squared_total_array*numpixel
    #print "an example total noise (not squarerooted) is: %s"%(detsignal + B_sky_array_total + D_array_total + R_squared_total_array*numpixel)[0]
    #print "an example total noise (squarerooted) is: %s"%sigma[0]
    
    
    # Now add noise from the detector
    
    noiseadded_signal = noiseadded_signal + noise_from_detector
    
    return noiseadded_signal

def plotfilament(SBdata,ax,colmap='viridis',onlyyellow=False,contours=False,colorbar=False,mockobs=False,labelaxes=False,label=''):
    # setting up the plot
    if mockobs:
        clabel = r'log signal (photons)'
    else:
        clabel = r'log photons/cm$^2$/s/sr'
    Vmin = None
    Vmax= None
    #fig = plt.figure(figsize = (7.5, 8.))
    #ax = plt.subplot(121)
    fontsize=13

    if labelaxes:
        ax.set_xlabel(r'X [cMpc]',fontsize=fontsize)
        ax.set_ylabel(r'Y [cMpc]',fontsize=fontsize)
    
        ax.tick_params(labelsize=fontsize) #,top=True,labeltop=True)
        ax.xaxis.set_label_position('top') 
        ax.xaxis.tick_top()
    else:
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        
    #colmap = 'viridis'#'gist_gray'#'plasma'#'viridis' #'afmhot'
    ax.patch.set_facecolor(cm.get_cmap(colmap)(0.)) # sets background color to lowest color map value

    
    ## If you only want to plot the SB greater than 1 photon/s/cm^2/arcsec^2 then do the following
    if onlyyellow:
        SBonlyyellow = SBdata
        SBonlyyellow[SBdata<0.] = -3.
        img = ax.imshow(SBonlyyellow.T,origin='lower', cmap=cm.get_cmap(colmap), vmin = Vmin, vmax=Vmax,interpolation='nearest')
        levels = [0,1,2]
        colours = ['yellow','cyan','purple']
    else:
        img = ax.imshow(SBdata.T,origin='lower',extent=(0,3.7,0,0.7), cmap=cm.get_cmap(colmap), vmin = Vmin, vmax=Vmax,interpolation='nearest')
        levels = np.array([-2,-1,0,1,2,3])
        colours = ('red','orange','yellow','cyan','purple','pink')
        #levels = np.array([-2,-1.5,-1,-0.5,0,0.3,1,1.5,2,2.5,3])
        #colours = ('red','black','orange','black','yellow','black','cyan','black','purple','black','pink')
    
    # plot contours
    cmap = cm.PRGn
    if contours:
        ax.contour(SBdata.T,levels,colors=colours)#,cmap=cm.get_cmap(cmap, len(levels) - 1),)

    div = axgrid.make_axes_locatable(ax)
    # plot colorbar
    if colorbar:
        cax = div.append_axes("bottom",size="15%",pad=0.1)
        cbar = plt.colorbar(img, cax=cax,orientation='horizontal')
        cbar.solids.set_edgecolor("face")
        cbar.ax.set_xlabel(r'%s' % (clabel), fontsize=fontsize)
        cbar.ax.set_ylabel(r'%s' % (clabel), fontsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)
    
    font = {'family': 'serif',
        'color':  'yellow',
        'weight': 'bold',
        'size': 12,
        }
    
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
    
    # Top right
    # plt.text(3.6,0.58,label,fontdict=font,horizontalalignment='right',backgroundcolor='black')
    
    # Bottom middle
    plt.text(1.8,0.8,label,fontdict=font,horizontalalignment='center',backgroundcolor='white')

    plt.tight_layout()


def loaddata():
    sl = [slice(None,None,None), slice(None,None,None)]
    if machine=='chinook':
        homedir='/Users/lokhorst/Eagle/'
    elif machine=='coho':
        homedir='/Users/deblokhorst/eagle/SlicesFromNastasha/'

    # Simulation snapnum 28 (z = 0), xy box size: 100Mpc, z slice width: 5Mpc,
    files_SF_28 = [homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5__fromSFR.npz',
                   homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5__fromSFR.npz',
                   homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5__fromSFR.npz',
                   homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5__fromSFR.npz']

    files_noSF_28 = [homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_noSFR.npz',
                     homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5_noSFR.npz',
                     homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5_noSFR.npz',
                     homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5_noSFR.npz']
                 
    # Load a 5Mpc slice of data
    print('data1 ('+files_noSF_28[0]+')...')
    data1 = (np.load(files_noSF_28[0])['arr_0'])[sl]
    data1 = get_halpha_SB.imreduce(data1, round(factor), log=True, method = 'average')
    print('data11 ('+files_SF_28[0]+')...')
    data11 = (np.load(files_SF_28[0])['arr_0'])[sl]
    data11 = get_halpha_SB.imreduce(data11, round(factor), log=True, method = 'average')
    print('5 Mpc slice...')
    data_5 = np.log10(10**data1+10**data11)
    print('delete data1, data11...')
    del data1
    del data11
    
    return data_5

def changeres(distance,resolution,data):
    pixscale =  {'50Mpc': 0.237/1000.*(1.+0.0115), '100Mpc': 0.477/1000.*(1.+0.0235),'200Mpc': 0.928/1000.*(1.+0.047) , '500Mpc': 2.178/1000.*(1.+0.12)} ### Mpc / arcsec (comoving)
    simpixsize = 100./32000. ### Mpc / pixel is resolution of raw data 
    factor = round(pixscale[distance]*resolution/simpixsize)
    size = 32000.
    # LATER determine the current resolution of the data. FOR NOW assume current resolution is 100 Mpc/ 32000 pixels ~ 3 kpc/pixel

    # If the factors are not integer multiples of 32000., I'll trim the data first and then imreduce it
    if 32000.%((factor)) != 0.:
        times_factor_fits_in = int(32000./factor)
        newsize = times_factor_fits_in * factor
        print("Before reducing resolution, the original data was trimmed to size %s."%newsize)
        datanew = data[0:int(newsize),0:int(newsize)]
    else:
        datanew = data
        newsize = size

    return get_halpha_SB.imreduce(datanew, round(factor), log=True, method = 'average'), newsize, factor

def defineboxes(data,size=100.):
    # size in Mpc = total box size of data
    pixlength = float(data.shape[0])
    
    # Define boxes around the filaments (snapnum 28)
    xbox_3 = np.array([53,53,56,56])*pixlength/size
    ybox_3 = (np.array([9.2,10,8.5,7.7])-0.2)*pixlength/size
    xbox_2 = np.array([44.5,44.5,46,46])*pixlength/size
    #    xbox_2 = np.array([43,43,46,46])*pixlength/size
    ybox_2 = (np.array([(7.9+6.9)/2.,(8.05+7.05)/2.,7.05,6.9])-0.05+0.2)*pixlength/size
    #    ybox_2 = (np.array([7.9,8.05,7.05,6.9])-0.05+0.2)*pixlength/size
    ##    ybox_2 = (np.array([7.8,8.1,7.1,6.8])-0.05+0.2)*pixlength/size
    xbox_1 = (np.array([47.4,46.2,46.9,48.1])+0.5)*pixlength/size
    ybox_1 = np.array([10.5,14,14,10.5])*pixlength/size
    xboxes = {'1':xbox_1,'2':xbox_2,'3':xbox_3}
    yboxes = {'1':ybox_1,'2':ybox_2,'3':ybox_3}
    
    return xboxes, yboxes

def extractdata(xfull,yfull,data):
    SBdata = np.zeros(xfull.shape)
    for i in range(yfull.shape[0]):
        for j in range(yfull.shape[1]):
                SBdata[i,j]  = data[xfull[i,j],yfull[i,j]]
    return SBdata

def getSBatfilament(data,resolution,distance):
### DOESN'T WORK YET ###
    datares, newsize, factor = changeres(distance,resolution,data) # change data to required resolution at selected distance
    xboxes, yboxes = defineboxes(datares)
    xfull, yfull= get_halpha_SB.indices_region(xboxes[boxnum].astype(int),yboxes[boxnum].astype(int)) 
    SBdata = extractdata(xfull,yfull,datares)
    return SBdata

#-------------------------------------- BODY OF PROGRAM STARTS HERE ---------------------------------------------#

if __name__ == "__main__":
    
    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    debugging       = arguments['--debug']
    saveloc         = arguments['--saveloc']

    resolution      = arguments['--res']
    exptime         = arguments['--exptime']
    mockobs         = arguments['--mockobs']
    cmos            = arguments['--cmos']
    
    if verbose:
        print arguments
    
    distance = '50Mpc'  ### '50Mpc' '100Mpc' '200Mpc' '500Mpc'
    boxnum = '1' ### which filament (there are 3)
    factor = 1
    machine='chinook'

    # load a slice of data
    data_5 = loaddata() # load in data at full resolution



    if mockobs:
        SBdata = getSBatfilament(data_5,resolution,distance)
        if saveloc:
            np.savez(saveloc+'/SBdata_%sarcsec_res.npz',SBdata)
        skynoise,mockobs = addnoise(SBdata,resolution,exptime=exptime,CMOS=cmos)
        # make the data look better (scale)
        mockobs_sub = mockobs - (int(np.min(mockobs)*100)/100.)
        # plot the data
        fig = plt.figure(figsize = (9.5, 5.))
        ax1 = plt.subplot(111)
        plotfilament(mockobs_sub**0.2,ax1,contours=False,mockobs=True,colmap='gist_gray')
        plt.savefig('mockobs_res%sas_exptime%shr.png'%(resolution,round(exptime/3600.)))
        
        
"""
        
    fig = plt.figure(figsize = (10.5, 5.))
    ax1 = plt.subplot(111)
    plotfilament(SBdata_full,ax1,contours=False,labelaxes=True)
    plt.show()

    fig = plt.figure(figsize = (7.5, 8.))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)

    plotfilament(SBdata_full,ax1,contours=False,labelaxes=True)
    plotfilament(SBdata_exp1_sub**0.2,ax2,colmap='gist_gray',contours=False,mockobs=True)
    plotfilament(SBdata_100_withnoise_sub**0.2,ax3,colmap='gist_gray',contours=False,mockobs=True)
    plt.show()
    
        
        
    resolution = 500.
    SBdata_500arcsec = getSBatfilament(data_5,resolution,distance)
    noise,SBdata_500arcsec_200hr = addnoise(SBdata_500arcsec,resolution,exptime=10**2*3600.,CMOS=True)
    SBdata_500arcsec_200hr_sub = SBdata_500arcsec_200hr-5.248
    fig = plt.figure(figsize = (9.5, 10.))
    ax1 = plt.subplot(111)
    plotfilament(SBdata_500arcsec_200hr_sub**0.2,ax1,contours=False,mockobs=True,colmap='gist_gray')
        

    #----------------------------------------- Add noise to a filament and then plot --------------------------#
    # first try adding noise to the SBdata in a filament and then plotting
    resolution = 500. ### arcsec
    data_50Mpc_500arcsec, newsize, factor = changeres(distance,resolution,data_5) # change data to required resolution at selected distance
    
    xboxes, yboxes = defineboxes(data_50Mpc_500arcsec)
    xfull, yfull= get_halpha_SB.indices_region(xboxes[boxnum].astype(int),yboxes[boxnum].astype(int)) 
    SBdata_50Mpc_500arcsec = extractdata(xfull,yfull,data_50Mpc_500arcsec)
    SBdata_50Mpc_500arcsec_withnoise = addnoise(SBdata_50Mpc_500arcsec,resolution,exptime=10**4*3600.,CMOS=True)
    fig = plt.figure(figsize = (7.5, 8.))
    ax = plt.subplot(121)
    print('SBdata_50Mpc away, 500arcsec per pix, %s Mpc per pix'%(newsize/32000.*100./SBdata_50Mpc_500arcsec.shape[0]))
    plotfilament(SBdata_50Mpc_500arcsec_withnoise,ax,contours=False,mockobs=True)
    ax2 = plt.subplot(122)
    plotfilament(SBdata_50Mpc_500arcsec,ax2,contours=False)
    plt.show()
    
    # Same distance, same exposure time, three different resolutions
    
    def plot_diffres():
        fig = plt.figure(figsize = (7.5, 8.))
        SBdata_100 = getSBatfilament(data_5,100,distance)
        SBdata_100_withnoise = addnoise(SBdata_100,100,exptime=10**4*3600.,CMOS=True)
        ax1 = plt.subplot(322)
        plotfilament(SBdata_100_withnoise,ax1,contours=False,mockobs=True)
        ax2 = plt.subplot(321)
        plotfilament(SBdata_100,ax2,contours=False)
    
        SBdata_500 = getSBatfilament(data_5,500,distance)
        SBdata_500_withnoise = addnoise(SBdata_500,500,exptime=10**4*3600.,CMOS=True)
        ax3 = plt.subplot(324)
        plotfilament(SBdata_500_withnoise,ax3,contours=False,mockobs=True)
        ax4 = plt.subplot(323)
        plotfilament(SBdata_500,ax4,contours=False)
    
        SBdata_1000 = getSBatfilament(data_5,1000,distance)
        SBdata_1000_withnoise = addnoise(SBdata_1000,1000,exptime=10**4*3600.,CMOS=True)
        ax5 = plt.subplot(326)
        plotfilament(SBdata_1000_withnoise,ax5,contours=False,mockobs=True)
        ax6 = plt.subplot(325)
        plotfilament(SBdata_1000,ax6,contours=False)
    
        plt.show()
    
    # Same distance, same resolution, different exposure times
    
    def plot_diffexptime(resolution,distance):
        SBdata = getSBatfilament(data_5,resolution,distance)
        noise,SBdata_exp0 = addnoise(SBdata,resolution,exptime=10**2*3600.,CMOS=True)    
        noise,SBdata_exp1 = addnoise(SBdata,resolution,exptime=10**3*3600.,CMOS=True)
        noise,SBdata_exp2 = addnoise(SBdata,resolution,exptime=10**4*3600.,CMOS=True)
        noise,SBdata_exp3 = addnoise(SBdata,resolution,exptime=10**5*3600.,CMOS=True)

        noise,SBdata_exp0 = addnoise(SBdata,resolution,exptime=10**2*3600.,CMOS=False)    

        
        fig = plt.figure(figsize = (7.5, 8.))
        ax1 = plt.subplot(221)
        plotfilament(SBdata,ax1)
        ax0 = plt.subplot(223)
        plotfilament(SBdata_exp0,ax0,contours=False,mockobs=True)
        #ax2 = plt.subplot(223)
        #plotfilament(SBdata_exp1,ax2,contours=False,mockobs=True)    
        ax3 = plt.subplot(224)
        plotfilament(SBdata_exp2,ax3,contours=False,mockobs=True)
        #ax4 = plt.subplot(224)
        #plotfilament(SBdata_exp3,ax4,contours=False,mockobs=True)
        plt.show()
        
        # 100 hours (adjust so pretty!)
        noise,SBdata_exp0 = addnoise(SBdata,resolution,exptime=10**2*3600.,CMOS=True)  
        SBdata_100hr_subtract = SBdata_exp0 - (int(np.min(SBdata_exp0)*100)/100.)
        fig = plt.figure(figsize = (9.5, 10.))
        ax1 = plt.subplot(111)
        plotfilament(SBdata_100hr_subtract**0.2,ax1,contours=False,mockobs=True,colmap='gist_gray')
        
        
    plot_diffexptime(500,distance)
    plot_diffexptime(100,distance)
    
    print('SBdata_50Mpc away, 500arcsec per pix, %s Mpc per pix'%round(31996.0/32000.*/SBdata.shape[0],2))
    
    #SBdata_1000hr_subtract=np.log10(10**SBdata_exp1-10**4.83)
    #SBdata_subtract = np.log10(10**SBdata_exp2 - 10**5.33)
    SBdata_1000hr_subtract=SBdata_exp1-5.89
    SBdata_subtract = SBdata_exp2 - 6.39
    
    
    fig = plt.figure(figsize = (9.5, 10.))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    plotfilament((SBdata_1000hr_subtract)**0.1,ax1,colmap='gist_gray',contours=False,mockobs=True)
    plotfilament((SBdata_subtract)**0.1,ax2,colmap='gist_gray',contours=False,mockobs=True) #10000 hr
    plt.show()
    
    
    fig = plt.figure(figsize = (9.5, 10.))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    plotfilament((SBdata_1000hr_subtract)**0.1,ax1,colmap='gist_gray',contours=False,mockobs=True)      # with flat noise added! (uniform)
    plotfilament((SBdata_1000hr_noisy_subtract)**0.1,ax2,colmap='gist_gray',contours=False,mockobs=True)               # with noisy noise added! (gaussian dist)
    plt.show()
    
    fig = plt.figure(figsize = (9.5, 10.))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    plotfilament(10**(SBdata_1000hr_subtract),ax1,colmap='gist_gray',contours=False,mockobs=True)
    plotfilament(10**(SBdata_subtract),ax2,colmap='gist_gray',contours=False,mockobs=True) #10000 hr
    plt.show()
    
    fig = plt.figure(figsize = (9.5, 10.))
    ax1 = plt.subplot(111)
    plotfilament((SBdata_1000hr_subtract)**0.2,ax1,colmap='gist_gray',contours=False,mockobs=True)
    plt.show()

    
    #----------------------------------------- Plot original data (check filament plot) ------------------------#
    # Plot the original data around the region we pulled out to do a cross-check
    #fig = plt.figure(figsize = (16.5, 15.))
    ax1 = plt.subplot(122)

    factor = 1. ## if you are plotting raw data (un-reduced in resolution)

    xbox_1 = (np.array([47.4,46.2,46.9,48.1])+0.5)*pixlength/size
    ybox_1 = np.array([10.5,14,14,10.5])*pixlength/size
    xystarts = [45.2,10.]    
    size=4.
    get_halpha_SB.makemap(data_5[(xystarts[0]/100.*3200./factor):((xystarts[0]+size)/100.*3200./factor),(xystarts[1]/100.*3200./factor):((xystarts[1]+size)/100.*3200./factor)],size,ax1,xystarts = xystarts)
    ax1.plot(np.append(xbox*100./3200.*factor,xbox[0]*100./3200.*factor),np.append(ybox*100./3200.*factor,ybox[0]*100./3200.*factor),color='r')
    plt.show()

"""