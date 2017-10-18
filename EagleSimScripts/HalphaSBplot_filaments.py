"""
Script streamlining the process of plotting the filaments with surface brightness contours and binned correctly.

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

def addnoise(data)
### DOESN'T WORK YET ###
    exptime = 1. # just one second since we are still in detections per second mode

    # Dragonfly info
    area_lens = np.pi*(14.3/2)**2 * 48.               # cm^2, 48 * 14.3 cm diameter lenses
    pix_size = 2.8                                    # arcsec
    ang_size_pixel  = (pix_size * (1./206265.))**2    # rad^2, the pixel size of the CCD

    tau_l = 0.85  # transmittance of the Dragonfly lens
    QE_new = 0.70  # quantum efficiency of the CMOS detector
    QE_old = 0.48     # quantum efficiency of the CCDs
    tau_f = 1.    # transmittance of the Halpha filter -- assumed for now
    B = HalphaSBplot_addnoise.getBackground(656.3,659.3) # *u.photon/u.second/u.arcsec**2/u.m**2  ****already multiplied by the bandwidth***
    D = 0.04  # *u.photon/u.second                             # dark current (electrons / s) 
    R_squared_old = 10.**2 # * u.photon                           # read noise (electrons)
    R_squared_new = 1.**2 # * u.photon                           # read noise (electrons)

    binpix_size = 100. # arcsec
    numpixel = round((binpix_size/pix_size)**2)
    detsignal = np.log10(10**totsignal * QE_old * tau_l * tau_f * area_lens * ang_size_pixel * numpixel)

    B_sky = B * QE_old * tau_l * tau_f * area_lens*(1/100)**2 * pix_size**2
    sigma = np.log10(np.sqrt(10**detsignal + B_sky*exptime*numpixel + D*exptime*numpixel + R_squared_old*numpixel))

def plotfilament(SBdata,ax,onlyyellow=False):
    # setting up the plot
    clabel = r'log photons/cm$^2$/s/sr'
    Vmin = None
    Vmax= None
    #fig = plt.figure(figsize = (7.5, 8.))
    #ax = plt.subplot(121)
    fontsize=13
    #ax.set_xlabel(r'X [cMpc]',fontsize=fontsize)
    #ax.set_ylabel(r'Y [cMpc]',fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    colmap = 'viridis' #'afmhot'
    ax.patch.set_facecolor(cm.get_cmap(colmap)(0.)) # sets background color to lowest color map value

    ## If you only want to plot the SB greater than 1 photon/s/cm^2/arcsec^2 then do the following
    #if onlyyellow:
    #    SBonlyyellow = SBdata
    #    SBonlyyellow[SBdata<0.] = -3.
    #    img = ax.imshow(SBonlyyellow.T,origin='lower', cmap=cm.get_cmap(colmap), vmin = Vmin, vmax=Vmax,interpolation='nearest')
    #    levels = [0,1,2]
    #    colours = ['yellow','cyan','purple']
    #else:
    img = ax.imshow(SBdata.T,origin='lower', cmap=cm.get_cmap(colmap), vmin = Vmin, vmax=Vmax,interpolation='nearest')
    levels = np.array([-2,-1,0,1,2,3])
    colours = ('red','orange','yellow','cyan','purple','pink')
    levels = np.array([-2,-1.5,-1,-0.5,0,0.3,1,1.5,2,2.5,3])
    colours = ('red','black','orange','black','yellow','black','cyan','black','purple','black','pink')
    
    # plot contours
    cmap = cm.PRGn
    ax.contour(SBdata.T,levels,colors=colours)#,cmap=cm.get_cmap(cmap, len(levels) - 1),)

    div = axgrid.make_axes_locatable(ax)
    cax = div.append_axes("bottom",size="15%",pad=0.1)
    cbar = plt.colorbar(img, cax=cax,orientation='horizontal')
    cbar.solids.set_edgecolor("face")
    cbar.ax.set_xlabel(r'%s' % (clabel), fontsize=fontsize)
    #cbar.ax.set_ylabel(r'%s' % (clabel), fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    
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
        print("the new size of the data array is %s."%newsize)
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


#-------------- pick your distance and desired resolution --------------#
resolution = 100.  ### arcsec
distance = '50Mpc'  ### '50Mpc' '100Mpc' '200Mpc' '500Mpc'
boxnum = '1' ### which filament (there are 3)
factor = 1
machine='coho'
#-----------------------------------------------------------------------------------#

data_5 = loaddata() # load in data at full resolution

# pull out the pixel limits for boxes that surround the three filaments
xboxes, yboxes = defineboxes(data_5)

# takes in pixel limits that bound a box and create arrays of x and y pixel values to pick out SB 
xfull, yfull= get_halpha_SB.indices_region(xboxes[boxnum].astype(int),yboxes[boxnum].astype(int)) 

# use pixel arrays to extract SB data in a box from the data array
SBdata_5 = extractdata(xfull,yfull,data_5)

SBdata_average = np.log10(np.mean(10**SBdata_5))
SBdata_median  = np.median(SBdata_5)

# plot the filament with contours
fig = plt.figure(figsize = (7.5, 8.))
ax = plt.subplot(121)
plotfilament(SBdata_5,ax)

# repeat with different resolution
data_50Mpc_100arcsec, newsize, factor = changeres(distance,resolution,data_5) # change data to required resolution at selected distance
xboxes, yboxes = defineboxes(data_50Mpc_100arcsec)
xfull, yfull= get_halpha_SB.indices_region(xboxes[boxnum].astype(int),yboxes[boxnum].astype(int)) 

SBdata_50Mpc_100arcsec = extractdata(xfull,yfull,data_50Mpc_100arcsec)
fig = plt.figure(figsize = (7.5, 8.))
ax = plt.subplot(121)
plotfilament(SBdata_50Mpc_100arcsec,ax)

# repeat with different resolution
resolution = 500. ### arcsec
data_50Mpc_500arcsec, newsize, factor = changeres(distance,resolution,data_5) # change data to required resolution at selected distance
xboxes, yboxes = defineboxes(data_50Mpc_500arcsec)
xfull, yfull= get_halpha_SB.indices_region(xboxes[boxnum].astype(int),yboxes[boxnum].astype(int)) 

SBdata_50Mpc_500arcsec = extractdata(xfull,yfull,data_50Mpc_500arcsec)
fig = plt.figure(figsize = (7.5, 8.))
ax = plt.subplot(121)
print('SBdata_50Mpc away, 500arcsec per pix, %s Mpc per pix'%(newsize/32000.*100./SBdata_50Mpc_500arcsec.shape[0]))
plotfilament(SBdata_50Mpc_500arcsec,ax)
plt.title('SBdata_50Mpc_500arcsec %s Mpc per pixel'%(newsize/32000.*100./SBdata_50Mpc_500arcsec.shape[0]))

plt.show()


# Plot the original data around the region we pulled out to do a cross-check
#fig = plt.figure(figsize = (16.5, 15.))
ax1 = plt.subplot(122)
factor = 1. ## only comment this out if you are using the  reduced resolution

get_halpha_SB.makemap(data_5[(xystarts[0]/100.*3200./factor):((xystarts[0]+size)/100.*3200./factor),(xystarts[1]/100.*3200./factor):((xystarts[1]+size)/100.*3200./factor)],size,ax1,xystarts = xystarts)
ax1.plot(np.append(xbox*100./3200.*factor,xbox[0]*100./3200.*factor),np.append(ybox*100./3200.*factor,ybox[0]*100./3200.*factor),color='r')
plt.show()

