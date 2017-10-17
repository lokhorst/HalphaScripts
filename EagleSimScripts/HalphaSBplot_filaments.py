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

#-------------- pick your distance and size of binning along one side --------------#
binsize = 100.  ### arcsec
distance = '50Mpc'  ### '50Mpc' '100Mpc' '200Mpc' '500Mpc'
boxnum = '1' ### which filament (there are 3)
#-----------------------------------------------------------------------------------#


def plotfilament(onlyyellow=False):
    clabel = r'log photons/cm$^2$/s/sr'
    Vmin = None
    Vmax= None
    fig = plt.figure(figsize = (7.5, 8.))
    ax = plt.subplot(121)
    fontsize=13
    #ax.set_xlabel(r'X [cMpc]',fontsize=fontsize)
    #ax.set_ylabel(r'Y [cMpc]',fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    colmap = 'viridis' #'afmhot'
    ax.patch.set_facecolor(cm.get_cmap(colmap)(0.)) # sets background color to lowest color map value

    ## If you only want to plot the SB greater than 1 photon/s/cm^2/arcsec^2 then do the following
    if onlyyellow:
        SBonlyyellow = SBdata_5
        SBonlyyellow[SBdata_5<0.] = -3.
        img = ax.imshow(SBonlyyellow.T,origin='lower', cmap=cm.get_cmap(colmap), vmin = Vmin, vmax=Vmax,interpolation='nearest')
        levels = [0,1,2]
        colours = ['yellow','cyan','purple']
    else:
        img = ax.imshow(SBdata_5.T,origin='lower', cmap=cm.get_cmap(colmap), vmin = Vmin, vmax=Vmax,interpolation='nearest')
        levels = np.array([-2,-1,0,1,2,3])
        colours = ('red','orange','yellow','cyan','purple','pink')

    cmap = cm.PRGn
    plt.contour(SBdata_5.T,levels,colors=colours)#,cmap=cm.get_cmap(cmap, len(levels) - 1),)

    div = axgrid.make_axes_locatable(ax)
    cax = div.append_axes("bottom",size="15%",pad=0.1)
    cbar = plt.colorbar(img, cax=cax,orientation='horizontal')
    cbar.solids.set_edgecolor("face")
    cbar.ax.set_xlabel(r'%s' % (clabel), fontsize=fontsize)
    #cbar.ax.set_ylabel(r'%s' % (clabel), fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    

pixscale =  {'50Mpc': 0.237/1000.*(1.+0.0115), '100Mpc': 0.477/1000.*(1.+0.0235), '200Mpc': 0.928/1000.*(1.+0.047) , '500Mpc': 2.178/1000.*(1.+0.12)} ### Mpc / arcsec (comoving)
simpixsize = 100./32000. ### Mpc / pixel
factor = pixscale[distance]*100./simpixsize
sl = [slice(None,None,None), slice(None,None,None)]

# Simulation snapnum 28 (z = 0), xy box size: 100Mpc, z slice width: 5Mpc,
files_SF_28 = ['/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5__fromSFR.npz']

files_noSF_28 = ['/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5_noSFR.npz']
                 
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

# Define boxes around the filaments (snapnum 28)
xbox_3 = np.array([53,53,56,56])*3200./100.
ybox_3 = (np.array([9.2,10,8.5,7.7])-0.2)*3200./100.
xbox_2 = np.array([44.5,44.5,46,46])*3200./100.
#    xbox_2 = np.array([43,43,46,46])*3200./100.
ybox_2 = (np.array([(7.9+6.9)/2.,(8.05+7.05)/2.,7.05,6.9])-0.05+0.2)*3200./100.
#    ybox_2 = (np.array([7.9,8.05,7.05,6.9])-0.05+0.2)*3200./100.
##    ybox_2 = (np.array([7.8,8.1,7.1,6.8])-0.05+0.2)*3200./100.
xbox_1 = (np.array([47.4,46.2,46.9,48.1])+0.5)*3200./100.
ybox_1 = np.array([10.5,14,14,10.5])*3200./100.
xboxes = {'1':xbox_1,'2':xbox_2,'3':xbox_3}
yboxes = {'1':ybox_1,'2':ybox_2,'3':ybox_3}

# Pull out the SB data from the boxes defined above
xfull, yfull= get_halpha_SB.indices_region(xboxes[boxnum].astype(int),yboxes[boxnum].astype(int)) 
SBdata_5 = np.zeros(xfull.shape)
for i in range(yfull.shape[0]):
    for j in range(yfull.shape[1]):
        SBdata_5[i,j]  = data_5[xfull[i,j],yfull[i,j]]
SBdata_average = np.log10(np.mean(10**SBdata_5))
SBdata_median  = np.median(SBdata_5)

# Plot the original data around the region we pulled out to do a cross-check
#fig = plt.figure(figsize = (16.5, 15.))
ax1 = plt.subplot(122)
#factor = 1 ## only comment this out if you are using the  reduced resolution
factor = 2.
get_halpha_SB.makemap(data_5[(xystarts[0]/100.*3200./factor):((xystarts[0]+size)/100.*3200./factor),(xystarts[1]/100.*3200./factor):((xystarts[1]+size)/100.*3200./factor)],size,ax1,xystarts = xystarts)
ax1.plot(np.append(xbox*100./3200.*factor,xbox[0]*100./3200.*factor),np.append(ybox*100./3200.*factor,ybox[0]*100./3200.*factor),color='r')
plt.show()

