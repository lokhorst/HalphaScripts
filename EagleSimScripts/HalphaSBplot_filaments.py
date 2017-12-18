#!/usr/bin/env python

""" Plots a filament from EAGLE simulation (contours overlaid or noise added for mock observation)

Usage: HalphaSBplot_filaments.py [-h] [-s] [-v] [-m] [-g] [--debug] [-c] [-r RESOLUTION] [-t EXPTIME]

Options:
    -h, --help                          Show this screen.
    -v, --verbose                       Show extra information [default: False]  
    -s, --saveloc                       Save intermediate array products [default: False]
    --debug                             Show extra extra information [default: False]  

    -p, --plotchecks                    Plot intermediate checks [default: False]

    -g, --maskgal                       Mask the galaxies before plotting [default: False]
    -c, --cmos                          Use specs for new cameras in mock observation [default: False] 
    -m, --mockobs                       Plot a mock observation (rather than just the EAGLE simulation). [default: False] 
    -r RESOLUTION, --res RESOLUTION     The desired resolution of the filament image in arcsec. [default: 500]
    -t EXPTIME, --exptime EXPTIME       The desired exposure time to integrate the filament over in seconds. [default: 10**3*3600.]

Examples:

Longer description:
Script streamlining the process of plotting the EAGLE filaments (with surface brightness contours and binned correctly).
Takes inputs for distance to filament and desired angular resolution of filament, then bins EAGLE data array read in from 
npz files obtained from Nastasha Wijers at Leiden Observatory.  Plots filament with SB contours overlaid.
Adds noise to create mock Dragonfly observations.
"""

import docopt
import os
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

import eagleSqlTools as sql

pixscale =  {'50Mpc': 0.237/1000.*(1.+0.0115), '100Mpc': 0.477/1000.*(1.+0.0235),'200Mpc': 0.928/1000.*(1.+0.047) , '500Mpc': 2.178/1000.*(1.+0.12)} ### Mpc / arcsec (comoving)

xbox_3 = np.array([53,53,56,56])
ybox_3 = (np.array([9.2,10,8.5,7.7])-0.2)
xbox_2 = np.array([44.5,44.5,46,46])
ybox_2 = (np.array([(7.9+6.9)/2.,(8.05+7.05)/2.,7.05,6.9])-0.05+0.2)
xbox_1 = (np.array([47.4,46.2,46.9,48.1])+0.5)
ybox_1 = np.array([10.5,14,14,10.5])
    
def maskgals2(xgal,ygal,rhgas,data,xystarts,size,distance):
    xsize = data.shape[0]
    ysize = data.shape[1]
    xlen  = size[0]
    ylen  = size[1]

    data_masked = np.array([[i for i in line] for line in data])
    for i in range(len(xgal)):
        xind = (xgal[i]-xystarts[0])*xsize/xlen
        yind = (ygal[i]-xystarts[1])*ysize/ylen
        if rhgas[i] > (pixscale[distance]*resolution*1000.):  # if the gas scale length is greater than the pixel size
            print("rhgas, %s, is greater than %s kpc"%(rhgas[i],round(pixscale[distance]*resolution*1000.)))
            xind = int(round(xind))
            yind = int(round(yind))
            if xind+2<data.shape[0] and yind+2<data.shape[1]:
                data_sub = np.mean([data[xind+2,yind+2],data[xind-2,yind-2],data[xind-2,yind+2],data[xind+2,yind-2]])
            elif xind+2<data.shape[0]:
                data_sub = np.mean([data[xind-2,yind-2],data[xind-2,yind],data[xind+2,yind-2]])
            elif yind+2<data.shape[0]:
                data_sub = np.mean([data[xind-2,yind-2],data[xind-2,yind+2],data[xind,yind-2]])
            # data_sub = -3
            data_masked[xind,yind]=data_sub
            # mask the pixels around the central pixel
            data_masked[xind+1,yind]=data_sub
            data_masked[xind-1,yind]=data_sub
            data_masked[xind,yind+1]=data_sub
            data_masked[xind,yind-1]=data_sub
            data_masked[xind+1,yind+1]=data_sub
            data_masked[xind-1,yind-1]=data_sub
            data_masked[xind+1,yind-1]=data_sub
            data_masked[xind-1,yind+1]=data_sub
        else:
            xind = int(round(xind))
            yind = int(round(yind))
           # print(xind,yind)
            
            # move it out of the corner if it is in the corner
            if yind > data.shape[1]-2:
                yind = data.shape[1]-2
                #print('the new yind is %s'%yind)
            if xind > data.shape[0]-2:
                xind = data.shape[0]-2
                #print('the new xind is %s'%xind)
            if yind<1:
                yind=1
                #print('the new yind is %s'%yind)
            if xind<1:
                xind=1
                #print('the new xind is %s'%xind)
                
            # now check for the highest number around
            start=-100; maxx=0; maxy=0; tomedian=[]
            for x in [xind-1,xind,xind+1]:
                for y in [yind-1,yind,yind+1]:
                    if data_masked[x,y]>start:
                        start = data_masked[x,y]
                        maxx=x
                        maxy=y
                    tomedian.append(data_masked[x,y])
                    
            data_sub = np.median(tomedian)
            data_masked[maxx,maxy]=data_sub
            print('masked %s, %s'%(maxx,maxy))
            
    return data_masked
    
def plotgals(xgal,ygal,rhgas,rhstar,ax1):
    for i in range(len(xgal)):
        circle1 = plt.Circle((xgal[i],ygal[i]), radius=rhgas[i]/1000., color='red',fill=False)
        ax1.add_artist(circle1)
        circle1 = plt.Circle((xgal[i],ygal[i]), radius=rhstar[i]/1000., color='blue',fill=False)
        ax1.add_artist(circle1)
        if verbose:
            print('%s and %s and the radius of stars is %s '%(xgal[i],ygal[i],rhstar[i]))
        
def getgalaxies_inFOV(distance):

    xmin = x_center-x_FOV[distance]/2.
    ymin = y_center-y_FOV[distance]/2.
    xmax = x_center+x_FOV[distance]/2.
    ymax = y_center+y_FOV[distance]/2.
    print xmin,xmax,ymin,ymax
        
    zmin = 10.; zmax = 15.  # using the mid z 12.5 sims
    
    return searchgals(xmin,xmax,ymin,ymax,zmin,zmax)
    
    
def getgalaxies_filament(number=1):

    if number==1:
        xmin = 46.2; xmax = 48.1; ymin = 10.5; ymax = 14.
    elif number==2:
        xmin = 44.5; xmax = 46.0; ymin = 6.9; ymax = 8.2
    elif number==3:
        xmin = 53.; xmax = 56.; ymin = 7.5; ymax = 10.  
    
    zmin = 10.; zmax = 15.  # using the mid z 12.5 sims
    
    return searchgals(xmin,xmax,ymin,ymax,zmin,zmax)
    
def searchgals(xmin,xmax,ymin,ymax,zmin,zmax,strict=False):
    
    mySim = ('RefL0100N1504',100.)
    con   = sql.connect("dlokhorst",password="mxdPB54Y")  

    myQuery  = "SELECT "+"SH.GalaxyID, \
                SH.StarFormationRate as SFR, \
                SH.CentreOfPotential_x, \
                SH.CentreOfPotential_y, \
                SH.CentreOfPotential_z, \
                SH.SubGroupNumber, \
                SH.MassType_Star, \
                SH.HalfMassProjRad_Gas, \
                SH.HalfMassProjRad_Star \
            FROM \
                %s_SubHalo as SH \
            WHERE \
                SH.SnapNum = 28 and \
                SH.CentreOfPotential_x >= %s and \
                SH.CentreOfPotential_x <= %s and \
                SH.CentreOfPotential_y >= %s and \
                SH.CentreOfPotential_y <= %s and \
                SH.CentreOfPotential_z >= %s and \
                SH.CentreOfPotential_z <= %s and \
                SH.MassType_Star > 0 "%('RefL0100N1504',xmin,xmax,ymin,ymax,zmin,zmax)
 
    if strict:
        myQuery  = "SELECT "+"SH.GalaxyID, \
                    SH.StarFormationRate as SFR, \
                    SH.CentreOfPotential_x, \
                    SH.CentreOfPotential_y, \
                    SH.CentreOfPotential_z, \
                    SH.SubGroupNumber, \
                    SH.MassType_Star, \
                    SH.HalfMassProjRad_Gas, \
                    SH.HalfMassProjRad_Star \
                FROM \
                    %s_SubHalo as SH \
                WHERE \
                    SH.SnapNum = 28 and \
                    SH.CentreOfPotential_x >= %s and \
                    SH.CentreOfPotential_x <= %s and \
                    SH.CentreOfPotential_y >= %s and \
                    SH.CentreOfPotential_y <= %s and \
                    SH.CentreOfPotential_z >= %s and \
                    SH.CentreOfPotential_z <= %s and \
                    SH.MassType_Star > 0 and \
                    SH.StarFormationRate > 0.00001 "%('RefL0100N1504',xmin,xmax,ymin,ymax,zmin,zmax)
        
                #  and \
              #  SH.StarFormationRate > 0.00001 

    if verbose:
        print myQuery
    
    myData = sql.execute_query(con,myQuery)

    xgal = myData['CentreOfPotential_x'][:]   # cMpc
    ygal = myData['CentreOfPotential_y'][:]   # cMpc
    #z = myData['CentreOfMass_z'][:]
    mgal = myData['MassType_Star'][:]         # M_solar
    rhgas = myData['HalfMassProjRad_Gas'][:]  # pkpc
    rhstar= myData['HalfMassProjRad_Star'][:] # pkpc
    
    return xgal,ygal,mgal,rhgas,rhstar
    
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
    print("Will reduce resolution by a factor of %s."%factor)
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

def define_filament_boxes(data,size=100.):
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
    xboxes, yboxes = define_filament_boxes(datares)
    xfull, yfull= get_halpha_SB.indices_region(xboxes[boxnum].astype(int),yboxes[boxnum].astype(int)) 
    SBdata = extractdata(xfull,yfull,datares)
    return SBdata

#----------------------------------------------- PLOTTING FUNCTIONS  ---------------------------------------------------------#

def plotfilamentboxes(ax1,snapnum=28):
    if snapnum == 27:
        # Pick out regions along the filaments in the map (snapnum 27)
        ax1.plot([53,53,56,56,53],[9.2,10,8.5,7.7,9.2],color='r', label='Region 3')
        ax1.plot(np.array([46.2,47.2,48.2,47.2,46.2]),[14,14,10,10,14],color='r', label='Region 1')
        ax1.plot(np.array([43,43,46,46,43]),[7.5,8.2,7.2,6.5,7.5],color='r', label = 'Region 2')
    
        xbox_3 = np.array([53,53,56,56])*3200./100.
        ybox_3 = np.array([9.2,10,8.5,7.7])*3200./100.

        xbox_2 = np.array([43,43,46,46])*3200./100.
        ybox_2 = (np.array([7.8,8.,7.1,6.8])-0.05)*3200./100.

        xbox_1 = np.array([47.4,46.2,46.9,48.1])*3200./100.
        ybox_1 = np.array([10.5,14,14,10.5])*3200./100.
    
    if snapnum == 28:
        # Pick out regions along the filaments in the map (snapnum 28)

        ax1.plot([53,53,56,56,53],np.array([9.2,10,8.5,7.7,9.2])-0.2,color='r', label='Region 3')
        ax1.plot(np.array([46.2,47.2,48.2,47.2,46.2])+0.5,[14,14,10,10,14],color='r', label='Region 1')
        ax1.plot(np.array([43,43,46,46,43]),np.array([7.5,8.2,7.2,6.5,7.5])+0.2,color='r', label = 'Region 2')
    
        xbox_3 = np.array([53,53,56,56])*3200./100.
        ybox_3 = (np.array([9.2,10,8.5,7.7])-0.2)*3200./100.

        xbox_2 = np.array([43,43,46,46])*3200./100.
        ybox_2 = (np.array([7.8,8.,7.1,6.8])-0.05+0.2)*3200./100.

        xbox_1 = (np.array([47.4,46.2,46.9,48.1])+0.5)*3200./100.
        ybox_1 = np.array([10.5,14,14,10.5])*3200./100.

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


def plotit(exptime, ax1, mymap='gist_gray', res=500, label=''):
    
    addnoisesqrt = False

    #  Add the noise to the EAGLE SB data
    try:
        if res>99 and res<101:
            resolution = 100. #arcsec
            SBdata = np.load('SBdata_100arcsec.npz')['arr_0']
        elif res>499 and res < 501:
            resolution = 500. #arcsec
            SBdata = np.load('SBdata_500arcsec.npz')['arr_0']   
        elif res>999 and res < 1001:
            resolution = 1000. #arcsec
            SBdata = np.load('SBdata_1000arcsec.npz')['arr_0']
        else:
            print "Unsupported resolution"
    except:
        print "Unsupported resolution"
        
    SBdata_exp0 = HalphaSBplot_addnoise.addnoise(SBdata,resolution,exptime=exptime,CMOS=True)

    # Plot the subtracted noiseadded data
    #fig = plt.figure(figsize = (9.5, 5.))
    #ax1 = plt.subplot(111)

    # Plot the data nicely
    median = np.median(SBdata_exp0);
    sig = np.sqrt(median)

    mymax = median + 40*sig
    mymin = median - 5*sig

    SBdata_clipped = SBdata_exp0
    SBdata_clipped[SBdata_clipped < mymin] = mymin
    SBdata_clipped[SBdata_clipped > mymax] = mymax
    SBdata_clipped = SBdata_clipped - mymin

    get_halpha_SB.makemapfilament(np.log10(SBdata_clipped**0.5),ax1,contours=False,mockobs=True,colmap=mymap,label=label,labelaxes=True)

def plotit_filament(SBdata,exptime, ax1, mymap='gist_gray', label=''):
    
    addnoisesqrt = False
        
    SBdata_exp0 = HalphaSBplot_addnoise.addnoise(SBdata,resolution,exptime=exptime,CMOS=True)

    # Plot the subtracted noiseadded data
    #fig = plt.figure(figsize = (9.5, 5.))
    #ax1 = plt.subplot(111)

    # Plot the data nicely
    median = np.median(SBdata_exp0);
    sig = np.sqrt(median)

    mymax = median + 40*sig
    mymin = median - 5*sig

    SBdata_clipped = SBdata_exp0
    SBdata_clipped[SBdata_clipped < mymin] = mymin
    SBdata_clipped[SBdata_clipped > mymax] = mymax
    SBdata_clipped = SBdata_clipped - mymin

    get_halpha_SB.makemapfilament(np.log10(SBdata_clipped**0.5),ax1,contours=False,mockobs=True,colmap=mymap,label=label,labelaxes=True)


def plotit_general(data, size, xystarts, ax1, exptime, mymap='gist_gray', res=100, label=''):
    
    data_obs = HalphaSBplot_addnoise.addnoise(data,resolution,exptime=exptime,CMOS=True)
    
    # Plot the subtracted noiseadded data
    #fig = plt.figure(figsize = (9.5, 5.))
    #ax1 = plt.subplot(111)
    
    # Plot the data nicely
    median = np.median(data_obs);
    sig = np.sqrt(median)
    
    mymax = median + 40*sig
    mymin = median - 5*sig
    
    SBdata_clipped = data_obs
    SBdata_clipped[SBdata_clipped < mymin] = mymin
    SBdata_clipped[SBdata_clipped > mymax] = mymax
    SBdata_clipped = SBdata_clipped - mymin

    get_halpha_SB.makemap(np.log10(SBdata_clipped**0.5),size,ax1,xystarts=xystarts,colmap=mymap,label=label,colorbar=False)


#------------------------------------------- BODY OF PROGRAM STARTS HERE ---------------------------------------------#

if __name__ == "__main__":
    
    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    debugging       = arguments['--debug']
    saveloc         = arguments['--saveloc']

    plotchecks      = arguments['--plotchecks']
    
    resolution      = arguments['--res']
    exptime         = arguments['--exptime']
    mockobs         = arguments['--mockobs']
    cmos            = arguments['--cmos']
    
    maskgal         = arguments['--maskgal']
    
    if verbose:
        print arguments
    
    distance = '50Mpc'  ### '50Mpc' '100Mpc' '200Mpc' '500Mpc'
    boxnum = '1' ### which filament (there are 3)
    factor = 1
    machine='coho'

    def plotforpaper():
        #label = 'Warm/Hot gas filament in the cosmic web from the EAGLE simulation'
        #map='bone'
        #plotit(float(exptime), map, float(resolution), label = label)
        #plt.show()
        
        
        ##### figure with full resolution EAGLE SB in filament in first plot, then mock observations at different resolutions following #####
        fig,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,figsize = (8.5, 7.))
        hour = 3600; map = 'viridis'; label = 'Warm/Hot gas filament in the cosmic web from the EAGLE simulation'
        SBdata = np.load('SBdata_full.npz')['arr_0']
        get_halpha_SB.makemapfilament(SBdata,ax1,colmap=map,labelaxes=True)#,label=label)
        map = 'bone'; label = '100" resolution'
        exptime = 10**3 * 3600
        plotit(exptime,ax2, map, 100)
        #ax2.locator_params(axis='y', nticks=3)
        ax2.set_xticklabels([])
        ax2.set_xlabel('')
        #ax2.set_xticks([])
        label = '500" resolution'
        plotit(exptime,ax3, map, 500)
        ax3.set_xticklabels([])
        ax3.set_xlabel('')
        #ax3.set_xticks([])
        label = '1000" resolution'
        plotit(exptime,ax4, map, 1000)
        ax4.set_xticklabels([])
        ax4.set_xlabel('')
        #ax4.set_xticks([])
        #plt.tight_layout()
        plt.subplots_adjust( hspace=0.2)
    
    plotforpaper()    
    plt.show()

    def maskgalaxies_infilament(resolution,distance):
        # load the full 100Mpc x 100Mpc box size data first
        data_tuple = (np.load('data_%s_%sarcsec.npz'%(distance,resolution))['arr_0'])
        factor = data_tuple[2]; newsize = data_tuple[1]; data = data_tuple[0];
        
        # select out a region around the filament 1
        xmin = np.min(xbox_1)-1. ; xmax = np.max(xbox_1)+1. ; 
        ymin = np.min(ybox_1)-1. ; ymax = np.max(ybox_1)+1.
        xystarts = [xmin-pixscale[distance]*resolution/2.,ymin-pixscale[distance]*resolution/2.]  ### 
        size     = [xmax-xmin,ymax-ymin]
        data_aroundfilament1 = data[(xmin/100.*(newsize/factor)):(xmax/100.*(newsize/factor)),(ymin/100.*(newsize/factor)):(ymax/100.*(newsize/factor))]
        
        if plotchecks:
            fig = plt.figure(figsize=[5,5])
            axis = plt.subplot(111)
            get_halpha_SB.makemap(data_aroundfilament1,size,axis,xystarts = xystarts)
            plt.show()
        
        # load the galaxy catalogue
        zmin = 10.; zmax = 15.
        xgal,ygal,mgal,rhgas,rhstar = searchgals(xmin,xmax,ymin,ymax,zmin,zmax)
        
        if plotchecks:
            fig = plt.figure(figsize=[5,5])
            axis = plt.subplot(111)
            get_halpha_SB.makemap(data_aroundfilament1,size,axis,xystarts = xystarts)
            plotgals(xgal,ygal,rhgas,rhstar,axis)
            plt.show()
            
        # mask the galaxies in the EAGLE data
        data_masked = maskgals2(xgal,ygal,rhgas,data,[0,0],[100.,100.],distance)  # using the full 100Mpcx100Mpc data for the next steps
        
        if plotchecks:
            fig = plt.figure(figsize=[5,10])
            axis = plt.subplot(121)
            get_halpha_SB.makemap(data_masked[(xmin/100.*(newsize/factor)):(xmax/100.*(newsize/factor)),(ymin/100.*(newsize/factor)):(ymax/100.*(newsize/factor))],
                                    size,axis,xystarts = xystarts)
            plotgals(xgal,ygal,rhgas,rhstar,axis)
            
            axis2 = plt.subplot(122)
            get_halpha_SB.makemap(data[(xmin/100.*(newsize/factor)):(xmax/100.*(newsize/factor)),(ymin/100.*(newsize/factor)):(ymax/100.*(newsize/factor))],
                                    size,axis2,xystarts = xystarts)
            plotgals(xgal,ygal,rhgas,rhstar,axis2)
            plt.show()
        
        # select out the filament region
        xboxes, yboxes = define_filament_boxes(data_masked)
        xfull, yfull= get_halpha_SB.indices_region(xboxes[boxnum].astype(int),yboxes[boxnum].astype(int)) 
        SBdata = extractdata(xfull,yfull,data_masked)
        
        if plotchecks:
            SBdata_notmasked = np.load('SBdata_%sarcsec.npz'%resolution)['arr_0']
            fig,(ax1,ax2) = plt.subplots(2,1,figsize = (8.5, 7.))
            map='viridis'
            get_halpha_SB.makemapfilament(SBdata_notmasked,ax1,colmap=map,labelaxes=True)
            get_halpha_SB.makemapfilament(SBdata,ax2,colmap=map,labelaxes=True)
            plt.show()
        
        np.savez('SBdata_%sarcsec_%s_masked.npz'%(resolution,distance))
        
        return SBdata
    
    
    def loadmaskeddata(resolution,distance):    
        # load the data if it exists or else make it and save it for later
        SB_fname = 'SBdata_%sarcsec_%s_masked.npz'%(resolution,distance)
        if os.path.isfile(SB_fname):
            print("data tuple exists, loading %s now..."%fname)
            SBdata = (np.load(SB_fname)['arr_0'])
        else:
            data_fname = 'data_%s_%sarcsec.npz'%(distance,resolution)
            print("data not saved, loading %s to make it now..."%data_fname)
            if os.path.isfile(data_fname):
                print("%s is saved, using it now to make SBdata..."%data_fname)
                SBdata = maskgalaxies_infilament(resolution,distance)
            else:
                print("data_%s_%sarcsec.npz does not exist, need to make that first..."%(distance,resolution))
                print("loading from total res, 5Mpc slice, file now...")
                total_fname = 'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_total.npz'
                if os.path.isfile(total_fname):
                    print("data exists, loading %s now..."%total_fname)
                    sl = [slice(None,None,None), slice(None,None,None)]
                    data = (np.load(total_fname)['arr_0'])[sl]
                else:
                    print("data not saved, loading from original files now...")
                    data = loaddata()
                    np.savez(total_fname,data)
                
                data_tuple = changeres(distance,resolution,data)
                np.savez(data_fname,data_tuple)
                
                SBdata = maskgalaxies_infilament(resolution,distance)
        
        return SBdata
    
    if maskgal:
        
        ## Masking after the binning has been done
        
        plotchecks=False
        verbose=False
        
        resolution = 100; distance = '50Mpc'
        # load the data
        SBdata_100 = maskgalaxies_infilament(resolution,distance)
        # plot the data
        fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize = (8.5, 7.))
        map = 'bone'; label = '100" resolution'; exptime = 10**4 * 3600
        plotit_filament(SBdata_100,exptime, ax1, mymap='gist_gray', label='')
        
        resolution = 500; distance = '50Mpc'
        # load the data if it exists or else make it and save it for later
        SBdata_500 = maskgalaxies_infilament(resolution,distance)
        # plot the data, now masked
        plotit_filament(SBdata_500,exptime, ax2, mymap='gist_gray', label='')
        
        resolution = 1000; distance = '50Mpc'
        # load the data if it exists or else make it and save it for later
        SBdata_1000 = maskgalaxies_infilament(resolution,distance)
        
        # plot the data, now masked
        plotit_filament(SBdata_1000,exptime, ax3, mymap='gist_gray', label='')
        plt.show()

        
    testplot=False
    if testplot:
        # load a slice of data
        data_5 = loaddata() # load in data at full resolution
        
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
    
    xboxes, yboxes = define_filament_boxes(data_50Mpc_500arcsec)
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