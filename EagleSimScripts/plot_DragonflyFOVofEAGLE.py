#!/usr/bin/env python

""" plots the EAGLE simulation with the Dragonfly FOV at specified distance overlaid on top.

Usage: plot_DragonflyFOVofEAGLE.py [-h] [-s] [-v] [-f] [-m] [-c] [-r RESOLUTION] [-t EXPTIME] 

Options:
    -h, --help                          Show this screen.
    -v, --verbose                       Show extra information [default: False]
    -s, --save                          Save intermediate array products [default: False]

    -f, --filament                      Want to just plot the filament
    -c, --cmos                          Use specs for new cameras in mock observation
    -m, --mockobs                       Plot a mock observation (rather than just the EAGLE simulation). 
    -o, --oneplotFOVs                   Plot one EAGLE image with Dragonfly FOVs at 50 100 200 500 Mpc on top
    -a, --allFOVs                       Plot four different EAGLE images size of Dragonfly FOV, resolution 100"x100", either raw EAGLE or mockobs

    -r RESOLUTION, --res RESOLUTION     The desired resolution of the EAGLE data image in arcsec. [default: 500]
    -t EXPTIME, --exptime EXPTIME       The desired exposure time to integrate the  EAGLEdata over in seconds. [default: 3600e3]
    -d DISTANCE, --distance DISTANCE    The desired distance for the EAGLE simulation to be 'placed'. [default: 50 Mpc]

Examples: python plot_DragonflyFOVofEAGLE.py -v -f -m -r 100 -t 3600e8

"""

import os
import numpy as np
import docopt
#import eagle_constants_and_units as c
#import cosmo_utils as csu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
#from astropy import constants as const
#from astropy import units as u

import get_halpha_SB
import HalphaSBplot_addnoise

pixscale =  {'50Mpc': 0.237/1000.*(1.+0.0115), '100Mpc': 0.477/1000.*(1.+0.0235),'200Mpc': 0.928/1000.*(1.+0.047) , '500Mpc': 2.178/1000.*(1.+0.12)} ### Mpc / arcsec (comoving)
x_center = 47.5
y_center = 12.
x_angFOV = 3.*60.*60. # " 
y_angFOV = 2.*60.*60. # "  
x_FOV = {distance: pixscale[distance]*x_angFOV for distance in ['50Mpc','100Mpc','200Mpc','500Mpc']}  # cMpc
y_FOV = {distance: pixscale[distance]*y_angFOV for distance in ['50Mpc','100Mpc','200Mpc','500Mpc']}  # cMpc

def changeres(distance,resolution,data):
    print('changing the resolution of the data')
    simpixsize = 100./32000. ### Mpc / pixel is resolution of raw data 
    factor = round(pixscale[distance]*resolution/simpixsize)
    print('factor is %s'%factor)
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

def loaddata(factor=1):
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
    #data1 = get_halpha_SB.imreduce(data1, round(factor), log=True, method = 'average')
    print('data11 ('+files_SF_28[0]+')...')
    data11 = (np.load(files_SF_28[0])['arr_0'])[sl]
    #data11 = get_halpha_SB.imreduce(data11, round(factor), log=True, method = 'average')
    print('5 Mpc slice...')
    data_5 = np.log10(10**data1+10**data11)
    print('delete data1, data11...')
    del data1
    del data11
    
    return data_5

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
    

#-------------------------------------- BODY OF PROGRAM STARTS HERE ---------------------------------------------#

if __name__ == "__main__":
    
    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    debugging       = arguments['--debug']
    saveloc         = arguments['--save']

    # required inputs
    resolution      = arguments['--res']
    exptime         = arguments['--exptime']
    distance        = arguments['--distance']
    
    # options
    mockobs         = arguments['--mockobs']
    cmos            = arguments['--cmos']
    filament        = arguments['--filament']
    oneplotFOVs     = arguments['--oneplotFOVs']
    allFOVs         = arguments['--allFOVs']
    
    if verbose:
        print arguments
    
    machine='coho'
        
    
    if mockobs and filament:
        print("congrats, you picked both mockobs and filament!")
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
        plt.show()
        
    elif filament:
        print("congrats, you picked only filament!")

        if resolution=='100' or resolution=='500' or resolution=='1000':
            SBdata = np.load('SBdata_'+resolution+'arcsec.npz')['arr_0']
        else:
            print("no data saved, need to load original data")
            
        fig = plt.figure(figsize = (9.5, 5.))
        ax1 = plt.subplot(111)
        
        label='EAGLE simulation'
        get_halpha_SB.makemapfilament(SBdata,ax1,colmap='viridis',onlyyellow=False,contours=False,colorbar=True,mockobs=False,labelaxes=True,label='')
        plt.show()
    elif oneplotFOVs:
        print("congrats, you picked oneplotFOVs!")
      
        resolution ='total'
        
        fname = 'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_'+resolution+'.npz'
        if os.path.isfile(fname):
            print("data exists, loading %s now..."%fname)
            sl = [slice(None,None,None), slice(None,None,None)]
            data = (np.load(fname)['arr_0'])[sl]
        else:
            print("data not saved, loading from original files now...")
            data = loaddata()
            np.savez(fname,data)
            
        xystarts = [38.,0.]
        size     = [30.,18.]
        #factor = int(round(factor_50Mpc))
        fig = plt.figure(figsize = (6.5, 5.))
        ax1 = plt.subplot(111)
        get_halpha_SB.makemap(data[(38./100.*(32000./factor)):(68./100.*(32000./factor)),(0./100.*(32000./factor)):(18./100.*(32000./factor))],size,ax1,xystarts = xystarts)
        
        ax1.set_xlim(38.,68.)
        ax1.set_ylim(0.,18.)
        
        for distance in ['50Mpc','100Mpc','200Mpc','500Mpc']:
            ax1.plot([x_center-x_FOV[distance]/2.,x_center+x_FOV[distance]/2.,x_center+x_FOV[distance]/2.,x_center-x_FOV[distance]/2.,x_center-x_FOV[distance]/2.],
                     [y_center+y_FOV[distance]/2.,y_center+y_FOV[distance]/2.,y_center-y_FOV[distance]/2.,y_center-y_FOV[distance]/2.,y_center+y_FOV[distance]/2.],color='gray')
            #ax1.text(x_center-x_FOV[distance]/2.,y_center+y_FOV[distance]/2.,distance+' away')
            ax1.text(x_center-x_FOV[distance]/2.,y_center,distance+' away',rotation='vertical',va='center')
        plt.show()
    
    elif allFOVs:
        print("congrats, you picked allFOVs!  These are the all the same resolution, different distances plots.")

        # load the data:
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

        # plot the data:
        def plot1():
            """
            Only four plots at different distances, either all raw EAGLE or all mock obs
            """
            # plot them all together
            #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (11.5, 10.))
            mockobs=False
            f, ((ax1, ax2, ax3, ax4)) = plt.subplots(1, 4, figsize = (15.5, 15.))
        
            for distance,axis in zip(['50Mpc','100Mpc','200Mpc','500Mpc'],[ax1,ax2,ax3,ax4]):
                data_tuple = data_dict[distance]
                factor = data_tuple[2]; newsize = data_tuple[1]; data = data_tuple[0]; resolution = 100.
            
                xystarts = [x_center-x_FOV[distance]/2.,y_center-y_FOV[distance]/2.]
                size     = [x_FOV[distance], y_FOV[distance]]
            
                x1 = ((x_center-x_FOV[distance]/2.)/100.*(newsize/factor))
                x2 = ((x_center+x_FOV[distance]/2.)/100.*(newsize/factor))
                y1 = ((y_center-y_FOV[distance]/2.)/100.*(newsize/factor))
                y2 = ((y_center+y_FOV[distance]/2.)/100.*(newsize/factor))
                data_FOV = data[int(x1):int(x2),int(y1):int(y2)]
            
                if mockobs:
                    label = 'Mock observation at '+distance; exptime=10**7*3600.; map='bone'
                    plotit_general(data_FOV,size,xystarts,axis,float(exptime), map, float(resolution), label = label)
                else:
                    #get_halpha_SB.
                    makemap(data_FOV,size,axis,xystarts = xystarts)
        
        def plot2():
            """
            Only eight plots, half raw EAGLE at different distances and half mock obs at different distances
            """
            f, ((ax1, ax2, ax3, ax4),(ax11, ax22, ax33, ax44)) = plt.subplots(2, 4, figsize = (15.5, 15.))
        
            for distance,axistop,axisbot in zip(['50Mpc','100Mpc','200Mpc','500Mpc'],[ax1,ax2,ax3,ax4],[ax11,ax22,ax33,ax44]):
                data_tuple = data_dict[distance]
                factor = data_tuple[2]; newsize = data_tuple[1]; data = data_tuple[0]; resolution = 100.
            
                xystarts = [x_center-x_FOV[distance]/2.,y_center-y_FOV[distance]/2.]
                size     = [x_FOV[distance], y_FOV[distance]]
            
                x1 = ((x_center-x_FOV[distance]/2.)/100.*(newsize/factor))
                x2 = ((x_center+x_FOV[distance]/2.)/100.*(newsize/factor))
                y1 = ((y_center-y_FOV[distance]/2.)/100.*(newsize/factor))
                y2 = ((y_center+y_FOV[distance]/2.)/100.*(newsize/factor))
                data_FOV = data[int(x1):int(x2),int(y1):int(y2)]
            
                get_halpha_SB.makemap(data_FOV,size,axistop,xystarts = xystarts)
            
                label = 'Mock observation at '+distance; exptime=10**5*3600.; map='bone'
                plotit_general(data_FOV,size,xystarts,axisbot,float(exptime), map, float(resolution), label = label)
                
                # put little boxes where zoomins are pulling out
                xystarts = [47.4,11.2]
                size = [1.0,1.2]
                axisbot.plot([47.4,47.4,47.4+1.0,47.4+1.0,47.4],[11.2,11.2+1.2,11.2+1.2,11.2,11.2],'r-')
                
                axisbot.set_xlim((x_center-x_FOV[distance]/2.),(x_center-x_FOV[distance]/2.)+x_FOV[distance])
                axisbot.set_ylim((y_center-y_FOV[distance]/2.),(y_center-y_FOV[distance]/2.)+y_FOV[distance])
                
        
        def plot3():
            """
            12 plots (eight like plot2), plus also zoomins on filaments for each distance.
            """
            f, ((ax1, ax2, ax3, ax4),(ax11, ax22, ax33, ax44),(ax111, ax222, ax333, ax444)) = plt.subplots(3, 4, figsize = (10.5, 10.))
            
            for distance,axistop,axismid,axisbot in zip(['50Mpc','100Mpc','200Mpc','500Mpc'],[ax1,ax2,ax3,ax4],[ax11,ax22,ax33,ax44],[ax111,ax222,ax333,ax444]):
                data_tuple = data_dict[distance]
                factor = data_tuple[2]; newsize = data_tuple[1]; data = data_tuple[0]; resolution = 100.
            
                xystarts = [x_center-x_FOV[distance]/2.,y_center-y_FOV[distance]/2.]
                size     = [x_FOV[distance], y_FOV[distance]]
            
                x1 = ((x_center-x_FOV[distance]/2.)/100.*(newsize/factor))
                x2 = ((x_center+x_FOV[distance]/2.)/100.*(newsize/factor))
                y1 = ((y_center-y_FOV[distance]/2.)/100.*(newsize/factor))
                y2 = ((y_center+y_FOV[distance]/2.)/100.*(newsize/factor))
                data_FOV = data[int(x1):int(x2),int(y1):int(y2)]
            
                get_halpha_SB.makemap(data_FOV,size,axistop,xystarts = xystarts)
            
                label = 'Mock observation at '+distance; exptime=10**7*3600.; map='bone'
                plotit_general(data_FOV,size,xystarts,axismid,float(exptime), map, float(resolution), label = label)
                
                # do the zoom-ins:                
                xzoomin = np.array([47.4,48.2])/100.*(newsize/factor)
                yzoomin = np.array([11.2,12.4])/100.*(newsize/factor)
                data_zoomin = data[int(xzoomin[0]):int(xzoomin[1]),int(yzoomin[0]):int(yzoomin[1])]
                
                xystarts = [47.4,11.2]
                size = [1.0,1.2]
                
                label = 'Mock observation at '+distance; exptime=10**7*3600.; map='bone'
                plotit_general(data_zoomin,size,xystarts,axisbot,float(exptime), map, float(resolution), label = label)
            
            
                
        def plot4():
            """
            4 plots, just the zoomins (because couldn't get the size ratios right in the above plot.. needs more tweaking)
            """
            f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (10.5, 10.))
            
            for distance,axisbot in zip(['50Mpc','100Mpc','200Mpc','500Mpc'],[ax1, ax2, ax3, ax4]):
                data_tuple = data_dict[distance]
                factor = data_tuple[2]; newsize = data_tuple[1]; data = data_tuple[0]; resolution = 100.
               
                # do the zoom-ins:
                xystarts = [47.4,11.2]
                size = [1.0,1.2]
                
                xzoomin = np.array([xystarts[0],(xystarts[0]+size[0])])/100.*(newsize/factor)
                yzoomin = np.array([xystarts[1],(xystarts[1]+size[1])])/100.*(newsize/factor)
                data_zoomin = data[int(xzoomin[0]):int(xzoomin[1]),int(yzoomin[0]):int(yzoomin[1])]
            
                label = 'Mock observation at '+distance; exptime=10**5*3600.; map='bone'
                plotit_general(data_zoomin,size,xystarts,axisbot,float(exptime), map, float(resolution), label = label)
            
        plot2()
        plt.tight_layout()
        plt.show()
        
        plot4()
        plt.tight_layout()
        plt.show()