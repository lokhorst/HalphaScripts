import numpy as np
import eagle_constants_and_units as c
import cosmo_utils as csu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
from astropy import constants as const
from astropy import units as u
#%matplotlib inline

def dothis():
    print("it worked")

def makemapfilament(SBdata,ax,colmap='viridis',onlyyellow=False,contours=False,colorbar=False,mockobs=False,labelaxes=False,label=''):
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
    
        #img = ax.imshow(data.T,extent=(xystarts[0],xystarts[0]+xsize,xystarts[1],xystarts[1]+ysize),origin='lower', cmap=cm.get_cmap(colmap),interpolation='nearest') # vmin = None, vmax=Vmax,
    
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
        #cbar.ax.set_ylabel(r'%s' % (clabel), fontsize=fontsize)
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

def imreduce(img, factor, log=True, method = 'average'):
    """
        img: 2D image array
        factor: factor by which to reduce the number of array elements along each axis
        log: whether or not the array contains log data values
        """
    if log:
        inimg = 10**img
    else:
        inimg = img
    inshape = np.array(img.shape)

    if np.sum(inshape%factor) != 0:
        print('Output grid must have a integer number of cells: cannot reduce image pixels by a factor %i'%factor)
        return None
    inimg = np.array(np.split(inimg,inshape[0]/factor,axis=0))
    inimg = np.array(np.split(inimg,inshape[1]/factor,axis=-1))

    inimg = np.sum(inimg,axis=-1)
    inimg = np.sum(inimg,axis=-1)
    
    if method == 'average':
        inimg = inimg/np.float(factor**2)
        #outimg = np.average(inimg[])
    if log:
        inimg = np.log10(inimg)
    return inimg.T

def makemap(data,size,ax,colmap='viridis',xystarts = [0.,0.],title = '',colorbar=True,mockobs=False,labelaxes=True,label=''):
    fontsize=13
    #xystarts = [0.,0.] # lower left origin of the plot
    Vmin = None
    Vmax = None
    
    if mockobs:
        clabel = r'log signal (photons)'
    else:
        clabel = r'log photons/cm$^2$/s/sr'
    
    if type(size) == float or type(size) == int:
        print('The type of size is '+str(type(size)))
        xsize = size
        ysize = size
    else:
        print('The type of size is '+str(type(size)))
        xsize = size[0]
        ysize = size[1]
    
    #fig = plt.figure(figsize = (5.5, 5.)) # large size just as a trick to get higher resolution
    #fig = plt.figure(figsize = (11., 10.))
    #ax = plt.subplot(111)
    
    if labelaxes:
        ax.set_xlabel(r'X [cMpc]',fontsize=fontsize)
        ax.set_ylabel(r'Y [cMpc]',fontsize=fontsize)
        ax.tick_params(labelsize=fontsize) #,top=True,labeltop=True)
        #ax.xaxis.set_label_position('top') 
        ax.xaxis.set_label_position('bottom') 
        #ax.xaxis.tick_top()
        #ax.minorticks_on()
    else:
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
    
    #colmap = 'viridis' #'afmhot'
    ax.patch.set_facecolor(cm.get_cmap(colmap)(0.)) # sets background color to lowest color map value
    
    # nearest neighbour interpolation does not do any averaging, it just picks the nearest point and uses that as the value for a specific section in the image
    img = ax.imshow(data.T,extent=(xystarts[0],xystarts[0]+xsize,xystarts[1],xystarts[1]+ysize),origin='lower', cmap=cm.get_cmap(colmap),interpolation='nearest') # vmin = None, vmax=Vmax,
    
    #plt.title(label,fontsize=fontsize)
    div = axgrid.make_axes_locatable(ax)
        
    if colorbar:
        #cax = div.append_axes("right",size="5%",pad=0.1)
        
        # bottom color bar:
        #cax = div.append_axes("bottom",size="15%",pad=0.1)
        #cbar.ax.set_xlabel(r'%s' % (clabel), fontsize=fontsize)
        
        # top color bar:
        cax = div.append_axes("top",size="15%",pad=0.1)
        cbar = plt.colorbar(img, cax=cax,orientation='horizontal')

        cbar.ax.set_xlabel(r'%s' % (clabel), fontsize=fontsize)
        cbar.ax.xaxis.set_label_position('top')      
        cbar.ax.xaxis.set_ticks_position('top')
        
        cbar.solids.set_edgecolor("face")
        #cbar.ax.set_ylabel(r'%s' % (clabel), fontsize=fontsize)
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
    #plt.text(1.8,0.8,label,fontdict=font,horizontalalignment='center',backgroundcolor='white')

    plt.tight_layout()
    
    
    #return img
    #plt.show()

def indices_region(xbox,ybox):
    # Tested in RegionSelection_testing.ipynb on coho
    
    # Create 2d arrays for both x and y indices values to pick out the data in the regions
    #(not sure, there's probably a better way to do this...)
    
    #[53,53,56,56,53],[9.2,10,8.5,7.7,9.2]
    #xbox = [53, 53, 56,56] # lower left, upper left, upper right, lower right
    #ybox = [9.2,10,8.5,7.7] # lower left, upper left, upper right, lower right
    
    #    xbox = [1,10,10,1] # lower left, upper left, upper right, lower right
    #    ybox = [20,15,25,30]
    
    #xacross_length = np.abs(xbox[2] - xbox[1])
    #xdown_length = np.abs(xbox[1]-xbox[0])
    #yacross_length = np.abs(ybox[2] - ybox[1])
    #ydown_length = np.abs(ybox[1]-ybox[0])
    xacross_length = (xbox[2] - xbox[1])
    xdown_length = (xbox[1]-xbox[0])
    yacross_length = (ybox[2] - ybox[1])
    ydown_length = (ybox[1]-ybox[0])
    print('xacross_length,yacross_length,xdown_length,ydown_length:  '+str(xacross_length)+', '+
          str(yacross_length)+', '+str(xdown_length)+', '+str(ydown_length))
        
    if xacross_length > yacross_length:
        iterable = np.round(np.arange(xacross_length+1)*float(yacross_length)/float(xacross_length))
        yacross = [ybox[1]+y for y in iterable.astype(int)]
        xacross = [xbox[1]+x for x in xrange(xacross_length+1)]
        print('xacross_length > yacross_length:')
        print('iterable: '+str(iterable))
    elif xacross_length < yacross_length:
        iterable = np.round(np.arange(yacross_length+1)*float(xacross_length)/float(yacross_length))
        xacross = [xbox[1]+x for x in iterable.astype(int)]
        yacross = [ybox[1]+y for y in xrange(yacross_length+1)]
        print('xacross_length < yacross_length:')
        print('iterable: '+str(iterable))
    else:
        xacross = [xbox[1]+x for x in xrange(xacross_length+1)]
        yacross = [ybox[1]+y for y in xrange(yacross_length+1)]
        print('xacross_length = yacross_length')
                                                                  
    if xdown_length > ydown_length:
        iterable = np.round(np.arange(xdown_length+1)*float(ydown_length)/float(xdown_length))
        yfull = [np.array(yacross)-y for y in iterable.astype(int)]
        xfull = [np.array(xacross)-x for x in xrange(xdown_length+1)]
        print('xdown_length > ydown_length:')
        print('iterable: '+str(iterable))
    elif xdown_length < ydown_length:
        iterable = np.round(np.arange(ydown_length+1)*float(xdown_length)/float(ydown_length))
        xfull = [np.array(xacross)-x for x in iterable.astype(int)]
        yfull = [np.array(yacross)-y for y in xrange(ydown_length+1)]
        print('xdown_length < ydown_length:')
        print('iterable: '+str(iterable))
    else:
        xfull = [np.array(xacross)-x for x in xrange(xdown_length+1)]
        yfull = [np.array(yacross)-y for y in xrange(ydown_length+1)]
        print('xdown_length = ydown_length:')

    xfull = np.array(xfull)
    yfull = np.array(yfull)
        
    return xfull, yfull

def extractdata(xfull,yfull,data):
    SBdata = np.zeros(xfull.shape)
    for i in range(yfull.shape[0]):
        for j in range(yfull.shape[1]):
                SBdata[i,j]  = data[xfull[i,j],yfull[i,j]]
    return SBdata

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
