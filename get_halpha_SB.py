import numpy as np
import eagle_constants_and_units as c
import cosmo_utils as csu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
from astropy import constants as const
from astropy import units as u
#%matplotlib inline
from astropy.stats import sigma_clip

def imreduce_masked(img,mask,factor):
    """
    reduces the resolution of an image by taking the mean of individual elements
    takes in a mask to mask out values to be not included in the mean
        img: 2D image array
        mask: mask for the 2D image array
        factor: factor by which to reduce the number of array elements along each axis
    examples:
    for testing: 
        image = np.array([[1,1,2,2],[1,201,2,2],[3,3,200,4],[3,3,4,4]])
    mask your mask like this:  
        clipped = sigma_clip(image,sig=3,iters=2)
        mask = clipped.mask
    for testing:
        image = np.array([[201,201,1,1],[201,201,1,1],[2,2,3,3],[2,2,3,3]])
        clipped = sigma_clip(image,sig=1,iters=1)
        mask = clipped.mask
    """
    
    inshape = np.array(img.shape)
    
    inimg = img
    inmask = mask
    
    if np.sum(inshape%factor) != 0:
        print('Output grid must have a integer number of cells: cannot reduce image pixels by a factor %i'%factor)
        return None
    
    # split along axes into groups that will be binned
    inimg = np.array(np.split(inimg,inshape[0]/factor,axis=0))
    inimg = np.array(np.split(inimg,inshape[1]/factor,axis=-1))
    # do the same for the masks
    inmask = np.array(np.split(inmask,inshape[0]/factor,axis=0))
    inmask = np.array(np.split(inmask,inshape[1]/factor,axis=-1))
    
    # make the masked array
    x = np.ma.array(inimg, mask=inmask)
    
    # take the mean along different axes
    x = np.ma.mean(x,axis=-1)
    x = np.ma.mean(x,axis=-1)
    
    # BUT WHAT IF THERE IS ONLY MASKED DATA WITHIN A BIN...
    if True in x.mask:
        print "WARNING: At least one bin contains only masked data - will be filled to -999."
        x = x.filled(-999)
        return x.T
    
    outimg = x.data
    
    return outimg.T

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

def makemap(data,size,ax,xystarts = [0.,0.],title = ''):
    fontsize=13
    #xystarts = [0.,0.] # lower left origin of the plot
    Vmin = None
    Vmax = None
    
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
    ax.set_xlabel(r'X [cMpc]',fontsize=fontsize)
    ax.set_ylabel(r'Y [cMpc]',fontsize=fontsize)
    ax.minorticks_on()
    ax.tick_params(labelsize=fontsize)
    
    colmap = 'viridis' #'afmhot'
    ax.patch.set_facecolor(cm.get_cmap(colmap)(0.)) # sets background color to lowest color map value
    
    # nearest neighbour interpolation does not do any averaging, it just picks the nearest point and uses that as the value for a specific section in the image
    img = ax.imshow(data.T,extent=(xystarts[0],xystarts[0]+xsize,xystarts[1],xystarts[1]+ysize),origin='lower', cmap=cm.get_cmap(colmap),interpolation='nearest') # vmin = None, vmax=Vmax,
    
    plt.title(title,fontsize=fontsize)
    div = axgrid.make_axes_locatable(ax)
    cax = div.append_axes("right",size="5%",pad=0.1)
    cbar = plt.colorbar(img, cax=cax)
    
    clabel = r'log photons/cm$^2$/s/sr'
    cbar.solids.set_edgecolor("face")
    cbar.ax.set_ylabel(r'%s' % (clabel), fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    
    return img
    #plt.show()

def indices_region(xbox,ybox):
    # Tested in Eagle.ipynb on coho
    
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
