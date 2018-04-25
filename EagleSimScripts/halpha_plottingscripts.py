### mask_galaxies.ipynb

def pltimg(data_FOV,xystarts,size,ax=None):
    if ax is None:
        fig, (ax) = plt.subplots(1, 1, figsize=(10, 10))
    colmap = 'viridis'
    ax.patch.set_facecolor(cm.get_cmap(colmap)(0.)) # sets background color to lowest color map value
    img = ax.imshow(data_FOV.T,extent=(xystarts[0],xystarts[0]+size[0],xystarts[1],xystarts[1]+size[1]),\
                    origin='lower', cmap=cm.get_cmap(colmap),interpolation='nearest')
    div = axgrid.make_axes_locatable(ax)
    cax = div.append_axes("top",size="5%",pad=0.1)
    cbar = plt.colorbar(img, cax=cax,orientation='horizontal')

def pltcutout(data_FOV,xystarts,size,ax=None):
    """
    useful if have a large array but only want to plot a small section to look at it briefly
    INPUT   data_FOV:     2D np array
            xystarts:     [x,y] lower starting values (to define axes limits) (Mpc)
            size:         [x,y] length of each side (to define axes limits) (Mpc)
    OUTPUT  selects out a 1Mpcx1Mpc box from the center of the input 2D array
            call to pltimg - creates image of the 1Mpcx1Mpc box with axes labelled etc.
    """
    xpixsize,ypixsize = data_FOV.shape
    xlength,ylength = size
    midx=xpixsize/2; midy=ypixsize/2
    boxlength = 1. #Mpc
    pixsizex=xpixsize/xlength*boxlength
    data_cutout = data_FOV[midx-pixsizex/2:midx+pixsizex/2,midy-pixsizex/2:midy+pixsizex/2]
        
    xystarts_x = xystarts[0] + (midx-pixsizex/2)/(xpixsize/xlength)
    xystarts_y = xystarts[1] + (midy-pixsizex/2)/(ypixsize/ylength)
    size_new = [boxlength,boxlength]

    pltimg(data_cutout,[xystarts_x,xystarts_y],size_new,ax=ax)

def plotgals(xgal,ygal,rhgas,rhstar,mgal,ax1,verbose):
    for i in range(len(xgal)):
        colour = 'green'
        if mgal[i]>10**8:
            colour = 'yellow'
        if mgal[i]>10**9:
            colour = 'orange'
        if mgal[i]>10**10:
            colour = 'red'
        circle1 = plt.Circle((xgal[i],ygal[i]), radius=rhgas[i]/1000., color=colour,fill=False)
        ax1.add_artist(circle1)
        circle1 = plt.Circle((xgal[i],ygal[i]), radius=rhstar[i]/1000., color='blue',fill=False)
        ax1.add_artist(circle1)
        if verbose:
            Mpcperpix = 0.477/1000.*(1.+0.0235) * 6.4
            if (rhstar[i]*5.) > (Mpcperpix*1000.) and (rhgas[i]*5.) > (Mpcperpix*1000.):
                print("5*rhstar, %.1f, is greater than %s kpc, and has galaxy mass of %s, and a 5*rhgas of %.1f."%\
                      ((rhstar[i]*5.),(Mpcperpix*1000.),mgal[i],(rhgas[i]*5.)))
                circle1 = plt.Circle((xgal[i],ygal[i]), radius=rhstar[i]/1000.*5., color='purple',fill=False)
                ax1.add_artist(circle1)

def plotdata(data,ax=None,bounds=None,colorbar=False,colmap='viridis'):
    """
    General use plotting of image
    """
    if ax is None:
        fig = plt.figure(figsize=(6, 3.2))
        ax = fig.add_subplot(111)
        oneplot=True
    if bounds is None:
        img = ax.imshow(data,origin='lower',cmap=cm.get_cmap(colmap),interpolation='nearest')
    else:
        img = ax.imshow(data,origin='lower',cmap=cm.get_cmap(colmap),vmin=bounds[0],vmax=bounds[1],interpolation='nearest')
    ax.set_aspect('equal')
    
    ax.patch.set_facecolor(cm.get_cmap(colmap)(0.)) # sets background color to lowest color map value
    ax.patch.set_facecolor('black')
    
    if colorbar:
        div = axgrid.make_axes_locatable(ax)
        cax = div.append_axes("right",size="10%",pad=0.15)
        cbar = plt.colorbar(img,cax=cax,orientation='vertical')#,boundaries=np.linspace(0,90000))
        cbar.ax.tick_params()

### extract_filament.ipynb

def plotfilament(SBdata_5,ax):
    clabel = r'log photons/cm$^2$/s/sr'; Vmin = None; Vmax= None
    fontsize=13
    ax.set_xlabel(r'X [cMpc]',fontsize=fontsize)
    ax.set_ylabel(r'Y [cMpc]',fontsize=fontsize)
    
    ax.tick_params(labelsize=fontsize)
    colmap = 'viridis' #'afmhot'
    ax.patch.set_facecolor(cm.get_cmap(colmap)(0.)) # sets background color to lowest color map value

    img = ax.imshow(SBdata_5,origin='lower', cmap=cm.get_cmap(colmap), vmin = Vmin, vmax=Vmax,interpolation='nearest')

    div = axgrid.make_axes_locatable(ax)
    cax = div.append_axes("right",size="15%",pad=0.1)
    cbar = plt.colorbar(img, cax=cax)
    cbar.solids.set_edgecolor("face")
    cbar.ax.set_ylabel(r'%s' % (clabel), fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
