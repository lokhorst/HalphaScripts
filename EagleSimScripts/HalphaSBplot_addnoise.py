"""
Mock observations of the EAGLE simulation.

Taking an exposure time

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

def getBackground(start,end,plot=True):
    # Returns the total background flux in the wavlength interval supplied i.e. returns (flux)*(wavlength interval) 
    wavelength = []
    flux = []
    with open('/Users/lokhorst/Documents/Eagle/Gemini_skybackground.dat','r') as f:  #wavelength in nm, flux in phot/s/nm/arcsec^2/m^2
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
        print index,index+1
        print total
        total = total+(flux[index]*(wavelength[index+1]-wavelength[index]))
        
    # if spacings are even, can just take the average of the flux array and times it by the total bandwidth
    np.mean(flux[start_ind:end_ind])*(wavelength[end_ind]-wavelength[start_ind])
    
    
    print('start index and end index: %s and %s'%(start_ind,end_ind))
    print(wavelength[start_ind:end_ind]-wavelength[start_ind+1:end_ind+1])
    if plot==True:
        plt.plot(wavelength[start_ind-100:end_ind+100],flux[start_ind-100:end_ind+100])
        top = max(flux[start_ind-100:end_ind+100])
        plt.plot([start,start,end,end,start],[0,top,top,0,0])
        plt.show()
        
    return total



factor = 10
sl = [slice(None,None,None), slice(None,None,None)]

# Simulation snapnum 27 (z = 0.1), xy box size: 100Mpc, z slice width: 5Mpc, zloc: 12.5Mpc
if machine == 'chinook':
    homedir = '/Users/lokhorst/Eagle/'
elif machine == 'coho':
    homeidr = 'Users/deblokhorst/SlicesfromNastasha/'

files_SF_27 = ['/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5__fromSFR.npz']

files_noSF_27 = ['/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5_noSFR.npz']


# Simulation snapnum 28 (z = 0), xy box size: 100Mpc, z slice width: 5Mpc,
files_SF_28 = [homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5__fromSFR.npz',
               homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5__fromSFR.npz',
               homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5__fromSFR.npz',
               homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5__fromSFR.npz']

files_noSF_28 = [homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_noSFR.npz',
                 homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5_noSFR.npz',
                 homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5_noSFR.npz',
                 homedir+'emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5_noSFR.npz']

#  For a 3nm filter width, that corresponds to a redshift slice of 20 Mpc (big!!)
data_20=[]   # log( photons / cm^2 / s / sr )
ind=0
print('Create data_20 (20Mpc wide slice) from 5Mpc slices at redshift of 0...')
for fname in files_SF_28+files_noSF_28:
    ind=ind+1
    print('loading data ('+fname+')...')
    data = (np.load(fname)['arr_0'])[sl]
    data = get_halpha_SB.imreduce(data, factor, log=True, method = 'average')
    if data_20 ==[]:
        print('(1/%s) first addition to data_20...'%(len(files_SF_28+files_noSF_28)))
        data_20 = data
    else:
        print('(%s/%s) adding data to data_20...'%(ind,len(files_SF_28+files_noSF_28)))
        data_20 = np.log10(10**data_20+10**data)
    del data

if __name__ == "__main__":

    # Dragonfly info
    area_lens = np.pi*(14.3/2)**2 * 48.               # cm^2, 48 * 14.3 cm diameter lenses
    pix_size = 2.8                                    # arcsec
    ang_size_pixel  = (pix_size * (1./206265.))**2    # rad^2, the pixel size of the CCD

    tau_l = 0.85  # transmittance of the Dragonfly lens
    QE_new = 0.70  # quantum efficiency of the CMOS detector
    QE_old = 0.48     # quantum efficiency of the CCDs
    tau_f = 1.    # transmittance of the Halpha filter -- assumed for now
    B = getBackground(656.3,659.3) # *u.photon/u.second/u.arcsec**2/u.m**2  ****already multiplied by the bandwidth***
    D = 0.04  # *u.photon/u.second                             # dark current (electrons / s) 
    R_squared_old = 10.**2 # * u.photon                           # read noise (electrons)
    R_squared_new = 1.**2 # * u.photon                           # read noise (electrons)


    # calculate the signal in a set exposure time
    exptime = 60.*60.*10**5  # seconds  (10^5 hours)
    totsignal = np.log10(10**data_20 * exptime) # log( photons / cm^2 /sr )

    print('Plotting without the noise added, small region (20 Mpc box)...')
    # Plotting parameters
    xystarts = [40.,0.]
    size     = 20.
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (10, 10))
    get_halpha_SB.makemap(totsignal[(xystarts[0]/100.*3200.):((xystarts[0]+size)/100.*3200.),(xystarts[1]/100.*3200.):((xystarts[1]+size)/100.*3200.)],size,ax1,xystarts = xystarts)
    ax1.set_title('sky signal, no noise')

    # calculate the number of detected electrons
    # how much are we going to bin?  100" square, so that is 
    binpix_size = 100. # arcsec
numpixel = round((binpix_size/pix_size)**2)

detsignal = np.log10(10**totsignal * QE_old * tau_l * tau_f * area_lens * ang_size_pixel * numpixel)
print('The total detected signal per pixel (taking into account transmittance of the lens, quantum efficiency, and the pixel size) is %s.'%detsignal)
get_halpha_SB.makemap(detsignal[(xystarts[0]/100.*3200.):((xystarts[0]+size)/100.*3200.),(xystarts[1]/100.*3200.):((xystarts[1]+size)/100.*3200.)],size,ax2,xystarts = xystarts)
ax2.set_title('detected signal, no noise')

# calculate the detected noise and system noise
B_sky = B * QE_old * tau_l * tau_f * area_lens*(1/100)**2 * pix_size**2
sigma_nophotonnoise = np.sqrt( B_sky*exptime*numpixel + D*exptime*numpixel + R_squared_old*numpixel)
print('For reference, the total noise per pixel MINUS shot noise (background sky, dark current, and read noise) is %s.'%sigma_nophotonnoise)
sigma = np.log10(np.sqrt(10**detsignal + B_sky*exptime*numpixel + D*exptime*numpixel + R_squared_old*numpixel))
get_halpha_SB.makemap(sigma[(xystarts[0]/100.*3200.):((xystarts[0]+size)/100.*3200.),(xystarts[1]/100.*3200.):((xystarts[1]+size)/100.*3200.)],size,ax3,xystarts = xystarts)
ax3.set_title('noise')


#  add the detected signal (electrons) plus the noise together - this is what we see
alldata = np.log10(10**detsignal + 10**sigma)
get_halpha_SB.makemap(alldata[(xystarts[0]/100.*3200.):((xystarts[0]+size)/100.*3200.),(xystarts[1]/100.*3200.):((xystarts[1]+size)/100.*3200.)],size,ax4,xystarts = xystarts)
ax4.set_title('detected signal + noise')

plt.show()


##### **** need to do the convolution with the Dragonfly psf still... **** ####

#print('Plot an example region (20 Mpc box)...')
# Plotting parameters
#xystarts = [40.,0.]
#size     = 20.
#fig = plt.figure(figsize = (6.5, 5.))
#ax1 = plt.subplot(111)
#get_halpha_SB.makemap(data_5[(xystarts[0]/100.*3200.):((xystarts[0]+size)/100.*3200.),(xystarts[1]/100.*3200.):((xystarts[1]+size)/100.*3200.)],size,ax1,xystarts = xystarts)
#plt.savefig('../plots/HalphaSBplot_size%sMpc.pdf'%size) # match the #Mpcmap to whichever data (data_5,data_10,data_15 etc) you are plotting

