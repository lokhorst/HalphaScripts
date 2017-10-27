"""
HalphaSB_reduce_resolution.py

purpose:  reduce the resolution of the EAGLE simulations (at several radial distances) to match the pixel resolution 
of the Dragonfly telescope (binned to 100" pixel sizes).


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

# Factor by which we reduce the simulation, pixel scale values are from the online NED cosmology calculator
pixel_scale_50Mpc   = 0.237/1000. # Mpc per arcsec, z = 0.0115
pixel_scale_100Mpc  = 0.477/1000. # z = 0.0235
pixel_scale_200Mpc  = 0.928/1000. # z = 0.047
pixel_scale_1000Mpc = 2.178/1000. # z = 0.12

print('At 50, 100, 200, and 1000 Mpc, have %s, %s, %s, and %s Mpc per 100\" '%(pixel_scale_50Mpc*100.,pixel_scale_100Mpc*100.,pixel_scale_200Mpc*100.,pixel_scale_1000Mpc*100.))

simboxsize = 100. # Mpc
simpixsize = 100./32000. # Mpc / pixel
factor_50Mpc = (pixel_scale_50Mpc * 100.)/simpixsize  # will give the number of pixels of the simulation need to bin to reach size of 100" on the sky
factor_100Mpc = (pixel_scale_100Mpc * 100.)/simpixsize  # will give the number of pixels of the simulation need to bin to reach size of 100" on the sky
factor_200Mpc = (pixel_scale_200Mpc * 100.)/simpixsize  # will give the number of pixels of the simulation need to bin to reach size of 100" on the sky
factor_1000Mpc = (pixel_scale_1000Mpc * 100.)/simpixsize  # will give the number of pixels of the simulation need to bin to reach size of 100" on the sky

print('At 50, 100, 200, and 1000 Mpc, fit %s, %s, %s, and %s simulation pixels into 100\" '%(factor_50Mpc,factor_100Mpc,factor_200Mpc,factor_1000Mpc))

for factor in [factor_50Mpc,factor_100Mpc,factor_200Mpc,factor_1000Mpc]:
    if factor < 1.:
        print('100\" resolution is below the resolution of the simulation.')
    else:
        print('bin by %s resolution pixels'%int(round(factor)))

"""
bin by 7.584 resolution pixels    -----> ~8
bin by 15.264 resolution pixels   -----> ~15
bin by 29.696 resolution pixels   -----> ~30
bin by 69.696 resolution pixels   -----> ~70
"""


# Great, we round those so that the imreduce script will work (hopefully).
# We can load up the snapnum 28 full size.  It would probably be best to first show the plot with the full simulation (maybe 50 Mpc across in size since probably don't need the 100Mpc) with the Dragonfly FOV all shown on top.  Then individually replot for each distance only the Dragonfly FOV at that distance and binned appropriately.

# snapnum 28, z = 0, width = 5Mpc, centred at z = 12.5Mpc, no SFing included, 100Mpc box size, Halpha emission ---- just testing distances
file_snapnum28_noSF = '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_noSFR.npz'
sl = [slice(None,None,None), slice(None,None,None)] 
data = (np.load(file_snapnum28_noSF)['arr_0'])[sl]

data_50 = get_halpha_SB.imreduce(data, int(round(factor_50Mpc)), log=True, method = 'average')
# For the other size FOV, the factors are not integer multiples of 32000., therefore I'll trim the data first and then imreduce it
#data_trim = data[0:30000,0:30000]
#data_100 = get_halpha_SB.imreduce(data_trim, int(round(factor_100Mpc)), log=True, method = 'average')
#data_trim = data[0:31800,0:31800]
#data_200 = get_halpha_SB.imreduce(data_trim, int(round(factor_200Mpc)), log=True, method = 'average')
#data_trim = data[0:31500,0:31500]
#data_1000 = get_halpha_SB.imreduce(data_trim, int(round(factor_1000Mpc)), log=True, method = 'average')

#  Full Plot  #
# Use data_50 because if don't bin at all, end up with some weird bright pixels that are not physical

xystarts = [38.,0.]
size     = [30.,18.]
factor = int(round(factor_50Mpc))
fig = plt.figure(figsize = (16.5, 15.))
ax1 = plt.subplot(111)
get_halpha_SB.makemap(data_50[(38./100.*(32000./factor)):(68./100.*(32000./factor)),(0./100.*(32000./factor)):(18./100.*(32000./factor))],size,ax1,xystarts = xystarts)

# Plot the FOV of dragonfly on top
# snapnum 28:  
scale_50  = 0.206 #kpc/"  ## at z = 0.01 (about 50 Mpc)
scale_100 = 0.467 #kpc/"  ## at z = 0.023 (about 100 Mpc)
scale_200 = 0.928 #kpc/"  ## at z = 0.047. (about 200 Mpc)
scale_1000 = 2.178 #kpc/" ## at z = 0.12 (about 1000 Mpc)

x_angFOV = 3.*60.*60. #" 
y_angFOV = 2.*60.*60. #"  
x_FOV_50  = x_angFOV * scale_50 / 1000. # Mpc
y_FOV_50  = y_angFOV * scale_50 / 1000. # Mpc
x_FOV_100 = x_angFOV * scale_100 / 1000. # Mpc
y_FOV_100 = y_angFOV * scale_100 / 1000. # Mpc
x_FOV_200 = x_angFOV * scale_200 / 1000. # Mpc
y_FOV_200 = y_angFOV * scale_200 / 1000. # Mpc
x_FOV_1000 = x_angFOV * scale_1000 / 1000. # Mpc
y_FOV_1000 = y_angFOV * scale_1000 / 1000. # Mpc

x_center = 53.
y_center = 9.
  
ax1.plot([x_center-x_FOV_50/2.,x_center+x_FOV_50/2.,x_center+x_FOV_50/2.,x_center-x_FOV_50/2.,x_center-x_FOV_50/2.],
         [y_center+y_FOV_50/2.,y_center+y_FOV_50/2.,y_center-y_FOV_50/2.,y_center-y_FOV_50/2.,y_center+y_FOV_50/2.],color='gray')
ax1.plot([x_center-x_FOV_100/2.,x_center+x_FOV_100/2.,x_center+x_FOV_100/2.,x_center-x_FOV_100/2.,x_center-x_FOV_100/2.],
         [y_center+y_FOV_100/2.,y_center+y_FOV_100/2.,y_center-y_FOV_100/2.,y_center-y_FOV_100/2.,y_center+y_FOV_100/2.],color='gray')
ax1.plot([x_center-x_FOV_200/2.,x_center+x_FOV_200/2.,x_center+x_FOV_200/2.,x_center-x_FOV_200/2.,x_center-x_FOV_200/2.],
         [y_center+y_FOV_200/2.,y_center+y_FOV_200/2.,y_center-y_FOV_200/2.,y_center-y_FOV_200/2.,y_center+y_FOV_200/2.],color='gray')
ax1.plot([x_center-x_FOV_1000/2.,x_center+x_FOV_1000/2.,x_center+x_FOV_1000/2.,x_center-x_FOV_1000/2.,x_center-x_FOV_1000/2.],
         [y_center+y_FOV_1000/2.,y_center+y_FOV_1000/2.,y_center-y_FOV_1000/2.,y_center-y_FOV_1000/2.,y_center+y_FOV_1000/2.],color='gray')
ax1.text(x_center-x_FOV_50/2.,y_center+y_FOV_50/2.,'50 Mpc away')
ax1.text(x_center-x_FOV_50/2.,y_center+y_FOV_100/2.,'100 Mpc away')
ax1.text(x_center-x_FOV_50/2.,y_center+y_FOV_200/2.,'200 Mpc away')
ax1.text(x_center-x_FOV_50/2.,y_center+y_FOV_1000/2.,'1000 Mpc away')

plt.savefig('HalphaSB_DragonflyFOVs_together_v3.pdf')
plt.close()

# Individual Plots #

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (16.5, 15.))

# Need to test to see if this works -- CHANGED makemap in get_halpha_SB so need to send new version to chinook
xystarts = [x_center-x_FOV_50/2.,y_center-y_FOV_50/2.]
size     = [x_FOV_50, y_FOV_50]
factor = int(round(factor_50Mpc))
fig = plt.figure(figsize = (16.5, 15.))
ax1 = plt.subplot(111)
get_halpha_SB.makemap(data_50[((x_center-x_FOV_50/2.)/100.*(32000./factor)):((x_center+x_FOV_50/2.)/100.*(32000./factor)),((y_center-y_FOV_50/2.)/100.*(32000./factor)):((y_center+y_FOV_50/2.)/100.*(32000./factor))],size,ax1,xystarts = xystarts)
plt.savefig('HalphaSB_DragonflyFOV_50Mpc_v2.pdf')

del data_50

data_trim = data[0:30000,0:30000]
data_100 = get_halpha_SB.imreduce(data_trim, int(round(factor_100Mpc)), log=True, method = 'average')

xystarts = [x_center-x_FOV_100/2.,y_center-y_FOV_100/2.]
size     = [x_FOV_100, y_FOV_100]
factor = int(round(factor_100Mpc))
numpix = data_100.shape[0]
simsize = 30000./32000. * 100.
fig = plt.figure(figsize = (16.5, 15.))
ax1 = plt.subplot(111)
get_halpha_SB.makemap(data_100[((x_center-x_FOV_100/2.)/simsize*numpix):((x_center+x_FOV_100/2.)/simsize*numpix),((y_center-y_FOV_100/2.)/simsize*numpix):((y_center+y_FOV_100/2.)/simsize*numpix)],size,ax1,xystarts = xystarts)
plt.savefig('HalphaSB_DragonflyFOV_100Mpc_v2.pdf')

del data_100

data_trim = data[0:31800,0:31800]
data_200 = get_halpha_SB.imreduce(data_trim, int(round(factor_200Mpc)), log=True, method = 'average')

xystarts = [x_center-x_FOV_200/2.,y_center-y_FOV_200/2.]
size     = [x_FOV_200, y_FOV_200]
factor = int(round(factor_200Mpc))
numpix = data_200.shape[0]
simsize = 31800./32000. * 100.
fig = plt.figure(figsize = (16.5, 15.))
ax1 = plt.subplot(111)
get_halpha_SB.makemap(data_200[((x_center-x_FOV_200/2.)/simsize*numpix):((x_center+x_FOV_200/2.)/simsize*numpix),((y_center-y_FOV_200/2.)/simsize*numpix):((y_center+y_FOV_200/2.)/simsize*numpix)],size,ax1,xystarts = xystarts)
plt.savefig('HalphaSB_DragonflyFOV_200Mpc_v2.pdf')

del data_200

data_trim = data[0:31500,0:31500]
data_1000 = get_halpha_SB.imreduce(data_trim, int(round(factor_1000Mpc)), log=True, method = 'average')

xystarts = [x_center-x_FOV_1000/2.,y_center-y_FOV_1000/2.]
size     = [x_FOV_1000, y_FOV_1000]
factor = int(round(factor_1000Mpc))
numpix = data_1000.shape[0]
simsize = 31500./32000. * 100.
fig = plt.figure(figsize = (16.5, 15.))
ax1 = plt.subplot(111)
get_halpha_SB.makemap(data_1000[((x_center-x_FOV_1000/2.)/simsize*numpix):((x_center+x_FOV_1000/2.)/simsize*numpix),((y_center-y_FOV_1000/2.)/simsize*numpix):((y_center+y_FOV_1000/2.)/simsize*numpix)],size,ax1,xystarts = xystarts)
plt.savefig('HalphaSB_DragonflyFOV_1000Mpc_v2.pdf')

# After the above is working, we can run this for loop to do all four plots nicely -- FIRST check and make sure 1000Mpc away is working since it's gonna go overbounds
axes = [ax1, ax2, ax3, ax4]
x_FOV_list = [x_FOV_50, x_FOV_100, x_FOV_200, x_FOV_1000]
y_FOV_list = [y_FOV_50, y_FOV_100, y_FOV_200, y_FOV_1000]
data_list  = [data_50, data_100, data_200, data_1000]
for ax, x_FOV, y_FOV, data_this in zip(axes, x_FOV_list, y_FOV_list, data_list):  ## maybe change zip to izip (from itertools import izip)
    xystarts = [x_center-x_FOV/2.,y_center-y_FOV/2.]
    size     = [x_FOV, y_FOV]
    get_halpha_SB.makemap(data_this[((x_center-x_FOV/2.)/100.*3200.):((x_center+x_FOV/2.)/100.*3200.),((y_center-y_FOV/2.)/100.*3200.):((y_center+y_FOV/2.)/100.*3200.)],size,ax,xystarts = xystarts)
    
plt.savefig('HalphaSB_DragonflyFOVs.pdf')
plt.close()




"""
To Do:
- send new version of get_halpha_SB to chinook
- figure out how to get github to work
"""