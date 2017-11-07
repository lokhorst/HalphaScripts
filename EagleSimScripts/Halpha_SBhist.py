import numpy as np
import eagle_constants_and_units as c
import cosmo_utils as csu
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
from astropy import constants as const
from astropy import units as u

import get_halpha_SB

factor = 10
sl = [slice(None,None,None), slice(None,None,None)]

# Simulation snapnum 27 (z = 0.1), xy box size: 100Mpc, z slice width: 5Mpc
files_SF_27 = ['/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5__fromSFR.npz']

files_noSF_27 = ['/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_27_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5_noSFR.npz']


# Simulation snapnum 28 (z = 0), xy box size: 100Mpc, z slice width: 5Mpc,
files_SF_28 = ['/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5__fromSFR.npz',
               '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5__fromSFR.npz']

files_noSF_28 = ['/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen12.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen17.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen2.5_noSFR.npz',
                 '/Users/lokhorst/Eagle/emission_halpha_L0100N1504_28_test2_SmAb_C2Sm_32000pix_5.000000slice_zcen7.5_noSFR.npz']

# Read in data from snapsnot 27
print('Reading: '+files_noSF_27[0]+'...')
data1 = (np.load(files_noSF_27[0])['arr_0'])[sl]
data1 = get_halpha_SB.imreduce(data1, factor, log=True, method = 'average')

print('Reading: '+files_noSF_28[0]+'...')
data2 = (np.load(files_noSF_28[0])['arr_0'])[sl]
data2 = get_halpha_SB.imreduce(data2, factor, log=True, method = 'average')

print('Reading" "'+files_SF_27[0]+'...')
data11 = (np.load(files_SF_27[0])['arr_0'])[sl]
data11 = get_halpha_SB.imreduce(data11, factor, log=True, method = 'average')
print('5 Mpc slice...')
data_27 = np.log10(10**data1+10**data11)

print('Reading" "'+files_SF_28[0]+'...')
data22 = (np.load(files_SF_28[0])['arr_0'])[sl]
data22 = get_halpha_SB.imreduce(data22, factor, log=True, method = 'average')
print('5 Mpc slice...')
data_28 = np.log10(10**data2+10**data22)

print('delete %s, %s...'%(files_noSF_27[0],files_noSF_28[0]))
del data11
del data22

fig = plt.subplots(1,1,figsize=[6,5])
plt.hist(data1.flatten(),bins=50,log='True',normed='True',histtype='step',label='noSF 27')
plt.hist(data2.flatten(),bins=50,log='True',normed='True',histtype='step',label='noSF 28')
plt.hist(data_27.flatten(),bins=50,log='True',normed='True',histtype='step',label='with SF 27')
plt.hist(data_28.flatten(),bins=50,log='True',normed='True',histtype='step',label='with SF 28')
plt.ylabel(r'$N^{-1}_{pix}dN_{pix}/dlog_{10}S_{B}$')
plt.xlabel(r'$log_{10}S_{B}$ ($ph$ $s^{-1}$ $cm^{-2}$ $sr^{-1})$')
plt.legend()
plt.savefig('../plots/Eagle_SBhist.pdf')
plt.show()
