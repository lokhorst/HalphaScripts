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

# Simulation snapnum 27 (z = 0.1), xy box size: 100Mpc, z slice width: 5Mpc, zloc: 12.5Mpc
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


# The following is snapnum 27 data (will do snapnum 28 seperately -- can compare the two then)p
print('Load non-star-forming Halpha emission data...')
print('data1 ('+files_noSF_27[0]+')...')
data1 = (np.load(files_noSF_27[0])['arr_0'])[sl]
#print('data11 ('+files_SF_27[0]+')...')
#data11 = (np.load(files_SF_27[0])['arr_0'])[sl]

data1 = get_halpha_SB.imreduce(data1, factor, log=True, method = 'average')
#data11 = get_halpha_SB.imreduce(data11, factor, log=True, method = 'average')

print('5 Mpc slice...')
data_5 = np.log10(10**data1)#+10**data11)

print('delete data1, data11...')
del data1
#del data11

print('data2 ('+files_noSF_27[1]+')...')
data2 = (np.load(files_noSF_27[1])['arr_0'])[sl]
#print('data22 ('+files_SF_27[1]+')...')
#data22 = (np.load(files_SF_27[1])['arr_0'])[sl]

data2 = get_halpha_SB.imreduce(data2, factor, log=True, method = 'average')
#data22 = get_halpha_SB.imreduce(data22, factor, log=True, method = 'average')

print('10 Mpc slice...')
data_10 = np.log10(10**data_5+10**data2)#+10**data22)

print('delete data2, data22..')
del data2
#del data22

print('data3 ('+files_noSF_27[2]+')...')
data3 = (np.load(files_noSF_27[2])['arr_0'])[sl]
data3 = get_halpha_SB.imreduce(data3, factor, log=True, method = 'average')

print('15 Mpc slice...')
data_15 = np.log10(10**data_10+10**data3)#+10**data11+10**data22+10**data33)

del data3
#print('data33 ('+files_SF_27[2]+')...')
#data33 = (np.load(files_SF_27[2])['arr_0'])[sl]

print('data4 ('+files_noSF_27[3]+')...')
data4 = (np.load(files_noSF_27[3])['arr_0'])[sl]
data4 = get_halpha_SB.imreduce(data4, factor, log=True, method = 'average')

print('20 Mpc slice...')
data_20 = np.log10(10**data_15+10**data4)#+10**data11+10**data22+10**data33+10**data44)

del data4

#print('Load star-forming Halpha emission data')
#print('data44 ('+files_SF_27[3]+')...')
#data44 = (np.load(files_SF_27[3])['arr_0'])[sl]

#print('Reduce all the data that we just loaded...')
#data33 = get_halpha_SB.imreduce(data33, factor, log=True, method = 'average')
#data44 = get_halpha_SB.imreduce(data44, factor, log=True, method = 'average')

#print('Add slices together to create larger slices (also adding together the SFing and non-SFing data)..')



print('Plot an example region (20 Mpc box)...')
# Plotting parameters
xystarts = [40.,0.]
size     = 20.

fig = plt.figure(figsize = (16.5, 15.))
ax1 = plt.subplot(111)
get_halpha_SB.makemap(data_15[(40./100.*3200.):(60./100.*3200.),0:(20./100.*3200.)],size,ax1,xystarts = xystarts)

# Pick out regions along the filaments in the map
ax1.plot([53,53,56,56,53],[9.2,10,8.5,7.7,9.2],color='r', label='Region 3')
ax1.plot(np.array([46.2,47.2,48.2,47.2,46.2]),[14,14,10,10,14],color='r', label='Region 1')
ax1.plot(np.array([43,43,46,46,43]),[7.5,8.2,7.2,6.5,7.5],color='r', label = 'Region 2')

plt.savefig('HalphaSBregions_15Mpcmap.pdf')

# Pick out the data inside the regions

xbox_3 = np.array([53,53,56,56])*3200./100.
ybox_3 = np.array([9.2,10,8.5,7.7])*3200./100.

xbox_2 = np.array([43,43,46,46])*3200./100.
ybox_2 = (np.array([7.8,8.,7.1,6.8])-0.05)*3200./100.

xbox_1 = np.array([47.4,46.2,46.9,48.1])*3200./100.
ybox_1 = np.array([10.5,14,14,10.5])*3200./100.

xbox = xbox_1
ybox = ybox_1
print('Region: '+str(xbox)+' , '+str(ybox))
xfull, yfull= get_halpha_SB.indices_region(xbox.astype(int),ybox.astype(int))

# Cycle through the x and y arrays, and put the SB data in the region into new arrays
SBdata_5 = np.zeros(xfull.shape)
SBdata_10 = np.zeros(xfull.shape)
SBdata_15 = np.zeros(xfull.shape)
SBdata_20 = np.zeros(xfull.shape)
for i in range(yfull.shape[0]):
    for j in range(yfull.shape[1]):
        SBdata_5[i,j]  = data_5[xfull[i,j],yfull[i,j]]
        SBdata_10[i,j] = data_10[xfull[i,j],yfull[i,j]]
        SBdata_15[i,j] = data_15[xfull[i,j],yfull[i,j]]
        SBdata_20[i,j] = data_20[xfull[i,j],yfull[i,j]]

SB_list = ['SBdata_5','SBdata_10','SBdata_15','SBdata_20']
SBdata_average = np.array([np.log10(np.mean(10**SBdata_5)),np.log10(np.mean(10**SBdata_10)),np.log10(np.mean(10**SBdata_15)),np.log10(np.mean(10**SBdata_20))])
SBdata_median  = np.array([np.median(SBdata_5),np.median(SBdata_10),np.median(SBdata_15),np.median(SBdata_20)])
print(SB_list)
print('Average SB in the selected region: '+str(SBdata_average))
print('Median SB in the selected region:'+str(SBdata_median))


#['SBdata_5', 'SBdata_10', 'SBdata_15', 'SBdata_20']
#Average SB in the selected region: [ 2.63944308  2.63944324  2.63946358  2.66657444]
#Median SB in the selected region:[-1.37115943 -1.36929631 -1.36385357 -1.26499414]

#print(SBdata)
#print(data_noSF[(46.2/100.*3200.):(48.1/100.*3200.),(10.5/100.*3200.):(14./100.*3200.)])

#f = open("SBdata.txt", "w")
#for i in SBdata:
#    f.write( str(i) +'\n' )      # str() converts to string
#f.close()

#f = open("SBdata_original.txt", "w")
#for i in data_noSF[(46.2/100.*3200.):(48.1/100.*3200.),(10.5/100.*3200.):(14./100.*3200.)]:
#    f.write( str(i) +'\n' )      # str() converts to string
#f.close()

# If want to also use the Eagle galaxy catalogue (such as pointing out the SFing regions, and how SFing they are)
"""
import eagleSqlTools as sql

mySim = ('RefL0100N1504',100.)
con   = sql.connect("dlokhorst",password="mxdPB54Y")

myQuery  = "SELECT \
                SH.GalaxyID, \
                SH.CentreOfMass_x, \
                SH.CentreOfMass_y, \
                SH.CentreOfMass_z, \
                SH.StarFormationRate as SFR \
            FROM \
                %s_SubHalo as SH \
            WHERE \
                SH.SnapNum = 27 and \
                SH.StarFormationRate > 0.00001 "%('RefL0100N1504')

myData = sql.execute_query(con,myQuery)

x = myData['CentreOfMass_x'][:]
y = myData['CentreOfMass_y'][:]
z = myData['CentreOfMass_z'][:]
SFR = myData['SFR'][:]
x_inbox = x[(x < 60.) & (x > 40.) & (y < 20.) & (y > 0.) & (z < 15.) & (z > 10.)]
y_inbox = y[(x < 60.) & (x > 40.) & (y < 20.) & (y > 0.) & (z < 15.) & (z > 10.)]
z_inbox = z[(x < 60.) & (x > 40.) & (y < 20.) & (y > 0.) & (z < 15.) & (z > 10.)]
SFR_inbox = SFR[(x < 60.) & (x > 40.) & (y < 20.) & (y > 0.) & (z < 15.) & (z > 10.)]

print(len(SFR_inbox[SFR_inbox<1.]))
print(len(SFR_inbox[SFR_inbox<0.1]))
print(len(SFR_inbox[SFR_inbox<0.01]))

# Circle SFing regions
for i in range(len(x_inbox)):
    if SFR_inbox[i] > 1.0:
        color = 'r'
        zorder = 5
    elif SFR_inbox[i] > 0.1:
        color = 'y'
        zorder = 4
    elif SFR_inbox[i] > 0.01:
        color = 'g'
        zorder = 3
    else:
        color = 'b'
        zorder = 2
    circle1 = plt.Circle((x_inbox[i],y_inbox[i]), 0.5, color=color,fill=False,zorder=zorder)
    ax2.add_artist(circle1)
"""