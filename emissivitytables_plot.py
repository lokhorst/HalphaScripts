"""
plotting the emissivity tables

Emissivity Tables are saved on scinet.  Pieced together how to use the hdf5 from make_maps_perrgn and make_maps_opt_locs.

Halpha lineind is 9 (from make_maps_opt_locs)
Lyalpha linind is 1 (from make_maps_opt_locs)

In [41]: tablefile.keys()
Out[41]: [u'header', u'lambda', u'lines', u'logd', u'logt']

In [53]: tablefile.values()
Out[53]: 
[<HDF5 group "/header" (4 members)>,
 <HDF5 dataset "lambda": shape (35,), type "|S12">,
 <HDF5 dataset "lines": shape (131, 46, 35), type "<f8">,
 <HDF5 dataset "logd": shape (46,), type "<f4">,
 <HDF5 dataset "logt": shape (131,), type "<f4">]

logrhocm3 == lognHcm3 (presumed from make_maps_perrgn)

In [44]: lambdas
Out[44]: 
array(['H  1  1026A', 'H  1  1216A', 'H  1  3798A', 'H  1  3835A',
       'H  1  3889A', 'H  1  3970A', 'H  1  4102A', 'H  1  4340A',
       'H  1  4861A', 'H  1  6563A', 'H  1 1.094m', 'H  1 1.282m',
       'H  1 1.875m', 'H  1 4.051m', 'H  1 913.3A', 'H  1 913.5A',
       'H  1 913.6A', 'H  1 913.8A', 'H  1 914.0A', 'H  1 914.3A',
       'H  1 914.6A', 'H  1 914.9A', 'H  1 915.3A', 'H  1 915.8A',
       'H  1 916.4A', 'H  1 917.2A', 'H  1 918.1A', 'H  1 919.4A',
       'H  1 921.0A', 'H  1 923.2A', 'H  1 926.2A', 'H  1 930.8A',
       'H  1 937.8A', 'H  1 949.7A', 'H  1 972.5A'], 
      dtype='|S12')

In [5]: wavelen
Out[5]: 
array(['O  1  5577A', 'O  1  6300A', 'O  1  6363A', 'O  1 145.5m',
       'O  1 63.17m', 'O  2 355.6A', 'O  2 356.9A', 'O  2 358.8A',
       'O  2 360.6A', 'O  2 361.8A', 'O  2 366.1A', 'O  2 373.6A',
       'O  2 384.4A', 'O  2 388.4A', 'O  2 414.6A', 'O  2 425.2A',
       'O  2 530.8A', 'O  2 833.8A', 'O  3  1661A', 'O  3  1666A',
       'O  3  2321A', 'O  3  4959A', 'O  3  5007A', 'O  3 206.4A',
       'O  3 215.3A', 'O  3 230.2A', 'O  3 231.1A', 'O  3 232.4A',
       'O  3 233.3A', 'O  3 234.2A', 'O  3 235.5A', 'O  3 236.9A',
       'O  3 238.9A', 'O  3 241.1A', 'O  3 241.2A', 'O  3 244.7A',
       'O  3 248.6A', 'O  3 248.7A', 'O  3 255.4A', 'O  3 263.2A',
       'O  3 263.9A', 'O  3 265.0A', 'O  3 268.1A', 'O  3 276.9A',
       'O  3 280.4A', 'O  3 299.3A', 'O  3 300.5A', 'O  3 305.4A',
       'O  3 373.6A', 'O  3 374.1A', 'O  3 500.7A', 'O  3 51.80m',
       'O  3 698.2A', 'O  3 703.4A', 'O  3 835.0A', 'O  3 88.33m',
       'O  4  1397A', 'O  4  1400A', 'O  4  1401A', 'O  4  1405A',
       'O  4  1407A', 'O  4 145.8A', 'O  4 146.8A', 'O  4 148.1A',
       'O  4 150.1A', 'O  4 152.4A', 'O  4 153.6A', 'O  4 158.6A',
       'O  4 159.7A', 'O  4 165.4A', 'O  4 166.6A', 'O  4 167.1A',
       'O  4 168.1A', 'O  4 168.2A', 'O  4 168.9A', 'O  4 169.9A',
       'O  4 170.7A', 'O  4 171.6A', 'O  4 174.3A', 'O  4 174.7A',
       'O  4 180.3A', 'O  4 181.8A', 'O  4 182.7A', 'O  4 184.2A',
       'O  4 186.0A', 'O  4 188.2A', 'O  4 196.4A', 'O  4 196.6A',
       'O  4 200.8A', 'O  4 203.0A', 'O  4 206.6A', 'O  4 207.6A',
       'O  4 214.8A', 'O  4 239.1A', 'O  4 25.88m', 'O  4 250.3A',
       'O  4 252.5A', 'O  4 260.5A', 'O  4 265.5A', 'O  4 273.5A',
       'O  4 280.5A', 'O  4 549.3A', 'O  4 554.4A', 'O  4 609.4A',
       'O  4 617.0A', 'O  4 779.9A', 'O  4 789.0A', 'O  5  1218A',
       'O  5 112.6A', 'O  5 113.4A', 'O  5 114.6A', 'O  5 116.4A',
       'O  5 118.2A', 'O  5 119.4A', 'O  5 124.9A', 'O  5 135.8A',
       'O  5 139.1A', 'O  5 172.6A', 'O  5 630.0A', 'O  6  1032A',
       'O  6  1038A', 'O  6 105.0A', 'O  6 116.0A', 'O  6 150.0A',
       'O  6 93.24A', 'O  6 94.06A', 'O  6 95.23A', 'O  6 97.00A',
       'O  6 99.86A', 'O  7  1624A', 'O  7  1640A', 'O  7 120.3A',
       'O  7 128.5A', 'O  7 135.8A', 'O  7 17.77A', 'O  7 18.63A',
       'O  7 21.60A', 'O  7 21.80A', 'O  7 21.81A', 'O  7 22.10A',
       'O  7 382.0A', 'O  7 96.15A', 'O  8 102.5A', 'O  8 14.46A',
       'O  8 14.53A', 'O  8 14.64A', 'O  8 14.82A', 'O  8 15.18A',
       'O  8 16.01A', 'O  8 18.97A', 'O  8 75.89A'], 
      dtype='|S12')



author:lokhorst
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

tablefilename = '/Users/lokhorst/Documents/Eagle/scripts/EmmisivityTables/hydrogen.hdf5'
tablefile = h5py.File(tablefilename, "r")
logTK =   np.array(tablefile.get('logt'),dtype=np.float32)
logrhocm3 =   np.array(tablefile.get('logd'),dtype=np.float32)
lines =   np.array(tablefile.get('lines'),dtype=np.float32)
wavelen =   np.array(tablefile.get('lambda'))
tablefile.close()

Halpha_emtable = lines[:,:,9] 
Lyalpha_emtable = lines[:,:,1]
Halpha_nH_3 = Halpha_emtable[:,25]
Lyalpha_nH_3 = Lyalpha_emtable[:,25]

plt.plot(logTK,Halpha_nH_3,label=r'H$\alpha$')
plt.plot(logTK,Lyalpha_nH_3,label=r'Ly$\alpha$')
plt.title(r'Emissivity (log$_{10}$ n$_H$ = -3, solar abundance)')
plt.xlim(3,8)
plt.xlabel(r'Log$_{10}$ Temperature (K)')
plt.ylabel(r'Log$_{10}$ $\epsilon / n_H^2$ ($erg$ $s^{-1}$ $cm^3$)')
#plt.legend(fancybox = True, loc=1)
#plt.savefig('emissivitytable.pdf')
#plt.show()

# oxygen lines in the visible: 0,1,2, 21, 22, 
tablefilename = '/Users/lokhorst/Documents/Eagle/scripts/EmmisivityTables/oxygen.hdf5'
tablefile = h5py.File(tablefilename, "r")
logTK =   np.array(tablefile.get('logt'),dtype=np.float32)
logrhocm3 =   np.array(tablefile.get('logd'),dtype=np.float32)
lines =   np.array(tablefile.get('lines'),dtype=np.float32)
wavelen =   np.array(tablefile.get('lambda'))
tablefile.close()

### Solar abundances here -- need to take that into account...  ###

indices = [0,1,2,21,22]

for index in indices:
    emtable = lines[:,:,index]
    emtable_nH_3 = emtable[:,25]
    plt.plot(logTK,emtable_nH_3,label=wavelen[index])
    
plt.legend(fancybox=True,loc=1)
plt.ylim(-30,-20)
#plt.xlim(3,8)
#plt.xlabel(r'Log$_{10}$ Temperature (K)')
#plt.ylabel(r'Log$_{10}$ $\epsilon / n_H^2$ ($erg$ $s^{-1}$ $cm^3$)')

plt.savefig('emissivitytable.pdf')

