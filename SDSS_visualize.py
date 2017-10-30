import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import Distance
from astropy import units as u
from astropy.coordinates import SkyCoord
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits

specObjID=[]
z=[]
ra=[]
dec=[]
cx=[]
cy=[]
cz=[]
with open('../SDSS/Skyserver_SQL_1_xyz.csv','r') as f:
    for line in f:
        if line[0]!='#' and len(line)>5:
            tmp = line.split(',')
            specObjID.append(float(tmp[0]))
            z.append(float(tmp[1]))
            ra.append(float(tmp[2]))
            dec.append(float(tmp[3]))
            cx.append(float(tmp[4]))
            cy.append(float(tmp[5]))
            cz.append(float(tmp[6]))

with open('../SDSS/Skyserver_SQL_2_xyz.csv','r') as f:
    for line in f:
        if line[0]!='#' and len(line)>5:
            tmp = line.split(',')
            specObjID.append(float(tmp[0]))
            z.append(float(tmp[1]))
            ra.append(float(tmp[2]))
            dec.append(float(tmp[3]))
            cx.append(float(tmp[4]))
            cy.append(float(tmp[5]))
            cz.append(float(tmp[6]))
            
specObjID = np.array(specObjID,'d')
z = np.array(z,'d')
ra = np.array(ra,'d')
dec = np.array(dec,'d')
cx = np.array(cx)
cy = np.array(cy)
cz = np.array(cz)

c = SkyCoord(ra=ra[0]*u.degree, dec=dec[0]*u.degree, distance=z[0]*300000/70*u.Mpc)
print(c.cartesian.x,c.cartesian.y,c.cartesian.z)
print(cx[0],cy[0],cz[0])
normnum = np.sqrt(cx[0]**2+cy[0]**2+cz[0]**2)
## The result is that these are consistent

coord_array = SkyCoord(ra=ra*u.degree,dec=dec*u.degree, distance=z*300000/70*u.Mpc)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
coord_100 = coord_array[coord_array.distance.value<100.]
ax.scatter(coord_100.cartesian.x,coord_100.cartesian.y,coord_100.cartesian.z)#,c = (256. - RVel_val/max(RVel_val)*256.) )

# numbers at specific redshifts (since we will have specific redshift bins)
filterwidth = 1. #nm
z_filterwidth = filterwidth/660. # 1nm
numbins = 0.2/z_filterwidth
plt.hist(z,bins=round(numbins))
plt.xlabel('redshift')
plt.ylabel(r'number of SDSS galaxies (per %s nm filterwidth)'%(filterwidth))
#plt.show()
plt.ylim(0,12000)
plt.xlim(0,0.2)
z_50 = 0.0115
z_100 = 0.0235
z_200 = 0.047
z_500 = 0.12
plt.plot([z_50,z_50],[0,12000],'k--')
plt.plot([z_100,z_100],[0,12000],'k--')
plt.plot([z_200,z_200],[0,12000],'k--')
plt.plot([z_500,z_500],[0,12000],'k--')
plt.text(z_50-0.006,11000,'50 Mpc',rotation=90)
plt.text(z_100-0.006,11000,'100 Mpc',rotation=90)
plt.text(z_200-0.006,11000,'200 Mpc',rotation=90)
plt.text(z_500-0.006,11000,'500 Mpc',rotation=90)

plt.savefig('SDSS_galaxies.pdf')

"""
>>> from astropy.coordinates import Distance
>>> c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, distance=770*u.kpc)
>>> c.cartesian.x  
<Quantity 568.7128654235232 kpc>
>>> c.cartesian.y  
<Quantity 107.3008974042025 kpc>
>>> c.cartesian.z  
<Quantity 507.88994291875713 kpc
"""