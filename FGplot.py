import numpy as np
import matplotlib.pyplot as plt
import pickle

# read in all the data
with open('xdata.dat', 'rb') as f:
    xvalue_two = pickle.load(f)
    
with open('ydata.dat', 'rb') as f:
    yvalue_two = pickle.load(f)
    
with open('xdata_one.dat', 'rb') as f:
    xvalue_one = pickle.load(f)
    
with open('ydata_one.dat', 'rb') as f:
    yvalue_one = pickle.load(f)
    
# check
print xvalue_two

# function to get rid of nasty values
def process(xvalue,yvalue):
    xvalue = np.array(xvalue)
    yvalue = np.array(yvalue)

    # get rid of bad x values (on the edges)
    yvalue2 = yvalue[xvalue[:]<498]
    xvalue2 = xvalue[xvalue[:]<498]
    xvalue3 = xvalue2[xvalue2[:]>0]
    yvalue3 = yvalue2[xvalue2[:]>0]
    # get rid of bad y values
    yvalue4 = yvalue3[yvalue3[:]<363]
    xvalue4 = xvalue3[yvalue3[:]<363]
    xvalue5 = xvalue4[yvalue4[:]>0]
    yvalue5 = yvalue4[yvalue4[:]>0]
    
    yvalue5 = -1.*yvalue5 + 364.
    yvalue5 = yvalue5/364. * 0.35
    xvalue5 = xvalue5/499. * 50.
    
    return [xvalue5, yvalue5]

xvalue_test1, yvalue_test1 = process(xvalue_one,yvalue_one)
xvalue_test2, yvalue_test2 = process(xvalue_two,yvalue_two)

# plot the data
plt.plot(xvalue_test1,yvalue_test1,linewidth=2.0,label='with resonant scattering')
plt.plot(xvalue_test2,yvalue_test2,linewidth=2.0,label='without resonant scattering')
plt.legend()
plt.xlabel('R (proper kpc)')
plt.ylabel('$10^{18}$ x I(R)')
plt.show()

# want to ratio the two lines
# pick a series of x-points
radius = np.arange(0,50,0.5)
# find the yvalues at each point
indices1 = []
for x in radius:
    for y in range(len(xvalue_test1)):
        if x==round(xvalue_test1[y],1):
            indices1.append(y)
            break
print(len(indices1))

indices2 = []
for x in radius:
    for y in range(len(xvalue_test2)):
        if x==round(xvalue_test2[y],1):
            indices2.append(y)
            break
print(len(indices2))

# Two subplots sharing both x/y axes
f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(8,4))
ax1.plot(xvalue_test1[indices1],yvalue_test1[indices1],'b-',label='with resonant scattering')
ax1.plot(xvalue_test2[indices2],yvalue_test2[indices2],'g-',label='without resonant scattering')
ax1.set_ylabel(r'$10^{18}$ x I$_{avg}$(R)') # (erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$)')   # circularly averaged surface brightness profile
ax1.set_yticks([0,0.1,0.2,0.3])
ax1.legend(shadow=True, fancybox=True)
ax1.get_legend()
#ax1.set_title('Sharing both axes')
ax2.plot(xvalue_test1[indices1],((yvalue_test1[indices1]-yvalue_test2[indices2])/yvalue_test2[indices2]),'k')
ax2.plot([0,50],[0,0],'k--')
ax2.set_ylim(-2,2)
ax2.set_ylabel(r'$\frac{I_{w} - I_{wo}}{I_{wo}}$',fontsize=20) # (erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$)')   # circularly averaged surface brightness profile

#ax2.set_ylabel('(I$_{w}$ - I$_{wo}$)/I$_{wo}$')
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
plt.xlabel('R (proper kpc)')
ax2.set_yticks([-2.0,-1.0,0,1.])
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.show()

