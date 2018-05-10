import matplotlib.pyplot as plt
import numpy as np
im1 = plt.imread('one3.png')    # with resonant scattering
im2 = plt.imread('two4.png')    # without resonant scattering

plt.imshow(im1, zorder=1)
plt.plot(0,0,'k',label='with resonant scattering')
plt.imshow(im2, zorder=0)
plt.plot(0,0,'r',label='without resonant scattering')

plt.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    left = 'off',
    right = 'off')
   # labelbottom='off') # labels along the bottom edge are off
#plt.locator_params(axis='y',nbins=5)#to specify number of ticks on both or any single axes
ylabels = ['0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35']
ylabels.reverse()
plt.yticks(np.arange(0,364.+364./7,364./7),ylabels)
xlabels = ['0','10','20','30','40','50']
plt.xticks(np.arange(0,499.+499./5,499./5),xlabels)
plt.xlabel('R (proper kpc)')
plt.ylabel('$10^{18}$ x I(R)')
plt.ylim(364,0)
plt.xlim(0,499)
plt.legend()
plt.show()