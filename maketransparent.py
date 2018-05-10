# This function makes all the white transparent
# This function grabs all the black lines to a list and can save them to a file so that we have the values of the plots (see FGplot.py on how to recover the plot from the saved list).

import matplotlib.pyplot as plt
import numpy as np
import pickle
from PIL import Image

img = Image.open('original.png')
img = img.convert("RGBA")

pixdata = img.load()

xvalue = []
yvalue = []
print(img.size[1])
print(img.size[0])

for y in xrange(img.size[1]):
    for x in xrange(img.size[0]):
        if pixdata[x, y] == (255, 255, 254, 255): # used to be (255, 255, 255, 255) # white
            pixdata[x, y] = (255, 255, 254, 0)    # used to be (255, 255, 255, 0)
        if pixdata[x, y] == (27,25,25, 255):   # black
            pixdata[x, y] = (255,0,0,255)
            xvalue.append(x)
            yvalue.append(y)

img.save("original_trans.png", "PNG")

## Uncomment to save data in lists to files
with open('xdata_tmp.dat', 'wb') as f:
    pickle.dump(xvalue, f)
with open('ydata_tmp.dat', 'wb') as f:
    pickle.dump(yvalue, f)
