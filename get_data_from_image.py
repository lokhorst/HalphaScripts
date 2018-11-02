#!/usr/bin/env python

""" get_data_from_image.py -- makes all the white transparent in the image, then grabs all the black 
lines to a list and can save them to a file so that we have the values of the plots (see FGplot.py 
on how to recover the plot from the saved list).

Usage: get_data_from_image [-v] [-o outputimage] [-f filename] <inputimage> 

Options:
    -v, --verbose                               Show extra information [default: False]
    
    -o OUTPUTIMAGE, --outimage OUTPUTIMAGE      Output image name for transparent version of image [default: transparent.png]
    -f FILENAME, --outfile FILENAME             Output file name to contain coordinates of data [default: datafromimage.txt]

Examples:
    
Author: lokhorst

"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
import docopt
from PIL import Image

if __name__ == "__main__":
    
    arguments = docopt.docopt(__doc__)
    
    # Mandatory arguments
    imgname = arguments['<inputimage>']
    
    # Non-mandatory arguments
    verbose = arguments['--verbose']
    outimgname = arguments['--outimage']
    outfilename = arguments['--outfile']
    
    if verbose:
        print arguments
    
#    img = Image.open('original.png')
    img = Image.open(imgname)
    img = img.convert("RGBA")

    pixdata = img.load()

    xvalue = []
    yvalue = []
    if verbose:
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

#    img.save("original_trans.png", "PNG")
    img.save(outimgname, "PNG")

    ## Uncomment to save data in lists to files
    with open('xdata_tmp.dat', 'wb') as f:
        pickle.dump(xvalue, f)
    with open('ydata_tmp.dat', 'wb') as f:
        pickle.dump(yvalue, f)
