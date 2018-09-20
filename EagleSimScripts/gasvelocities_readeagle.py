#!/usr/bin/env python

"""gasvelocities_readeagle.py -- Uses the read_eagle routine to read in an EAGLE simulation snapshot, select out a region, then get the velocities of the gas particles within the region.

Usage:
    gasvelocities_readeagle [-h] [-v] [-o DIRECTORY]

Options:
    -h, --help                                  Show this screen
    -v, --verbose                               Show extra information [default: False]

    -o DIRECTORY, --outputdir DIRECTORY         Output directory name  [default: .]
    

Examples:
    python gasvelocities_readeagle.py -v ~/eagle/L0012N0188/REFERENCE/data/snapshot_028_z000p000/

Notes:
In order for this script to run, need to have read_eagle installed. (https://github.com/jchelly/read_eagle)
Also need to have the read_header script in the same directory. (http://icc.dur.ac.uk/Eagle/Scripts/read_header.py)

"""

import docopt

import os.path
import math
import sys
import itertools
import subprocess
import re

import numpy as np
import matplotlib.pyplot as plt 
from read_eagle import EagleSnapshot
from read_header import read_header
import h5py

import astropy.units as u

def print_verbose_string(printme):
    print >> sys.stderr, "VERBOSE: %s" % printme
    
class GasVelocities_ReadEagle:
    
    def __init__(self, gn, sgn, centre, velocity, load_region_length=2):

        # Load information from the header.
        self.a, self.h, self.boxsize = read_header()

        # Load data.
        self.gas = self.read_galaxy(0, gn, sgn, centre, velocity, load_region_length)

        # Plot.
        self.plot(gn, sgn)
        self.hist(gn, sgn)
        self.locplot(gn, sgn, centre)

    def read_galaxy(self, itype, gn, sgn, centre, velocity, load_region_length):
        """ For a given galaxy (defined by its GroupNumber and SubGroupNumber)
        extract the Temperature, Density and StarFormationRate of all gas particles
        using the read_eagle routine. Conversion factors are still loaded directly
        from the hdf5 files. """

        data = {}

        # Initialize read_eagle module.
        eagle_data = EagleSnapshot('./data/snap_028_z000p000.0.hdf5')

        # Put centre into cMpc/h units.
        centre *= self.h

        # Select region to load, a 'load_region_length' cMpc/h cube centred on 'centre'.
        region = np.array([
            (centre[0]-0.5*load_region_length), (centre[0]+0.5*load_region_length),
            (centre[1]-0.5*load_region_length), (centre[1]+0.5*load_region_length),
            (centre[2]-0.5*load_region_length), (centre[2]+0.5*load_region_length)
        ]) 
        eagle_data.select_region(*region)

        # Load data using read_eagle, load conversion factors manually.
        f = h5py.File('./data/snap_028_z000p000.0.hdf5', 'r')
        for att in ['GroupNumber', 'SubGroupNumber', 'StarFormationRate', 'Temperature', 'Density', 'Coordinates', 'Velocity']:
            tmp  = eagle_data.read_dataset(itype, att)
            cgs  = f['PartType%i/%s'%(itype, att)].attrs.get('CGSConversionFactor')
            aexp = f['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
            hexp = f['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')
            data[att] = np.multiply(tmp, cgs * self.a**aexp * self.h**hexp, dtype='f8')  # This puts all attributes into their proper values (not /h)
        f.close()

        # Mask to selected GroupNumber and SubGroupNumber.
        mask = np.logical_and(data['GroupNumber'] == gn, data['SubGroupNumber'] == sgn)
        for att in data.keys():
            data[att] = data[att][mask]
        
        # Put centre back into proper units
        centre /= self.h
        
        centre_cgs = centre * u.Mpc.to(u.cm)
        velocity_cgs = velocity * u.km.to(u.cm)
        
        data['rel_location'] = np.asarray([(partcoords-centre_cgs) for partcoords in data['Coordinates']])
        data['radius'] = np.asarray([float(np.sqrt(np.dot(relloc,relloc))) for relloc in data['rel_location']])
        
        data['rel_velocity'] = np.asarray([(partvelocity-velocity_cgs) for partvelocity in data['Velocity']])
        data['veldot'] = np.asarray([np.dot(relvel,relloc)/np.sqrt(np.dot(relloc,relloc)) for relvel, relloc in zip(data['rel_velocity'],data['rel_location'])])
        data['veldot_norm'] = np.asarray([np.dot(relvel,relloc)/(np.sqrt(np.dot(relloc,relloc))*np.sqrt(np.dot(relvel,relvel))) for relvel, relloc in zip(data['rel_velocity'],data['rel_location'])])

        for i in range(1):
            print data['Velocity'][i] * u.cm.to(u.km)
            print velocity_cgs* u.cm.to(u.km)
            print (data['Velocity'][i]-velocity_cgs) * u.cm.to(u.km)
            print ""
            print data['Coordinates'][i] * u.cm.to(u.Mpc)
            print centre_cgs * u.cm.to(u.Mpc)
            print (data['Coordinates'][i]-centre_cgs) * u.cm.to(u.Mpc)
            print ""
            print data['rel_location'][i] * u.cm.to(u.Mpc)
            print data['rel_velocity'][i] * u.cm.to(u.km)
            print ""
            print np.dot((data['Velocity'][i]-velocity_cgs),(data['Coordinates'][i]-centre_cgs))
            print np.dot(data['rel_velocity'][i],data['rel_location'][i])
            print np.sqrt(np.dot(data['rel_velocity'][i],data['rel_velocity'][i]))
            print np.sqrt(np.dot(data['rel_location'][i],data['rel_location'][i]))
            print ""
            print data['veldot'][i] * u.cm.to(u.km)
            print np.sqrt(np.dot(data['rel_velocity'][i],data['rel_velocity'][i])) * u.cm.to(u.km)
            print data['veldot_norm'][i]
            print ""
        print len(data['GroupNumber'])
        
        return data
        
    def plot(self, gn, sgn):
        """ Plot Relative Velocity -- Radius relation. (weight by mass? temperature?) """
        plt.figure()        
        
        ### would it be straightforward to add the emission calculation here, and then sort by H-alpha emission luminosity as well?  not really -- need to use a c function to interpolate the emissivity table -- doable but maybe save it for later if this stuff isn't clear (if split by temp, density etc)

        # Plot currently star forming gas red.
        mask = np.where(self.gas['StarFormationRate'] > 0)
        plt.scatter(np.log10(self.gas['Density'][mask]), np.log10(self.gas['Temperature'][mask]),
            c='red', s=3, edgecolor='none')

        # Plot currently non star forming gas blue.
        mask = np.where(self.gas['StarFormationRate'] == 0)
        plt.scatter(np.log10(self.gas['Density'][mask]), np.log10(self.gas['Temperature'][mask]),
            c='blue', s=3, edgecolor='none')

        # Save plot.
        plt.minorticks_on()
        plt.ylabel('log10 Temperature [K]'); plt.xlabel('log10 Density [g/cm**3]')
        plt.tight_layout()
        plt.savefig('gasvelocities_readeagle_gn%s_sgn%s.png'%(gn,sgn))
        plt.close()

    def locplot(self, gn, sgn, centre):
        """ Plot locations of the particles in x, y """
        plt.figure()    
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize = (7,7))        
            
        mask = np.where(self.gas['Temperature'] > 10**6)
        ax1.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc), color = 'red', alpha = 0.005,label = 'logT > 6')        
#        density,xedges,yedges = np.histogram2d(self.gas['Coordinates'][mask][:,0],self.gas['Coordinates'][mask][:,1], bins=60)
#        density = density.T
#        ax1.imshow(density,interpolation='nearest',origin='low',extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],cmap='Greys')
        
        mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5))
        ax2.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc), color = 'orange', alpha=0.3, label = '5 < logT < 6')
        
        mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4))
        ax3.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc), color = 'green', alpha=0.3, label = '4 < logT < 5')
        
        mask = np.where(self.gas['Temperature'] <= 10**4)
        ax4.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc), color = 'blue', alpha=0.3, label = 'logT < 4')
        
        for ax in [ax1,ax2,ax3,ax4]:
         #   ax.plot(centre[0]* u.Mpc.to(u.cm),centre[1]* u.Mpc.to(u.cm),'o',color='black')
            ax.plot(centre[0],centre[1],'o',color='black')
            ax.set_ylabel('Y (cMpc)'); ax.set_xlabel('X (cMpc)')
            ax.minorticks_on()
            
        for ax in [ax2,ax3,ax4]:
            ax.legend()
        leg = ax1.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(0.3)
            
        plt.tight_layout()
        plt.savefig('gasvelocities_locations_gn%s_sgn%s.png'%(gn,sgn))
        plt.close()
        

    def hist(self, gn, sgn):
        """ Plot Histograms of the veldot == costheta value  """
        fig, (ax2, ax1) = plt.subplots(2,1,figsize = (7,7))        
       # fig, (ax1) = plt.subplots(1,1,figsize = (5,5))        

        print self.gas['veldot_norm'] 

        for ax in [ax1,ax2]:
            mask = np.where(self.gas['Temperature'] > 10**6)
            ax.hist(self.gas['veldot_norm'][mask], color = 'red', alpha = 0.5, label = 'logT > 6')
        
            mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5))
            ax.hist(self.gas['veldot_norm'][mask], color = 'orange', alpha = 0.5, label = '5 < logT < 6')
        
            mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4))
            ax.hist(self.gas['veldot_norm'][mask], color = 'green', alpha = 0.5, label = '4 < logT < 5')
        
            mask = np.where(self.gas['Temperature'] <= 10**4)
            ax.hist(self.gas['veldot_norm'][mask], color = 'blue', alpha = 0.8, label = 'logT < 4')

        ax2.legend(framealpha=0.5, loc = 4)

#        # Save plot.
        ax1.minorticks_on()
        ax1.set_ylim(0,500)
        ax2.set_ylim(2000,9000)
        
        ax2.set_ylabel('number')
        ax1.set_ylabel('number'); ax1.set_xlabel('cos theta between relative velocity and relative location of gas and galaxy')
#        plt.ylabel('log10 Temperature [K]'); plt.xlabel('log10 Density [g/cm**3]')
        plt.tight_layout()
        
        ax2.set_xticklabels([])
        ax2.set_xlabel('')
        fig.subplots_adjust( hspace = 0.07)
        
        plt.savefig('gasvelocities_histogram_gn%s_sgn%s.png'%(gn,sgn))
        plt.close()
        

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":
    
    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    #snapshotdirectory = arguments['<snapshotdirectory>']

    # Non-mandatory options without arguments
    verbose = arguments['--verbose']
    
    # Non-mandatory options with arguments
    output_dir = arguments['--outputdir']
    
    # Centre is the COP for GN=1 SGN=0 taken from the database.  ### THIS IS FOR THE 12Mpc REFERENCE BOX ###
    centre = np.array([12.08808994,4.47437191,1.41333473])  # cMpc
    velocity = np.array([25.2508, 25.5319, 2.74462]) # km / s
    x = GasVelocities_ReadEagle(1, 0, centre,velocity, load_region_length=0.2)
    
    # Use the COP for the galaxies that we consider in the paper (print off COP coords in extract_FOV_and_cutout_galaxies)
    # 4067, 0, 50.487, 13.5846, 12.3391
    #centre = np.array([50.487, 13.5846, 12.3391])  # cMpc
    #x = GasVelocities_ReadEagle(4067, 0, centre)
    