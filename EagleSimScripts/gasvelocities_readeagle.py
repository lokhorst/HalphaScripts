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
import mpl_toolkits.axes_grid1 as axgrid

def print_verbose_string(printme):
    print >> sys.stderr, "VERBOSE: %s" % printme

class GasVelocities_ReadEagle:

    def __init__(self, gn, sgn, centre, velocity, load_region_length=2, fileloc='./data/', incSF=True, useNH=False):

        # Load information from the header.
        self.a, self.h, self.boxsize = read_header(fileloc=fileloc)

        # Load data.
        self.gas = self.read_galaxy(0, gn, sgn, centre, velocity, load_region_length,fileloc)

        # Plot.
        self.plot(gn, sgn)
        self.hist(gn, sgn)
        self.locplot(gn, sgn, centre)
        self.hist_byradius(gn, sgn)
        self.hist_relvel_byradius(gn, sgn)
        self.plotOabundance(gn, sgn, incSF)
        self.plotlocsdensity(gn, sgn, incSF, useNH)

    def read_galaxy(self, itype, gn, sgn, centre, velocity, load_region_length, fileloc):
        """ For a given galaxy (defined by its GroupNumber and SubGroupNumber)
        extract the Temperature, Density and StarFormationRate of all gas particles
        using the read_eagle routine. Conversion factors are still loaded directly
        from the hdf5 files. """

        data = {}

        # Initialize read_eagle module.
        eagle_data = EagleSnapshot(fileloc+'snap_028_z000p000.0.hdf5')

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
        f = h5py.File(fileloc+'snap_028_z000p000.0.hdf5', 'r')
        for att in ['GroupNumber', 'SubGroupNumber', 'StarFormationRate', 'Temperature', 'Density', 'Coordinates', 'Velocity','SmoothedElementAbundance/Oxygen','SmoothedElementAbundance/Hydrogen']:
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

        print 'data[\'SmoothedElementAbundance/Oxygen\']:'
        print data['SmoothedElementAbundance/Oxygen']

        centre_cgs = centre * u.Mpc.to(u.cm)
        velocity_cgs = velocity * u.km.to(u.cm)

        data['rel_location'] = np.asarray([(partcoords-centre_cgs) for partcoords in data['Coordinates']])
        data['radius'] = np.asarray([float(np.sqrt(np.dot(relloc,relloc))) for relloc in data['rel_location']])

        data['rel_velocity'] = np.asarray([(partvelocity-velocity_cgs) for partvelocity in data['Velocity']])
        data['rel_speed'] = np.asarray([np.sqrt(np.dot(relvel,relvel)) for relvel in data['rel_velocity']])
        data['veldot'] = np.asarray([np.dot(relvel,relloc)/np.sqrt(np.dot(relloc,relloc)) for relvel, relloc in zip(data['rel_velocity'],data['rel_location'])])
        data['veldot_norm'] = np.asarray([np.dot(relvel,relloc)/(np.sqrt(np.dot(relloc,relloc))*np.sqrt(np.dot(relvel,relvel))) for relvel, relloc in zip(data['rel_velocity'],data['rel_location'])])

        print 'radius'
        print data['radius']* u.cm.to(u.Mpc)
        print 'speed'
        print data['rel_speed']* u.cm.to(u.km)
        print ''

        # following adapted from the make maps perrgn script
        atomw_H          = 1.00794005e0
        uthing           = 1.66053892e-24
        data['nH'] = (data['SmoothedElementAbundance/Hydrogen'])*(data['Density'])/(atomw_H*uthing)

        # DEBUG PLOT
      #  fig, (ax) = plt.subplots(1,1,figsize = (7,7))
      #  ax.hist(data['radius']* u.cm.to(u.Mpc))
      #  plt.show()
        """
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
        """

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

     #   print self.gas['veldot_norm']

        for ax in [ax1,ax2]:
            mask = np.where(self.gas['Temperature'] > 10**6)
            ax.hist(self.gas['veldot_norm'][mask], color = 'red', alpha = 0.5, label = 'logT > 6')

            mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5))
            ax.hist(self.gas['veldot_norm'][mask], color = 'orange', alpha = 0.5, label = '5 < logT < 6')

            mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4))
            ax.hist(self.gas['veldot_norm'][mask], color = 'green', alpha = 0.5, label = '4 < logT < 5')

            mask = np.where(self.gas['Temperature'] <= 10**4)
            ax.hist(self.gas['veldot_norm'][mask], color = 'blue', alpha = 0.8, label = 'logT < 4')
        hist, bin_edges = np.histogram(self.gas['veldot_norm'][mask], density=False)

        ax2.legend(framealpha=0.5, loc = 4)

#        # Save plot.
        ax1.minorticks_on()
       # ax1.set_ylim(0,500)
        ax1.set_ylim(0,max(hist)+200)
        ax2.set_ylim(2000,9000)

        ax1.set_axisbelow(True)
        ax1.grid(linestyle='--')
        ax2.grid()

        ax2.set_ylabel('number')
        ax1.set_ylabel('number'); ax1.set_xlabel('cos theta between relative velocity and relative location of gas and galaxy')
#        plt.ylabel('log10 Temperature [K]'); plt.xlabel('log10 Density [g/cm**3]')
        plt.tight_layout()

        ax2.set_xticklabels([])
        ax2.set_xlabel('')
        fig.subplots_adjust( hspace = 0.07)

        plt.savefig('gasvelocities_histogram_gn%s_sgn%s.png'%(gn,sgn))
        plt.close()

    def hist_byradius(self, gn, sgn):
        """ Plot Histograms of the veldot == costheta value  """
        fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,1,figsize = (7,7))
       # fig, (ax1) = plt.subplots(1,1,figsize = (5,5))

       # print self.gas['veldot_norm']

        for ax,minr,maxr in zip([ax1,ax2,ax3,ax4,ax5],[0,0.05,0.1,0.15,0.2],[0.05,0.1,0.15,0.2,0.4]):
            mask = np.where((self.gas['Temperature'] > 10**6) & (self.gas['radius']* u.cm.to(u.Mpc) > minr ) & (self.gas['radius']* u.cm.to(u.Mpc) <= maxr ))
            ax.hist(self.gas['veldot_norm'][mask], color = 'red', alpha = 0.5, label = 'logT > 6')

            mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5) & (self.gas['radius']* u.cm.to(u.Mpc) > minr ) & (self.gas['radius']* u.cm.to(u.Mpc) <= maxr ))
            ax.hist(self.gas['veldot_norm'][mask], color = 'orange', alpha = 0.5, label = '5 < logT < 6')

            mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4) & (self.gas['radius']* u.cm.to(u.Mpc) > minr ) & (self.gas['radius']* u.cm.to(u.Mpc) <= maxr ))
            ax.hist(self.gas['veldot_norm'][mask], color = 'green', alpha = 0.5, label = '4 < logT < 5')

            mask = np.where(self.gas['Temperature'] <= 10**4 & (self.gas['radius']* u.cm.to(u.Mpc) > minr ) & (self.gas['radius']* u.cm.to(u.Mpc) <= maxr ))
            ax.hist(self.gas['veldot_norm'][mask], color = 'blue', alpha = 0.8, label = 'logT < 4')

            ax.text(0.01,0.85,' %s Mpc < R < %s Mpc'%(minr,maxr), transform=ax.transAxes)

       # hist, bin_edges = np.histogram(self.gas['veldot_norm'][mask], density=False)

        ax1.legend(framealpha=0.5, loc = 4)

        for ax in [ax1,ax2,ax3,ax4,ax5]:
            ax.minorticks_on()
            ax.set_axisbelow(True)
            ax.grid(linestyle='--')
            ax.set_xlim(-1.1,1.1)
        for ax in [ax1,ax2,ax3,ax4]:
            ax.set_ylabel('number')
            ax.set_xticklabels([])
            ax.set_xlabel('')

        ax5.set_ylabel('number'); ax5.set_xlabel('cos theta between relative velocity and relative location of gas and galaxy')

        plt.tight_layout()
        fig.subplots_adjust( hspace = 0.07)

        plt.savefig('gasvelocities_hist_radbins_gn%s_sgn%s.png'%(gn,sgn))
        plt.close()


    def hist_relvel_byradius(self, gn, sgn):
            """ Plot Histograms of the veldot == costheta value  """
            fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,1,figsize = (7,7))
           # fig, (ax1) = plt.subplots(1,1,figsize = (5,5))

            for ax,minr,maxr in zip([ax1,ax2,ax3,ax4,ax5],[0,0.05,0.1,0.15,0.2],[0.05,0.1,0.15,0.2,0.4]):
                mask = np.where((self.gas['Temperature'] > 10**6) & (self.gas['radius']* u.cm.to(u.Mpc) > minr ) & (self.gas['radius']* u.cm.to(u.Mpc) <= maxr ))
                print 'in mask'
                print self.gas['rel_speed'][mask]* u.cm.to(u.km)
                ax.hist(self.gas['rel_speed'][mask]* u.cm.to(u.km), color = 'red', alpha = 0.5, label = 'logT > 6')

                mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5) & (self.gas['radius']* u.cm.to(u.Mpc) > minr ) & (self.gas['radius']* u.cm.to(u.Mpc) <= maxr ))
                ax.hist(self.gas['rel_speed'][mask]* u.cm.to(u.km), color = 'orange', alpha = 0.5, label = '5 < logT < 6')

                mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4) & (self.gas['radius']* u.cm.to(u.Mpc) > minr ) & (self.gas['radius']* u.cm.to(u.Mpc) <= maxr ))
                ax.hist(self.gas['rel_speed'][mask]* u.cm.to(u.km), color = 'green', alpha = 0.5, label = '4 < logT < 5')

                mask = np.where(self.gas['Temperature'] <= 10**4 & (self.gas['radius']* u.cm.to(u.Mpc) > minr ) & (self.gas['radius']* u.cm.to(u.Mpc) <= maxr ))
                ax.hist(self.gas['rel_speed'][mask]* u.cm.to(u.km), color = 'blue', alpha = 0.8, label = 'logT < 4')

                ax.text(0.01,0.85,' %s Mpc < R < %s Mpc'%(minr,maxr), transform=ax.transAxes)

          #  hist, bin_edges = np.histogram(self.gas['rel_velocity'][mask], density=False)

            ax1.legend(framealpha=0.5, loc = 4)

            for ax in [ax1,ax2,ax3,ax4,ax5]:
                ax.minorticks_on()
                ax.set_axisbelow(True)
                ax.grid(linestyle='--')
               # ax.set_xlim(-1.1,1.1)
            for ax in [ax1,ax2,ax3,ax4]:
                ax.set_ylabel('number')
                ax.set_xticklabels([])
                ax.set_xlabel('')

            ax5.set_ylabel('number'); ax5.set_xlabel('relative velocity (km/s)')

            plt.tight_layout()
            fig.subplots_adjust( hspace = 0.07)

            plt.savefig('gasvelocities_hist_radbins_relvel_gn%s_sgn%s.png'%(gn,sgn))
            plt.close()

    def plotOabundance(self,gn, sgn, incSF):
        """Plot a map of abundance in the gas particles for Oxygen """

        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize = (7,7))

        'from make_maps_opts_locs: solar abundance is 4.8977835E-4'
        oxysol = 4.8977835E-4
        print 'solar oxygen abundance:'
        print oxysol
        print np.log10(oxysol)


        'grab the oxygen abundances into an array first, so can pick out limits for the plots'
        yy = self.gas['SmoothedElementAbundance/Oxygen']

        zmax = max(yy)
        print np.percentile(yy,[10.,25.,50.,75.,90.])
        zmin = np.percentile(yy,10.)
        print 'zmin,zmax'
        print zmin,zmax
        print 'log zmin,zmax'
        print np.log10(zmin),np.log10(zmax)

        print ""

        cm = plt.cm.get_cmap('RdYlBu')

        zmin = np.log10(zmin)
        zmax = np.log10(zmax)

        if incSF:
            mask = np.where(self.gas['Temperature'] > 10**6)
            print "Including SFing gas. Number of particles is:"
        else:
            mask = np.where((self.gas['Temperature'] > 10**6) & (self.gas['StarFormationRate'] == 0))
            print "Not including SFing gas. Number of particles is:"
      #  yy = self.gas['SmoothedElementAbundance/Oxygen'][mask]
        yy = np.log10(self.gas['SmoothedElementAbundance/Oxygen'][mask])
        print len(yy)
        sc = ax1.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),\
                         self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc),alpha=0.3,\
                         c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm,label = 'logT > 6')
        t = ax1.text(0.1,0.9,'logT > 6',transform=ax1.transAxes,color='black',fontsize=9)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))

        if incSF:
            mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5))
        else:
            mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5) & (self.gas['StarFormationRate'] == 0))
      #  yy = self.gas['SmoothedElementAbundance/Oxygen'][mask]
        yy = np.log10(self.gas['SmoothedElementAbundance/Oxygen'][mask])
        print len(yy)
        sc = ax2.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),\
                         self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc),alpha=0.3,\
                         c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm,label = '5 < logT < 6')
        t = ax2.text(0.1,0.9,'5 < logT < 6',transform=ax2.transAxes,color='black',fontsize=9)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))
        div = axgrid.make_axes_locatable(ax2)
        cax = div.append_axes("right",size="5%",pad=0.1)
        cbar = plt.colorbar(sc,label='', cax=cax,orientation='vertical')
        cbar.ax.set_ylabel(r'Oxygen Mass Fraction')#, fontsize=fontsize)

        if incSF:
            mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4))
        else:
            mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4) & (self.gas['StarFormationRate'] == 0))
     #   yy = self.gas['SmoothedElementAbundance/Oxygen'][mask]
        yy = np.log10(self.gas['SmoothedElementAbundance/Oxygen'][mask])
        print len(yy)
        sc = ax3.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),\
                         self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc),alpha=0.3,\
                         c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm,label = '4 < logT < 5')
        t = ax3.text(0.1,0.9,'4 < logT < 5',transform=ax3.transAxes,color='black',fontsize=9)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))

        if incSF:
            mask = np.where(self.gas['Temperature'] <= 10**4)
        else:
            mask = np.where((self.gas['Temperature'] <= 10**4) & (self.gas['StarFormationRate'] == 0))
     #   yy = self.gas['SmoothedElementAbundance/Oxygen'][mask]
        yy = np.log10(self.gas['SmoothedElementAbundance/Oxygen'][mask])
        print len(yy)
        sc = ax4.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),\
                         self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc),alpha=0.3,\
                         c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm,label = 'logT < 4')
        t = ax4.text(0.1,0.9,'logT < 4',transform=ax4.transAxes,color='black',fontsize=9)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))
        div = axgrid.make_axes_locatable(ax4)
        cax = div.append_axes("right",size="5%",pad=0.1)
        cbar = plt.colorbar(sc,label='', cax=cax,orientation='vertical')
        cbar.ax.set_ylabel(r'Oxygen Mass Fraction')#, fontsize=fontsize)

        for ax in [ax1,ax2,ax3,ax4]:
         #   ax.plot(centre[0]* u.Mpc.to(u.cm),centre[1]* u.Mpc.to(u.cm),'o',color='black')
           # ax.plot(centre[0],centre[1],'o',color='black')
            ax.set_ylabel('Y (cMpc)'); ax.set_xlabel('X (cMpc)')
            ax.minorticks_on()

        plt.tight_layout()
       # plt.show()
        SFind = '_noSF'
        if incSF:
            SFind=''
        plt.savefig('gasvelocities_oxygenabunds_gn%s_sgn%s%s.png'%(gn,sgn,SFind))
        plt.close()

    def plotlocsdensity(self,gn, sgn, incSF, useNH):
        """Plot a map of abundance in the gas particles for Oxygen """

        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize = (7,7))

        'grab the densities into an array first, so can pick out limits for the plots'
        if useNH:
            yy = self.gas['nH']
        else:
            yy = self.gas['Density']

        zmax = max(yy)
        print np.percentile(yy,[10.,25.,50.,75.,90.])
        zmin = np.percentile(yy,10.)
        print 'zmin,zmax'
        print zmin,zmax
        print 'log zmin,zmax'
        print np.log10(zmin),np.log10(zmax)
        print ""

        cm = plt.cm.get_cmap('RdYlBu')

        zmin = np.log10(zmin)
        zmax = np.log10(zmax)

        if incSF:
            mask = np.where(self.gas['Temperature'] > 10**6)
            print "Including SFing gas. Number of particles is:"
        else:
            mask = np.where((self.gas['Temperature'] > 10**6) & (self.gas['StarFormationRate'] == 0))
            print "Not including SFing gas. Number of particles is:"
      #  yy = self.gas['SmoothedElementAbundance/Oxygen'][mask]
        if useNH:
            yy = np.log10(self.gas['nH'][mask])
        else:
            yy = np.log10(self.gas['Density'][mask])
        print len(yy)
        sc = ax1.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),\
                         self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc),alpha=0.3,\
                         c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm,label = 'logT > 6')
        t = ax1.text(0.1,0.9,'logT > 6',transform=ax1.transAxes,color='black',fontsize=9)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))

        if incSF:
            mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5))
        else:
            mask = np.where((self.gas['Temperature'] <= 10**6) & (self.gas['Temperature'] > 10**5) & (self.gas['StarFormationRate'] == 0))
      #  yy = self.gas['SmoothedElementAbundance/Oxygen'][mask]
        if useNH:
            yy = np.log10(self.gas['nH'][mask])
        else:
            yy = np.log10(self.gas['Density'][mask])
        print len(yy)
        sc = ax2.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),\
                         self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc),alpha=0.3,\
                         c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm,label = '5 < logT < 6')
        t = ax2.text(0.1,0.9,'5 < logT < 6',transform=ax2.transAxes,color='black',fontsize=9)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))
        div = axgrid.make_axes_locatable(ax2)
        cax = div.append_axes("right",size="5%",pad=0.1)
        cbar = plt.colorbar(sc,label='', cax=cax,orientation='vertical')
        if useNH:
            cbar.ax.set_ylabel(r'log10 nH')#, fontsize=fontsize)
        else:
            cbar.ax.set_ylabel(r'log10 Density')#, fontsize=fontsize)

        if incSF:
            mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4))
        else:
            mask = np.where((self.gas['Temperature'] <= 10**5) & (self.gas['Temperature'] > 10**4) & (self.gas['StarFormationRate'] == 0))
     #   yy = self.gas['SmoothedElementAbundance/Oxygen'][mask]
        if useNH:
            yy = np.log10(self.gas['nH'][mask])
        else:
            yy = np.log10(self.gas['Density'][mask])
        print len(yy)
        sc = ax3.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),\
                         self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc),alpha=0.3,\
                         c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm,label = '4 < logT < 5')
        t = ax3.text(0.1,0.9,'4 < logT < 5',transform=ax3.transAxes,color='black',fontsize=9)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))

        if incSF:
            mask = np.where(self.gas['Temperature'] <= 10**4)
        else:
            mask = np.where((self.gas['Temperature'] <= 10**4) & (self.gas['StarFormationRate'] == 0))
     #   yy = self.gas['SmoothedElementAbundance/Oxygen'][mask]
        if useNH:
            yy = np.log10(self.gas['nH'][mask])
        else:
            yy = np.log10(self.gas['Density'][mask])
        print len(yy)
        sc = ax4.scatter(self.gas['Coordinates'][mask][:,0]* u.cm.to(u.Mpc),\
                         self.gas['Coordinates'][mask][:,1]* u.cm.to(u.Mpc),alpha=0.3,\
                         c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm,label = 'logT < 4')
        t = ax4.text(0.1,0.9,'logT < 4',transform=ax4.transAxes,color='black',fontsize=9)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='black'))
        div = axgrid.make_axes_locatable(ax4)
        cax = div.append_axes("right",size="5%",pad=0.1)
        cbar = plt.colorbar(sc,label='', cax=cax,orientation='vertical')
        if useNH:
            cbar.ax.set_ylabel(r'log10 nH')#, fontsize=fontsize)
        else:
            cbar.ax.set_ylabel(r'log10 Density')#, fontsize=fontsize)

        for ax in [ax1,ax2,ax3,ax4]:
         #   ax.plot(centre[0]* u.Mpc.to(u.cm),centre[1]* u.Mpc.to(u.cm),'o',color='black')
           # ax.plot(centre[0],centre[1],'o',color='black')
            ax.set_ylabel('Y (cMpc)'); ax.set_xlabel('X (cMpc)')
            ax.minorticks_on()

        plt.tight_layout()
        # plt.show()

        lognHind = ''
        if useNH:
            lognHind = '_nH'
        SFind = '_noSF'
        if incSF:
            SFind=''
        plt.savefig('gasvelocities_locsdensity_gn%s_sgn%s%s%s.png'%(gn,sgn,SFind,lognHind))
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
    """
    # Centre is the COP for GN=1 SGN=0 taken from the database.  ### THIS IS FOR THE 12Mpc REFERENCE BOX ###
    fileloc = '/Users/lokhorst/Data/EAGLE/RefL0012N0188/snapshot_028_z000p000/'
    centre = np.array([12.08808994,4.47437191,1.41333473])  # cMpc
    velocity = np.array([25.2508, 25.5319, 2.74462]) # km / s
    x = GasVelocities_ReadEagle(1, 0, centre,velocity, load_region_length=0.2, fileloc = fileloc)
    """
    # 25 Mpc REFERENCE BOX
    """galID, gn, subgn, x, y, z, vx, vy, vz
    296111, 12, 0, 12.7689, 7.78359, 5.0325, 17.1401, 72.837, 53.3104, 7.49494e+10"""
    ###fileloc = '/Users/lokhorst/Data/EAGLE/RefL0025N0376/snapshot_028_z000p000/'
    fileloc = '/Volumes/Cerulean/EAGLE/EagleParticleData/RefL0025N0376/snapshot_028_z000p000/'
    centre = np.array([12.7689, 7.78359, 5.0325])  # cMpc
    velocity = np.array([17.1401, 72.837, 53.3104]) # km / s
    x = GasVelocities_ReadEagle(12, 0, centre,velocity, load_region_length=0.2, fileloc = fileloc, incSF=False,useNH=True)

    # Use the COP for the galaxies that we consider in the paper (print off COP coords in extract_FOV_and_cutout_galaxies)
    """galID, gn, subgn, x, y, z, vx, vy, vz
    13738373, 134,  1, 50.0543327332, 14.3326339722, 11.5063114166, -185.157424927, -219.539093018,  25.2807445526
    13552548, 7631, 0, 50.2796440125, 13.9161596298, 11.1202640533, -98.7809524536,   64.9828720093, 472.819396973
    11911337, 4067, 0, 50.4869918823, 13.5846481323, 12.3390817642, -81.2321090698,   43.9984779358, 362.303588867 ## This one first
    17647764, 718,  0, 50.6496925354, 13.6481962204, 13.0192470551, -24.9452381134, -195.142959595,  131.791793823
    9958488,  1801, 0, 51.2028465271, 12.965212822,  13.6717453003, -43.2850875854, -199.545211792,  100.731079102
    13566041, 7701, 0, 51.2168502808, 14.3266563416, 13.995010376,  -16.8092041016, -294.485870361,  66.829208374
    """
    """
    ### THIS IS FOR THE 100Mpc REFERENCE BOX ###
  #  11911337, 4067, 0, 50.4869918823, 13.5846481323, 12.3390817642, -81.2321090698,   43.9984779358, 362.303588867 ## This one first
    fileloc = '/Users/lokhorst/Data/EAGLE/RefL0100N1504/snapshot_028_z000p000/'
    centre = np.array([50.4869918823, 13.5846481323, 12.3390817642])  # cMpc
    velocity = np.array([-81.2321090698, 43.9984779358, 362.303588867]) # km / s
    x = GasVelocities_ReadEagle(4067, 0, centre, velocity, load_region_length=0.2, fileloc = fileloc)


  #  13738373, 134,  1, 50.0543327332, 14.3326339722, 11.5063114166, -185.157424927, -219.539093018,  25.2807445526
    centre = np.array([50.0543327332, 14.3326339722, 11.5063114166])  # cMpc
    velocity = np.array([-185.157424927, -219.539093018,  25.2807445526]) # km / s
    x = GasVelocities_ReadEagle(134,  1, centre, velocity, load_region_length=0.2, fileloc = fileloc)

  #  13552548, 7631, 0, 50.2796440125, 13.9161596298, 11.1202640533, -98.7809524536,   64.9828720093, 472.819396973
    centre = np.array([50.2796440125, 13.9161596298, 11.1202640533])  # cMpc
    velocity = np.array([-98.7809524536,   64.9828720093, 472.819396973]) # km / s
    x = GasVelocities_ReadEagle(7631, 0, centre, velocity, load_region_length=0.2, fileloc = fileloc)

  #  17647764, 718,  0, 50.6496925354, 13.6481962204, 13.0192470551, -24.9452381134, -195.142959595,  131.791793823
    centre = np.array([50.6496925354, 13.6481962204, 13.0192470551])  # cMpc
    velocity = np.array([-24.9452381134, -195.142959595,  131.791793823]) # km / s
    x = GasVelocities_ReadEagle(718,  0, centre, velocity, load_region_length=0.2, fileloc = fileloc)

  #  9958488,  1801, 0, 51.2028465271, 12.965212822,  13.6717453003, -43.2850875854, -199.545211792,  100.731079102
    centre = np.array([51.2028465271, 12.965212822,  13.6717453003])  # cMpc
    velocity = np.array([-43.2850875854, -199.545211792,  100.731079102]) # km / s
    x = GasVelocities_ReadEagle(1801, 0, centre, velocity, load_region_length=0.2, fileloc = fileloc)

  #  13566041, 7701, 0, 51.2168502808, 14.3266563416, 13.995010376,  -16.8092041016, -294.485870361,  66.829208374
    centre = np.array([51.2168502808, 14.3266563416, 13.995010376])  # cMpc
    velocity = np.array([-16.8092041016, -294.485870361,  66.829208374]) # km / s
    x = GasVelocities_ReadEagle(7701, 0, centre, velocity, load_region_length=0.2, fileloc = fileloc)
    """
