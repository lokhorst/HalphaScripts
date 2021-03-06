# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:43:18 2017

@author: wijers

General cosmological utility functions; initially copied from make_maps to be
loaded without the entire read_Eagle machinery
"""

import numpy as np
import make_maps_opts_locs as ol # needed for some ion data
reload(ol)
import eagle_constants_and_units as c
import ctypes as ct
import h5py

def findemtables(element,zcalc):
    
    #### checks and setup
    
    if not element in ol.elements:
        print("There will be an error somewhere: %s is not included or misspelled. \n" % element)
    
    if zcalc < 0. and zcalc > 1e-4:
        zcalc = 0.0
        zname = ol.zopts[0]
        interp = False

    elif zcalc in ol.zpoints:
        # only need one table
        zname = ol.zopts[ol.zpoints.index(zcalc)]
        interp = False
    
    elif zcalc <= ol.zpoints[-1]:
        # linear interpolation between two tables
        zarray = np.asarray(ol.zpoints)
        zname1 = ol.zopts[len(zarray[zarray<zcalc])-1]
        zname2 = ol.zopts[-len(zarray[zarray>zcalc])]
        interp = True
    else:
        print("Chosen z value requires extrapolation. This has not been implemented. \n") 
        
    
    #### read in the tables; interpolate tables in z if needed and possible
    
    if not interp:
        tablefilename = ol.dir_emtab + zname + '/Tables/' + element + '.hdf5'
        tablefile = h5py.File(tablefilename, "r")
        #energies = np.array(tablefile.get('header/spectrum/logenergy_ryd'))
        #fluxes = np.array(tablefile.get('header/spectrum/logflux'))
        logTK =   np.array(tablefile.get('logt'),dtype=np.float32)  
        logrhocm3 =   np.array(tablefile.get('logd'),dtype=np.float32)
        lines =   np.array(tablefile.get('lines'),dtype=np.float32)  
        
        
        tablefile.close()
    
    if interp: #linear interpolation: 1./(a1-a0) * ( (a1-a)*f0 + (a-a0)*f1 )
        tablefilename1 = ol.dir_emtab + zname1 + '/Tables/' + element + '.hdf5'
        tablefile1 = h5py.File(tablefilename1, "r")
        #energies = np.array(tablefile.get('header/spectrum/logenergy_ryd'))
        #fluxes = np.array(tablefile.get('header/spectrum/logflux'))
        logTK1 =   np.array(tablefile1.get('logt'),dtype=np.float32)  
        logrhocm31 =   np.array(tablefile1.get('logd'),dtype=np.float32)
        lines1 =   np.array(tablefile1.get('lines'),dtype=np.float32) 
        
        tablefile1.close()
        
        tablefilename2 = ol.dir_emtab + zname2 + '/Tables/' + element + '.hdf5'
        tablefile2 = h5py.File(tablefilename2, "r")
        #energies = np.array(tablefile.get('header/spectrum/logenergy_ryd'))
        #fluxes = np.array(tablefile.get('header/spectrum/logflux'))
        logTK2 =   np.array(tablefile2.get('logt'),dtype=np.float32)  
        logrhocm32 =   np.array(tablefile2.get('logd'),dtype=np.float32)
        lines2 =   np.array(tablefile2.get('lines'),dtype=np.float32) 
        
        tablefile2.close()
        
        if (np.all(logTK1 == logTK2) and np.all(logrhocm31 == logrhocm32)):
            print("interpolating 2 emission tables")
            lines = 1./(float(zname2)-float(zname1)) * ( (float(zname2)-zcalc)*lines1 + (zcalc-float(zname1))*lines2 )
            logTK = logTK1
            logrhocm3 = logrhocm31
        else: 
            print("Temperature and density ranges of the two interpolation z tables don't match. \n")
            print("Using nearest z table in stead.")
            if abs(zcalc - float(zname1)) < abs(zcalc - float(zname2)):
                logTK = logTK1
                logrhocm3 = logrhocm31
                lines = lines1
            else:
                logTK = logTK2
                logrhocm3 = logrhocm32
                lines = lines2
    
    return lines, logTK, logrhocm3

           
# calculate emission using C function (interpolator)
def find_emdenssq(z,elt,lognH,logT,lineind):

    p_emtable, logTK, lognHcm3 = findemtables(elt,z)
    emtable = p_emtable[:,:,lineind]
    NumPart = len(lognH)
    inlogemission = np.zeros(NumPart,dtype=np.float32)
    
    if len(logT) != NumPart:
        print('logrho and logT should have the same length')
        return None

    # need to compile with some extra options to get this to work: make -f make_emission_only
    print("------------------- C interpolation function output --------------------------\n")
    cfile = ol.c_interpfile

    acfile = ct.CDLL(cfile)
    interpfunction = acfile.interpolate_emdenssq

    interpfunction.argtypes = [np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),\
                           ct.c_longlong , \
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(logTK)*len(lognHcm3),)), \
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(lognHcm3),)), \
                           ct.c_int,\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(logTK),)), \
                           ct.c_int,\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,))]

    # argument conversion

    res = interpfunction(lognH.astype(np.float32),\
               logT.astype(np.float32),\
               ct.c_longlong(NumPart),\
               np.ndarray.flatten(emtable.astype(np.float32)),\
               lognHcm3.astype(np.float32),\
               ct.c_int(len(lognHcm3)),\
               logTK.astype(np.float32),\
               ct.c_int(len(logTK)), \
               inlogemission \
              )
    
    print("-------------- C interpolation function output finished ----------------------\n")
    
    if res != 0:
        print('Something has gone wrong in the C function: output %s. \n',str(res))
        return None
        
    return inlogemission



def findiontables(ion,z):
    # README in dir_iontab:
    # files are hdf5, contain ionisation fraction of a species for rho, T, z
    
    
    #### checks and setup
    
    if not ion in ol.ions:
        print("There will be an error somewhere: %s is not included or misspelled. \n" % ion)
    
    tablefilename = ol.dir_iontab + ion + '.hdf5'      
    tablefile = h5py.File(tablefilename, "r")
    logTK =   np.array(tablefile.get('logt'),dtype=np.float32)  
    lognHcm3 =   np.array(tablefile.get('logd'),dtype=np.float32)
    ziopts = np.array(tablefile.get('redshift'),dtype=np.float32) # monotonically incresing, first entry is zero
    balance_d_t_z = np.array(tablefile.get('ionbal'),dtype=np.float32)
    tablefile.close()

    if z < 0. and z > 1e-4:
        z = 0.0
        zind = 0 
        interp = False
        
    elif z in ziopts:
        # only need one table
        zind = np.argwhere(z == ziopts)
        interp = False
    
    elif z <= ziopts[-1]:
        # linear interpolation between two tables
        zind1 = np.sum(ziopts<z)-1
        zind2 = -sum(ziopts>z)
        interp = True
    else:
        print("Chosen z value requires extrapolation. This has not been implemented. \n") 
        
    
    #### read in the tables; interpolate tables in z if needed and possible
    
    if not interp:
       balance = np.squeeze(balance_d_t_z[:,:,zind]) # for some reason, extra dimensions are tacked on 
    
    if interp: #linear interpolation: 1./(a1-a0) * ( (a1-a)*f0 + (a-a0)*f1 )
        balance1 = balance_d_t_z[:,:,zind1]
        balance2 = balance_d_t_z[:,:,zind2]
        
        print("interpolating 2 emission tables")
        balance = 1./( ziopts[zind2] - ziopts[zind1]) * ( (ziopts[zind2]-z)*balance1 + (z-ziopts[zind1])*balance2 )

    
    return balance, logTK, lognHcm3

def find_ionbal(z,ion,lognH,logT):
    
    # compared to the line emission files, the order of the nH, T indices in the balance tables is switched
    balance, logTK, lognHcm3 = findiontables(ion,z) #(np.array([[0.,0.],[0.,1.],[0.,2.]]), np.array([0.,1.,2.]), np.array([0.,1.]) ) 
    NumPart = len(lognH)
    inbalance = np.zeros(NumPart,dtype=np.float32)
    
    if len(logT) != NumPart:
        print('logrho and logT should have the same length')
        return None

    # need to compile with some extra options to get this to work: make -f make_emission_only
    print("------------------- C interpolation function output --------------------------\n")
    cfile = '/home/wijers/plot_sims/make_maps_emission_coupled/HsmlAndProject/emission.so'

    acfile = ct.CDLL(cfile)
    interpfunction = acfile.interpolate_emdenssq # just a linear interpolator; works for non-emission stuff too

    interpfunction.argtypes = [np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),\
                           ct.c_longlong , \
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(logTK)*len(lognHcm3),)), \
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(logTK),)), \
                           ct.c_int,\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(lognHcm3),)), \
                           ct.c_int,\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,))]

 
    res = interpfunction(logT.astype(np.float32),\
               lognH.astype(np.float32),\
               ct.c_longlong(NumPart),\
               np.ndarray.flatten(balance.astype(np.float32)),\
               logTK.astype(np.float32),\
               ct.c_int(len(logTK)), \
               lognHcm3.astype(np.float32),\
               ct.c_int(len(lognHcm3)),\
               inbalance \
              )

    print("-------------- C interpolation function output finished ----------------------\n")
    
    if res != 0:
        print('Something has gone wrong in the C function: output %s. \n',str(res))
        return None
        
    return inbalance


# cosmological basics

def comoving_distance_cm(z): # assumes Omega_k = 0
    if z < 1e-8:
        print('Using 0 comoving distance from z. \n')
        return 0.
        
    def integrand(zi):
        return (c.omega0*(1.+zi)**3 + c.omegalambda)**0.5
    zi_arr = np.arange(0,z+z/512.,z/512.)
    com = np.trapz(1./integrand(zi_arr),x=zi_arr)
    return com * c.c/(c.hubble*c.hubbleparam)

def ang_diam_distance_cm(z):
    return comoving_distance_cm(z)/(1.+z)

def lum_distance_cm(z):
    return comoving_distance_cm(z)*(1.+z)

def Hubble(z):
    return c.hubble*c.hubbleparam*(c.omega0*(1.+z)**3 + c.omegalambda)**0.5   #km/cm * cm/Mpc
      
def solidangle(alpha,beta): # alpha = 0.5 * pix_length_1/D_A, beta = 0.5 * pix_length_2/D_A
    #from www.mpia.de/~mathar/public/mathar20051002.pdf
    # citing  A. Khadjavi, J. Opt. Soc. Am. 58, 1417 (1968).
    # stored in home/papers
    # using the exact formula, with alpha = beta, 
    # the python exact formula gives zero for alpha = beta < 10^-3--10^-4
    # assuming the pixel sizes are not so odd that the exact formula is needed in one direction and gives zero for the other,
    # use the Taylor expansion to order 4 
    # testing the difference between the Taylor and python exact calculations shows that 
    # for square pixels, 10**-2.5 is a reasonable cut-off
    # for rectangular pixels, the cut-off seems to be needed in both values
    if alpha < 10**-2.5 or beta <10**-2.5:
        return 4*alpha*beta - 2*alpha*beta*(alpha**2+beta**2)
    else: 
        return 4*np.arccos(((1+alpha**2 +beta**2)/((1+alpha**2)*(1+beta**2)))**0.5)
