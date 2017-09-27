# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 16:46:28 2016

@author: wijers

Helper file for make_maps: containts lists, directories, etc. used in emission calculations
"""

import eagle_constants_and_units as c
import string as string
import ion_header as ion

# example call:
# mmap.emission_calc([6.25,6.25,6.25],12.5,12.5,2.5,800,800,'carbon',49,log=True,Wfig='test_CV_Sb_L0012N0188_snap28_zslice_0p2_box_center_GadSm_SmAb', Wtitle='Surface brightness', Wmin=-10, Wmax=None, Wlabel=r'$S_B$ [photons $s^{-1}$ $cm^{-2}$ $sr^{-1}$]', Wcol='gray',periodic = True, kernel = 'gadget')

#######################
#      choices        #
#######################
  
file_type = 'snap'
simdir='/disks/eagle/%s/' # use RECALIBRATED in L0025N0752, REFERENCE default

npzdir = '/net/luttero/data1/line_em_abs/temp/'
imgdir = ''

############################
#    setup and functions   #
############################


dir_emtab = '/disks/strw17/serena/Tables/%s/Tables/'
c_interpfile = '/home/wijers/plot_sims/make_maps_emission_coupled/HsmlAndProject/emission.so'
hsml_dir = '/home/wijers/plot_sims/make_maps_emission_coupled/HsmlAndProject/'
# redhift table monotonically increasing, starting at 0
dir_iontab = '/disks/strw17/serena/IonizationTables/HM01G+C/%s'

kernel_list = ['C2','gadget']
# desngb = 58 read out from sample hdf5 file (RunTimePars)
desngb = 58

# must be monotonic
zopts = \
['0.0000','0.2709','0.3988','0.4675','0.7778','0.9567','1.0006','1.2590',\
'1.4870','1.7370','2.0130','2.3160','2.4790','2.8290','3.0170','3.2140',\
'3.4210','3.6380','4.1050','4.6190','4.8950','5.4880']
elements = ['sulfur','silicon','oxygen','nitrogen','neon','magnesium','iron',\
'hydrogen','helium','carbon','calcium']
ions = ['al1', 'al2', 'al3',\
        'c1', 'c2', 'c3', 'c4', 'c5', 'c6',\
        'fe2', 'fe3', 'fe17',\
        'h1',\
        'he1', 'he2',\
        'mg1', 'mg2',
        'n2', 'n3', 'n4', 'n5', 'n6', 'n7',\
        'ne8', 'ne9', 'ne10',\
        'o1', 'o3', 'o4', 'o5', 'o6', 'o7', 'o8',\
        's5',\
        'si2', 'si3', 'si4', 'si13']
        
zpoints = [float(z) for z in zopts] 

# copied from Parameters/ChemicalElements in simulation hdf5 files. 
# Seems to be the same accross simulations (and it should be, if the same cooling tabels are used) 
# matches Wiersema, Schaye, Theuns et al. 2009 table 1 values
solar_abunds_ea = {'calcium':  6.435500108636916E-5,\
                'carbon':   0.002066543558612466,\
                'helium':   0.2805553376674652,\
                'hydrogen': 0.706497848033905,\
                'iron':     0.0011032151523977518,\
                'magnesium':5.907064187340438E-4,\
                'neon':     0.0014144604792818427,\
                'nitrogen': 8.356256294064224E-4,\
                'oxygen':   0.00549262436106801,\
                'silicon':  6.825873861089349E-4,\
                'sulfur':  4.0898521547205746E-4}
# Serena Bertone's element tables use different values (copied from calcium file): 
# 'element solar abundances in Cloudy. Number density relative to H.'
solar_abunds_sb = {'calcium':  2.290866E-6,\
                'carbon':   2.4547108E-4,\
                'helium':   0.077624775,\
                'hydrogen': 1.0,\
                'iron':     2.8183817E-5,\
                'magnesium':3.467368E-5,\
                'neon':     1.0E-4,\
                'nitrogen': 8.511386E-5,\
                'oxygen':   4.8977835E-4,\
                'silicon':  3.467368E-5,\
                'sulfur':  1.8197019E-5}

# use bertone abundances for comparison with table, converted from number density relative to hydrogen to mass fraction
solar_elts = solar_abunds_sb.keys()
totdens = sum(ion.atomw[string.capwords(elt)]*solar_abunds_sb[elt] for elt in solar_elts)
def abundconv(elt):
    return ion.atomw[string.capwords(elt)]*solar_abunds_sb[elt]/totdens
solar_abunds = {elt: abundconv(elt) for elt in solar_elts}

elements_ion = {'c1': 'carbon', 'c2': 'carbon', 'c3': 'carbon', 'c4': 'carbon', 'c5': 'carbon', 'c6': 'carbon',\
             'fe2': 'iron', 'fe3': 'iron', 'fe17': 'iron', \
             'h1': 'hydrogen', 'lyalpha': 'hydrogen', 'halpha': 'hydrogen',\
             'he1': 'helium', 'he2': 'helium',\
             'mg1': 'magnesium', 'mg2': 'magnesium',\
             'n2': 'nitrogen', 'n3': 'nitrogen', 'n4': 'nitrogen', 'n5': 'nitrogen', 'n6': 'nitrogen', 'n7': 'nitrogen',\
             'ne8': 'neon', 'ne9': 'neon', 'ne10': 'neon',\
             'o1': 'oxygen', 'o3': 'oxygen', 'o4': 'oxygen', 'o5': 'oxygen', 'o6': 'oxygen', 'o7': 'oxygen', 'o8': 'oxygen',\
             's5': 'sulfur',\
             'si2': 'silicon', 'si3': 'silicon', 'si4': 'silicon', 'si13': 'silicon'}

# for emission lines; not strictly ions (Serena Bertone and Freeke van de Voort line choices)            
line_nos_ion = {'c5': 49,'c6': 56, 'o6': 119,'n6': 65, 'n7': 69,'o7': 136, 'o8':149, 'halpha':9, 'lyalpha':1, 'fe17': 594}
line_eng_ion = {'c5': 307.88*c.ev_to_erg,'c6': 367.47*c.ev_to_erg, 'n6': 419.86*c.ev_to_erg, 'n7': 500.24*c.ev_to_erg, 'o6': 12.01*c.ev_to_erg, 'o7': 573.95*c.ev_to_erg, 'o8': 653.55*c.ev_to_erg, 'halpha':1.89*c.ev_to_erg, 'lyalpha':10.19*c.ev_to_erg, 'fe17': 726.97*c.ev_to_erg} # eV

snapzstr = {'024': '000p366', '025': '000p271','026': '000p183','027': '000p101', '028': '000p000'}
