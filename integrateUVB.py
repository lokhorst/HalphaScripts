

import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as u
from astropy.io import fits

def doIntegral(uvbspectrum,uvbfrequency,axis=None):
    """Integrate over the UVB spectrum
    Follow the functional form from Gould & Weinberg 1996, for optically thick gas
    illuminated by the UVB.
    integral of J(nu) / (nu) d(nu) from nu_naught to infinity
    what is nu_naught?  13.6 eV is ionization energy of hydrogen --> 91 nm = 911 AA
    """
    # reverse them so we get a positive value
    uvbspectrum = uvbspectrum[::-1]
    uvbfrequency =uvbfrequency[::-1]

    # running into floating point error
    uvbspectrum = uvbspectrum*10**15

    integratethis = uvbspectrum / (uvbfrequency)

    if axis is not None:
        print uvbspectrum
        print uvbfrequency
        print integratethis
        axis.plot(uvbfrequency, integratethis)


    print ""
    print min(integratethis)
    print max(integratethis)
    print ""
    total = np.trapz( integratethis, uvbfrequency, axis=-1)

    if np.isinf(total).any():
        print("WARNING: inf value from integration.")

    return total/10**15

if __name__ == "__main__":
    # load the uvb data
    uvbdata = np.genfromtxt('HMUVB_CUBA.txt',dtype='f',comments='#',skip_header=1)
    uvbdata = uvbdata.transpose()
    #print uvbdata.shape
    uvbwavelength = uvbdata[0][1:] # wavelength in Angstroms
    #print uvbwavelength
    uvbspectrum = uvbdata[1][1:] # J (in units of ergs/s/cm^2/Hz/sr)
    uvbspectrum2 = uvbdata[2][1:] # J (in units of ergs/s/cm^2/Hz/sr)
    uvbspectrum5 = uvbdata[5][1:] # J (in units of ergs/s/cm^2/Hz/sr)
    print len(uvbwavelength),len(uvbspectrum)
    #print uvbspectrum
    # convert to photons from ergs
#    photon_wavelength = 656.3e-9*u.m
#    SB_cgs = 10**-22 *u.erg / u.arcsec**2 / u.cm**2 / u.s
#    SB_ph = SB_cgs * photon_wavelength/(const.h.to('erg s') * const.c.to('m/s')) * (206265.*u.arcsec)**2/u.sr

    print const.h.to('eV s')
    print const.c.to('m/s')
    #print const.h.to('eV s')*const.c.to('m/s')/(uvbwavelength*10**-10)
    uvbenergy = const.h.to('eV s')*const.c.to('m/s')/(uvbwavelength*10**-10)
    #print uvbenergy.value

    uvbfrequency = const.c.to('m/s')/(uvbwavelength*10**-10)
    #print uvbfrequency.value

    print const.h.to('eV s') * const.c.to('m/s')/13.6
    print 13.6*u.eV / const.h.to('eV s')

    plotchecks=True
    if plotchecks:
        fig,ax=plt.subplots(1,1,figsize=(10,10))
        ax.plot(uvbwavelength, np.log10(uvbspectrum))
        ax.plot(uvbwavelength, np.log10(uvbspectrum2))
        ax.plot(uvbwavelength, np.log10(uvbspectrum5))
        ax.set_xscale('log')
        #ax.set_xlim(5,5000)
        #ax.set_ylim(-26,-20)
        ax.set_ylabel('J (ergs/s/cm^2/Hz/sr)')
        ax.set_xlabel('wavelength (Angstroms)')
        ax.grid()
        plt.savefig('integrateUVB_JvsWavelength.pdf')

    if plotchecks:
        fig,ax=plt.subplots(1,1,figsize=(10,10))
        ax.plot(np.log10(uvbenergy.value/1000), (uvbspectrum*uvbfrequency.value*10**-7*(100)**2*10**9))
        ax.set_yscale('log')
        #ax.set_xlim(5,5000)
        #ax.set_ylim(-26,-20)
        ax.set_ylabel(r'J$\nu$ (Watt/m^2/sr)')
        ax.set_xlabel(r'h$\nu$ (keV)')
        ax.grid()
        plt.savefig('integrateUVB_JFreqvsEnergy.pdf')

    if plotchecks:
        fig,ax=plt.subplots(1,1,figsize=(10,10))
        ax.plot(np.log10(uvbfrequency.value), (uvbspectrum*uvbfrequency.value*10**-7*(100)**2*10**9))
        ax.set_yscale('log')
        #ax.set_xlim(5,5000)
        #ax.set_ylim(-26,-20)
        ax.set_ylabel(r'J$\nu$ (Watt/m^2/sr)')
        ax.set_xlabel(r'$\nu$ (Hz)')
        ax.grid()
        plt.savefig('integrateUVB_JFreqvsFreq.pdf')

    uvbfrequency=uvbfrequency.value

    # do the integration
    total = doIntegral(uvbspectrum[uvbfrequency>=3.28846551546e+15],uvbfrequency[uvbfrequency>=3.28846551546e+15])
    total = doIntegral(uvbspectrum,uvbfrequency)
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    total = doIntegral(uvbspectrum[uvbfrequency>=3.2e+15],uvbfrequency[uvbfrequency>=3.2e+15],axis=ax)
    #plt.show()
    print total
    print total/const.h.to('erg s')
