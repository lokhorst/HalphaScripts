
"""
Usage: python integratUVB.py

Get UVB data from: http://www.ucolick.org/~pmadau/CUBA/DOWNLOADS.html
"""


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
    nu_naught:  13.6 eV is ionization energy of hydrogen --> 91 nm = 911 AA
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
        axis.plot(np.log10(uvbfrequency.value), np.log10(uvbspectrum.value))

    print "Min and max of what we are integrating:"
    print min(integratethis)
    print max(integratethis)
    print ""
    total = np.trapz( integratethis, uvbfrequency, axis=-1)

    if np.isinf(total).any():
        print("WARNING: inf value from integration.")

    return total/10**15

def doMeanIntegral(uvbspectrum,uvbfrequency,axis=None):
    """Integrate over the UVB spectrum
    Follow the functional form from Gould & Weinberg 1996, for optically thick gas
    illuminated by the UVB.
    integral of J(nu) d(nu) from nu_naught to infinity
    nu_naught:  13.6 eV is ionization energy of hydrogen --> 91 nm = 911 AA
    """
    # reverse them so we get a positive value
    uvbspectrum = uvbspectrum[::-1]
    uvbfrequency =uvbfrequency[::-1]

    # running into floating point error
    uvbspectrum = uvbspectrum*10**15

    integratethis = uvbspectrum

#    print uvbspectrum
#    print uvbfrequency
#    print integratethis

    if axis is not None:
        print uvbspectrum
        print uvbfrequency
        print integratethis
        axis.plot(np.log10(uvbfrequency.value), np.log10(uvbspectrum.value))

    print ""
    print "Min and max of what we are integrating:"
    print min(integratethis)
    print max(integratethis)
    print ""
    total = np.trapz( integratethis, uvbfrequency, axis=-1)

    if np.isinf(total).any():
        print("WARNING: inf value from integration.")

    return total/10**15

def doMeanIntegral2(uvbspectrum,uvbfrequency,axis=None):
    """Integrate over the UVB spectrum
    Follow the functional form from Gould & Weinberg 1996, for optically thick gas
    illuminated by the UVB.
    integral of J(nu) nu d(nu) from nu_naught to infinity
    nu_naught:  13.6 eV is ionization energy of hydrogen --> 91 nm = 911 AA
    """
    # reverse them so we get a positive value
    uvbspectrum = uvbspectrum[::-1]
    uvbfrequency =uvbfrequency[::-1]

    # running into floating point error
    uvbspectrum = uvbspectrum*10**15

    integratethis = uvbspectrum * uvbfrequency

#    print uvbspectrum
#    print uvbfrequency
#    print integratethis

    if axis is not None:
        print uvbspectrum
        print uvbfrequency
        print integratethis
        axis.plot(np.log10(uvbfrequency.value), np.log10(uvbspectrum.value))

    print "Min and max of what we are integrating:"
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

    uvbwavelength = uvbdata[0][1:]*u.Angstrom # wavelength in Angstroms
    uvbspectrum = uvbdata[1][1:]*u.erg /  u.s / u.cm**2 / u.Hz / u.sr  # J (in units of ergs/s/cm^2/Hz/sr) # redshift z = 0
    uvbspectrum2 = uvbdata[2][1:] # J (in units of ergs/s/cm^2/Hz/sr)
    uvbspectrum5 = uvbdata[5][1:] # J (in units of ergs/s/cm^2/Hz/sr)

    _, iuniq = np.unique(uvbwavelength, return_index=True)
    print iuniq
    print len(uvbwavelength)
    print len(iuniq)

    uvbwavelength = uvbwavelength[iuniq]
    uvbspectrum = uvbspectrum[iuniq]
    uvbspectrum2 = uvbspectrum2[iuniq]
    uvbspectrum5 = uvbspectrum5[iuniq]

#    #convert to photons from ergs
#    photon_wavelength = 656.3e-9*u.m
#    SB_cgs = 10**-22 *u.erg / u.arcsec**2 / u.cm**2 / u.s
#    SB_ph = SB_cgs * photon_wavelength/(const.h.to('erg s') * const.c.to('m/s')) * (206265.*u.arcsec)**2/u.sr

    uvbenergy = const.h.to('eV s')*const.c.to('m/s')/(uvbwavelength.to('m'))
    print uvbenergy.unit
    uvbfrequency = const.c.to('m/s')/(uvbwavelength.to('m'))
    print uvbfrequency.unit

    debug=False
    if debug:
        print len(uvbwavelength),len(uvbspectrum),len(uvbenergy),len(uvbfrequency)
        print const.h.to('eV s')
        print const.c.to('m/s')
        print const.h.to('eV s') * const.c.to('m/s')/13.6
        print 13.6*u.eV / const.h.to('eV s')

    plotchecks=True
    if plotchecks:
        fig,ax=plt.subplots(1,1,figsize=(10,10))
        ax.plot(uvbwavelength, np.log10(uvbspectrum.value))
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
        ax.plot(np.log10(uvbenergy.value/1000), (uvbspectrum.value*uvbfrequency.value*10**-7*(100)**2*10**9))
        ax.set_yscale('log')
        #ax.set_xlim(5,5000)
        #ax.set_ylim(-26,-20)
        ax.set_ylabel(r'J$\nu$ (Watt/m^2/sr)')
        ax.set_xlabel(r'h$\nu$ (keV)')
        ax.grid()
        plt.savefig('integrateUVB_JFreqvsEnergy.pdf')

    if plotchecks:
        fig,ax=plt.subplots(1,1,figsize=(10,10))
        ax.plot(np.log10(uvbfrequency.value), (uvbspectrum.value*uvbfrequency.value*10**-7*(100)**2*10**9))
        ax.set_yscale('log')
        #ax.set_xlim(5,5000)
        #ax.set_ylim(-26,-20)
        ax.set_ylabel(r'J$\nu$ (Watt/m^2/sr)')
        ax.set_xlabel(r'$\nu$ (Hz)')
        ax.grid()
        plt.savefig('integrateUVB_JFreqvsFreq.pdf')

    #uvbfrequency=uvbfrequency.value

    # do the integration
    nu_naught = 13.6*u.eV / const.h.to('eV s')
    total = doIntegral(uvbspectrum[uvbfrequency.value>=nu_naught.value],uvbfrequency[uvbfrequency.value>=nu_naught.value])
    #total = doIntegral(uvbspectrum,uvbfrequency)
    #fig,ax=plt.subplots(1,1,figsize=(10,10))
    #total = doIntegral(uvbspectrum[uvbfrequency.value>=3.2e+15],uvbfrequency[uvbfrequency.value>=3.2e+15],axis=ax)
    #plt.show()
    print "Integral of J(nu)/(nu) is: "
    print total
    print "Integral of J(nu)/(h nu) is: "
    print total/const.h.to('erg s')

    weighttotal = doMeanIntegral(uvbspectrum[uvbfrequency.value>=nu_naught.value],uvbfrequency[uvbfrequency.value>=nu_naught.value])
    print "Integral of J(nu) is: "
    print weighttotal
    print "[Integral of J(nu)]/[Integral of J(nu)/(nu)] is: "
    print weighttotal/total
    print "[Integral of J(nu)]/[Integral of J(nu)/(nu)] in units of nu_naught is: "
    print weighttotal/total/nu_naught

    print ""
    print weighttotal/(const.h.to('erg s') * const.c.to('m/s'))*121.6*10**-9*u.m

    weighttotal2 = doMeanIntegral2(uvbspectrum[uvbfrequency.value>=nu_naught.value],uvbfrequency[uvbfrequency.value>=nu_naught.value])
    print "Integral of J(nu)*nu is: "
    print weighttotal2
    print "[Integral of J(nu)*nu]/[Integral of J(nu)] is:"
    print weighttotal2/weighttotal
    print "[Integral of J(nu)*nu]/[Integral of J(nu)] in units of nu_naught is:"
    print weighttotal2/weighttotal/nu_naught
