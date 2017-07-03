
import numpy as np
import matplotlib.pylab as plt
import warnings

class PSF(object):

    """This class generates representations of various point-spread
    functions."""

    def __init__(self, flux, fwhm):
        self.flux = flux
        self.fwhm = fwhm

    def moffat_profile(self, beta=4.765, flux_fraction=1.0):
        """Profile with default beta suggested by Trujillo et al. 2001. This
            value of beta corresponds to Kolmogorov turbulence."""
        def moffat(theta):
            totalflux = flux_fraction*self.flux
            fwhm = self.fwhm
            ibeta = 1.0/beta
            alpha = fwhm/(2.0*np.sqrt(2.0**ibeta - 1.0))
            coefficient = totalflux*(beta-1.0)/(np.pi*(alpha**2))
            firstMoffat = coefficient * (1 + (theta/alpha)**2)**(-beta)
            return firstMoffat
        return moffat

    def two_moffat_profile(self, flux_fraction=1.0):
        """Simple analytical expression given by Racine (1996) to model the
        CCD-based PSF. Racine's analysis did not go far enough out to contain
        the full aureole component given by King (1971), which was determined
        photographically."""
        def two_moffat(theta):
            totalflux = flux_fraction*self.flux
            fwhm = self.fwhm

            # The first Moffat profile
            flux = 0.8*totalflux
            beta = 7.0
            ibeta = 1.0/beta
            alpha = fwhm/(2.0*np.sqrt(2.0**ibeta - 1.0))
            coefficient = flux*(beta-1.0)/(np.pi*(alpha**2))
            firstMoffat = coefficient * (1 + (theta/alpha)**2)**(-beta)

            # The second Moffat profile
            flux = 0.2*totalflux
            beta = 2.0
            ibeta = 1.0/beta
            alpha = fwhm/(2.0*np.sqrt(2.0**ibeta - 1.0))
            coefficient = flux*(beta-1.0)/(np.pi*(alpha**2))
            secondMoffat = coefficient * (1 + (theta/alpha)**2)**(-beta)

            # The total profile
            return firstMoffat + secondMoffat
        return two_moffat

    def aureole_profile(self, d0=50, flux_fraction=1.0):
        """Analytical function suggested by Racine (1996) for modelling the full
        PSF aureole. The default value of d0 matches the photographic data of
        King (1971). Note that d0 is specified in units of the FWHM. (It is 50
        and not 100 as in the Racine paper because it is in units of FWHM and
        not HW). Note that constructing a PSF with an Aureole will change the
        effective PSF (likely by a very small amount, assuming only a small
        fraction of the light is in the aureole) but this will not be correctly
        captured by the fwhm property of the PSF."""
        def aureole(theta):
            totalflux = flux_fraction*self.flux
            fwhm = self.fwhm
            d = d0 * fwhm
            coeff = totalflux/(8*np.pi*d**2)
            cosfac = np.cos(np.arctan(theta/(2*d)))
            return coeff*(cosfac**3.0)
        return aureole
