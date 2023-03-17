from math import *
import numpy as np
import sys
from scipy.optimize import minimize

def xg (Rns, Mns):
    """
    |

    Compute compactenss parameters of neutron star.

    :math:`x_g = 2 G M / (R c^2)`

    :param Rns: neutron star radius [km]
    :param Mns: neutron star mass [Solar mass]

    :returns: compactness parameter [uniteless]


    """

    G = 6.67430e-8   ## cgs
    Msol = 2e33      ## Solar mass in gramms
    R = Rns * 1e5    ## cm
    M = Mns * Msol
    c   = 2.998e+10  ## speed of light in cm/s

    res = 2.0 * G * M / R / c / c


    return res



def g14 (Rns, Mns):
    """
    |

    Calculate the **free fall acceleration** at the surface of neutron star.

    :math:`g_{14} = G M / (R^2 \sqrt{1 - x_g}) / 10^{14}`

    :param Rns: neutron star radius [km]
    :param Mns: neutron star mass [Solar mass]

    :returns: free fall acceleration at the surface of neutron star [g/1e14 cm/s/s]
    

    """

    G = 6.67430e-8   ## cgs
    Msol = 2e33      ## Solar mass in gramms
    R = Rns * 1e5    ## cm
    M = Mns * Msol
    c   = 2.998e+10  ## speed of light in cm/s

    xg = 2.0 * G * M / R / c / c
    g14_val = G*M / R / R / sqrt(1.0 - xg ) / 1e14 ## Gudmundsson et al. (1983), eq. (2-3)
    return g14_val

## Function to compute the luminosity of neutron star

def compute_L (theta, phi, Rns, Tmap):
    """
    |

    Calculate the **total thermal X-ray luminosity** of neutron star.

    :math:`L = \sigma_\mathrm{SB} R_\mathrm{NS}^2 \int_{4\pi} T_s^4 (\\theta, \\phi) \\sin \\theta d\\theta d\\phi`

    :param Tmap: Tmap is the surface thermal map [K]
    :param theta: list magnetic latitude [radians] where Tmap is provided 
    :param phi: list of magnetic longuitudets [radians] where Tmap is provided
    :param Rns: radius of neutron star [km]
           
    The thermal map can be one dimensional (in the case of axisymmetric temperature distribution) or two-dimensional.

    :returns: total thermal luminosity [erg/s]
    

    """

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    Bp = 1e14         ## G
    c   = 2.998e+10   ## speed of light in cm/s
    
    R = Rns * 1e5     ## cm
    
    dtheta = theta[1] - theta[0]
    dphi   = phi[1] - phi[0] 
    
    L = 0

    if len(Tmap.shape) == 1:

        for i in range (0, len(phi)):
            for j in range (0, len(theta)):
                L = L + sigma_SB * pow(R, 2.0) * sin(theta[j]) * pow(Tmap [j], 4) * dphi * dtheta 

    elif len(Tmap.shape) == 2:

        for i in range (0, len(phi)):
            for j in range (0, len(theta)):
                L = L + sigma_SB * pow(R, 2.0) * sin(theta[j]) * pow(Tmap [i,j], 4) * dphi * dtheta 

            
    return L

def compute_L_param (param, Teff, Rns, Mns):
    """
    |

    Calculate the **total thermal X-ray luminosity** produced by single/two or more parameter blackbody.

    :param param: list of areas and relative temperatures
    :param Teff: effective temperature
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]

    :returns: total thermal luminosity [erg/s]

    :example: compute_L_param ([1, 1], 1e6, 12, 1.4) computes luminosity of neutron star with surface temperature 1 MK, radius 12 km and mass 1.4 solar mass
    

    """

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    R = Rns * 1e5

    fc = sigma_SB * R*R * 4.0 * pi

    if len(param) == 2:
        res = param[0] * pow(param[1] * Teff, 4)
    elif len(param) == 4:
        res = param[0] * pow(param[2] * Teff, 4) + param[1] * pow(param[3] * Teff, 4)
    else:
        print ('Function works only with param size 2 and 4')
        sys.exit(0)



    res = fc * res

    return res


## Function to compute the effective temperature

def compute_Teff (theta, phi, Rns, Tmap):
    """
    |

    Calculate the effective temperature of neutron star.

    :param Tmap: Tmap is the surface thermal map [K]
    :param theta: list magnetic latitude [radians] where Tmap is provided 
    :param phi: list of magnetic longuitudets [radians] where Tmap is provided
    :param Rns: radius of neutron star [km]

    The thermal map can be one dimensional (in the case of axisymmetric temperature distribution) or two-dimensional.

    :returns: Effective temperature [K]
    

    """


    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.

    Lns = compute_L (theta, phi, Rns, Tmap)
    return (pow(Lns / (4.0 * pi * sigma_SB * pow(Rns*1e5, 2)) , 1.0 / 4.0))

## Lensing factor following the article Poutanen (2020) A&A 640, A24 (2020)

def D_factor (cospsi, Rns, Mns):
    """
    |

    Calculate the lensing factor following the article Poutanen (2020) A&A 640, A24 (2020)

    :param cospsi: cosine of propagation angle
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass] 

    :returns: Dimensionless lensing factor
    

    """



    G   = 6.67430e-8  ## cgs
    R   = Rns * 1e5   ## cm
    M   = Mns * 2e33  ## 1.4 Solar mass in g
    c   = 2.998e+10   ## speed of light in cm/s
    Rs  = 2.0 * G * M / c / c
    
    ## key parameters of the approximation
    u   = Rs / R
    y   = 1.0 - cospsi
    res = 1.0 + 3.0 * u*u*y*y / 112.0 - e / 100.0 * u*y* (2.0 * np.log(1.0 - y/2.0) + y*(1.0 - 3.0*y/4.0) / (1.0 - y/2.0))
    return res

## Function to compute the emission angle psi following the article Poutanen (2020) A&A 640, A24 (2020) 

def alpha (cospsi, Rns, Mns):
    G   = 6.67430e-8  ## cgs
    #RNS = 12e5        ## 12 km in cm
    R = Rns * 1e5
    M   = Mns * 2e33  ## 1.4 Solar mass in g
    c   = 2.998e+10   ## speed of light in cm/s
    Rs  = 2.0 * G * M / c / c
    
    ## key parameters of the approximation
    u   = Rs / R
    y   = 1.0 - cospsi
    x = (1.0 - u)*y*(1.0 + u*u*y*y/112.0 - e/100.0*u*y*(np.log(1-y/2.0) + y/2))
    res = 1.0 - x
    return np.arccos(res)



## Function to compute the emission angle psi following the article Poutanen (2020) A&A 640, A24 (2020) returns cos(alpha) 

def cos_alpha (cospsi, Rns, Mns):
    G   = 6.67430e-8  ## cgs
    #RNS = 12e5        ## 12 km in cm
    R = Rns * 1e5
    M   = Mns * 2e33  ## 1.4 Solar mass in g
    c   = 2.998e+10   ## speed of light in cm/s
    Rs  = 2.0 * G * M / c / c
    
    ## key parameters of the approximation
    u   = Rs / R
    y   = 1.0 - cospsi
    x = (1.0 - u)*y*(1.0 + u*u*y*y/112.0 - e/100.0*u*y*(np.log(1-y/2.0) + y/2))
    res = 1.0 - x
    return res


## Single blackbody spectra
## Result is computed in units: erg s^{-1} cm^{-2} keV^{-1}
def single_BB (Teff, Rns, Mns):

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    c   = 2.998e+10   ## speed of light in cm/s

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c
    g14 = G*Mns*Msol / R / R / sqrt(1.0 - xg ) / 1e14 ## Gudmundsson et al. (1983), eq. (2-3)


    Eph = np.logspace (-1.2, 1.62, 142) ## keV

    Teff_inf = Teff * sqrt(1 - xg)

    res = 15.0*sigma_SB / ( pow(pi, 4) * pow(kB, 4)) * np.power(Eph, 3) / (np.exp(Eph / (kB *Teff_inf)) - 1.0) / (1.0 - xg) 

    max_sp = np.max(res)

    for i in range (0, len(res)):
        if res[i] < 1: #max_sp / 1e8:
            res[i] = 0.0

    return [Eph, res]


def single_BB_photons (Teff, Rns, Mns, eph, nphot, integ = True):
    """
    |

    Calculate thermal spectra [number of photons per energy bin] emitted by neutron star with fixed surface temperature.

    :param Teff: effective temperature of neutron star [K]
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list of energies where spectra is to be calculated [keV]
    :param nphot: total number of received photons
    :param integ: if True it returns integer counts (0.3 means no counts in a bin), if False it returns fraction of counts as well.

    :returns: redshifted spectra [number of photons per energy bin]
    

    """

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    c   = 2.998e+10   ## speed of light in cm/s

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c
    g14 = G*Mns*Msol / R / R / sqrt(1.0 - xg ) / 1e14 ## Gudmundsson et al. (1983), eq. (2-3)

    Teff_inf = Teff * sqrt(1 - xg)

    spec = 15.0*sigma_SB / ( pow(pi, 4) * pow(kB, 4)) * np.power(eph, 2) / (np.exp(eph / (kB *Teff_inf)) - 1.0) / (1.0 - xg) 

    spec = spec / np.sum(spec) * nphot

    if integ:
        spec = np.asarray(spec, dtype=int)

    return spec

def examine_spectral_fit_1BB_photons (param, Teff, Rns, Mns, eph, nphot, L, integ = True):

    """
    |

    Function to examine the quality of single blackbody fit [number of photons per energy bin] in comparison the neutron star spectra.

    :param param: :sc: relative area of hot region, :pc: relative temperature of hot region
    :param Teff: effective temperature of neutron star [K]
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list of energies where spectra is to be calculated [keV]
    :param nphot: total number of received photons
    :param L: luminosity of neutron star
    :param integ: if True it returns integer counts (0.3 means no counts in a bin), if False it returns fraction of counts as well.

    :returns: redshifted spectra [number of photons per energy bin]


    """


    spec = single_BB_photons (param[1]*Teff, Rns, Mns, eph, nphot, integ = False)

    Lcomp = compute_L_param (param, Teff, Rns, Mns)

    #print ('Lcomp = ', Lcomp)

    spec =  spec * Lcomp / L

    if integ:
        spec = np.asarray(spec, dtype=int)

    return spec

def two_BB_photons (param, Teff, Rns, Mns, eph, nphot, integ=True):
    """
    |

    Calculate thermal spectra [number of photons per energy bin] composed of two blackbodies.

    :param param: array of four parameters [sc, sh, pc, ph]
    :sc: fraction of total surface area covered with first hot spot
    :sh: fraction of total surface area covered with second hot spot
    :pc: temperature of first hot spot as a fraction of Teff 
    :ph: temperature of second hot spot as a fraction of Teff 
    :param Teff: effective temperature of neutron star
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list of energies where spectra is to be calculated [keV]
    :param nphot: total number of received photons
    :param integ: integ = True correspond to integer result and integ = False corresponds to float result 

    :returns: :sp: redshifted spectra [number of photons per energy bin]
    

    """

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.

    R = Rns * 1e5

    sc, sh, pc, ph = param

    sc = abs(sc)
    sh = abs(sh)
    pc = abs(pc)
    ph = abs(ph)

    sp_Teff1 = single_BB_photons (pc*Teff, Rns, Mns, eph, nphot, integ=False)
    sp_Teff2 = single_BB_photons (ph*Teff, Rns, Mns, eph, nphot, integ=False)

    spec = sc * sp_Teff1 +  sh * sp_Teff2

    spec = spec / np.sum(spec) * nphot 

    if integ:
        spec = np.asarray(spec, dtype=int)

    return spec

def examine_spectral_fit_2BB_photons (param, Teff, Rns, Mns, eph, nphot, L, integ = True):

    """
    |

    Function to examine the quality of two blackbody fit [number of photons per energy bin] in comparison the neutron star spectra.

    :param param: :sc: relative area of hot region, :pc: relative temperature of hot region
    :param Teff: effective temperature of neutron star [K]
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list of energies where spectra is to be calculated [keV]
    :param nphot: total number of received photons
    :param L: luminosity of neutron star
    :param integ: if True it returns integer counts (0.3 means no counts in a bin), if False it returns fraction of counts as well.

    :returns: redshifted spectra [number of photons per energy bin]


    """


    spec = two_BB_photons (param, Teff, Rns, Mns, eph, nphot, integ = False)

    Lcomp = compute_L_param (param, Teff, Rns, Mns)

    spec = np.asarray(spec) * Lcomp / L

    if integ:
        spec = np.asarray(spec, dtype=int)

    return spec



def get_redshifted_spectra_pole_3D (theta, phi, Tmap, Rns, Mns):
    """
    |

    Calculate thermal spectra emitted by neutron star observed from its magnetic pole (top row of thermal map).

    :param Tmap: two dimensional array describing the surface thermal map [K]
    :param theta: list of magnetic latitude [radians] where Tmap is provided
    :param phi: list of magnetic longuitudets [radians] where Tmap is provided
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]

    :returns: :eph: logarithmic energy mesh [keV]
              :sp: redshifted spectra [erg s^-1 cm^-2 keV^-1]
              :visible_map: two-dimensional array which contains only the temperature distribution on visible hemisphere
    

    """


    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    c   = 2.998e+10   ## speed of light in cm/s

    map_of_visible = np.zeros ((len(phi), len(theta)))

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c

    Ts_inf = Tmap * sqrt(1 - xg)

    ## Here we prepare variables for integration over the visible hemisphere

    Eph = np.logspace (-1.2, 1.62, 142) ## keV
    #Eph = np.linspace (0.063, 41.68, 142)


    sp_red = np.zeros(142)

    dtheta = theta[1] - theta[0]
    dphi   = phi[1] - phi[0]

    en_red = 1.0 / (1.0 - xg)

    tst = 0

    for i in range (0, len(phi)):
        for j in range (0, len(theta)):
            al = alpha (cos(theta[j]), Rns, Mns)


            if (al < pi / 2.0) and (al > -pi / 2):
                Df = D_factor (cos(theta[j]), Rns, Mns)
                sp_red = sp_red +  Df * 15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(theta[j]) * np.cos(al)  * np.power(Eph, 3) / (np.exp(Eph / kB / Ts_inf[i,j]) - 1.0) * dtheta * dphi
                map_of_visible[i,j] = Ts_inf[i, j]

    max_sp = np.max(sp_red)

    for i in range (0, len(sp_red)):
        if sp_red[i] < 1.0: #max_sp / 1e8:
            sp_red[i] = 0.0

    return [Eph, sp_red, map_of_visible]

def spectra_pole (Tmap, Rns, Mns, eph):
    """
    |

    Calculate redshifted thermal spectra [erg s^-1 cm^-2 keV^-1] emitted by neutron star observed from its magnetic pole (top row of thermal map).

    :math:`H_E^\\infty = \\frac{15 \\sigma_{SB}}{\\pi^5 k_B^4} \\int_\\mathrm{viz} \\frac{E^3 \\cos \\alpha \\; \\mathcal{D} \\sin \\theta d\\theta d\\phi}{\\exp(E/(k_B T_s^\\infty)) - 1}`

    :param Tmap: member of Tmap class containing surface temperature map.
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list containing energy mesh :math:`E` [keV]

    :returns: :sp: redshifted spectra [erg s^-1 cm^-2 keV^-1]
              :visible_map: two-dimensional array which contains only the temperature distribution on visible hemisphere
    

    """

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    c   = 2.998e+10   ## speed of light in cm/s

    map_of_visible = np.zeros ((len(Tmap.phi), len(Tmap.theta)))

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c

    Ts_inf = Tmap.Ts * sqrt(1 - xg)

    sp_red = np.zeros(len(eph))

    dtheta = Tmap.theta[1] - Tmap.theta[0]
    dphi   = Tmap.phi[1] - Tmap.phi[0]

    en_red = 1.0 / (1.0 - xg)

    for i in range (0, len(Tmap.phi)):
        for j in range (0, len(Tmap.theta)):
            al = alpha (cos(Tmap.theta[j]), Rns, Mns)

            if al < pi / 2.0:
                Df = D_factor (cos(Tmap.theta[j]), Rns, Mns)
                sp_red = sp_red +  Df * 15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(Tmap.theta[j]) * np.cos(al)  * np.power(eph, 3) / (np.exp(eph / kB / Ts_inf[i,j]) - 1.0) * dtheta * dphi
                map_of_visible[i,j] = Ts_inf[i, j]

    return [sp_red, map_of_visible]






def get_redshifted_spectra_pole_photons (theta, phi, Tmap, Rns, Mns, eph, nphot):
    """
    |

    Calculate thermal spectra [integer number of photons per energy bin] emitted by neutron star observed from its magnetic pole (top row of thermal map).

    :param Tmap: two dimensional array describing the surface thermal map [K]
    :param theta: list of magnetic latitude [radians] where Tmap is provided
    :param phi: list of magnetic longuitudets [radians] where Tmap is provided
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list energies where spectra should be computed [keV]
    :param nphot: total number of photons to be generated

    :returns: :sp: number of photons per energy bin
              :visible_map: two-dimensional array which contains only the temperature distribution on visible hemisphere
    

    """

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    c   = 2.998e+10   ## speed of light in cm/s

    map_of_visible = np.zeros ((len(phi), len(theta)))

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c

    Ts_inf = Tmap * sqrt(1 - xg)

    ## Here we prepare variables for integration over the visible hemisphere

    sp_red = np.zeros(len(eph))

    dtheta = theta[1] - theta[0]
    dphi   = phi[1] - phi[0]

    for i in range (0, len(phi)):
        for j in range (0, len(theta)):
            al = alpha (cos(theta[j]), Rns, Mns)


            if al < pi / 2.0:
                Df = D_factor (cos(theta[j]), Rns, Mns)
                sp_red = sp_red +  Df * 15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(theta[j]) * np.cos(al)  * np.power(eph, 2) / (np.exp(eph / kB / Ts_inf[i,j]) - 1.0) * dtheta * dphi
                map_of_visible[i,j] = Ts_inf[i, j]

    coeff = nphot / np.sum(sp_red)

    sp_red_n = np.asarray(sp_red * coeff, dtype=int)


    return [sp_red_n, map_of_visible]


def get_redshifted_spectra_equator_3D (theta, phi, Tmap, Rns, Mns):
    
    map_of_visible = np.zeros ((len(phi), len(theta)))
    
    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    Bp = 1e14         ## G
    c   = 2.998e+10   ## speed of light in cm/s
    
    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c
    g14 = G*Mns*Msol / R / R / sqrt(1.0 - xg ) / 1e14 ## Gudmundsson et al. (1983), eq. (2-3)
    
    #print ('xg = ', xg, ' log_10 of  mean non-redshifted Ts = ', log10(np.mean(Tmap)))
    
    Ts_inf = Tmap * sqrt(1 - xg)
    
    ## Here we prepare variables for integration over the visible hemisphere
    
    Eph = np.logspace (-1.2, 1.62, 142) ## keV
    logEph = np.log10(Eph)
    
    sp_red = np.zeros(142)

    dtheta = theta[1] - theta[0]
    dphi   = phi[1] - phi[0] 
    
    en_red = 1.0 / (1.0 - xg)

    for i in range (0, len(phi)):
        for j in range (0, len(theta)):
            new_theta = sin(theta[j])*cos(phi[i])
            al = alpha (new_theta, Rns, Mns)
            
            
            
            #print (phi[i], theta[j], al)
            if al < pi / 2.0:
                Df = D_factor (new_theta, Rns, Mns)
                sp_red = sp_red +  Df * 15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(theta[j]) * np.cos(al) * np.power(Eph, 3) / (np.exp(Eph / kB / Ts_inf[i,j]) - 1.0) * dtheta * dphi            
                map_of_visible[i,j] = Ts_inf[i, j]
                
    max_sp = np.max(sp_red)
    
    for i in range (0, len(sp_red)):
        if sp_red[i] < 1:
            sp_red[i] = 0.0

    return [Eph, sp_red, map_of_visible]


def get_redshifted_spectra_equator_photons (theta, phi, Tmap, Rns, Mns, eph, nphot):
    """
    |

    Calculate thermal spectra [integer number of photons per energy bin] emitted by neutron star observed from its magnetic equator (top row of thermal map).

    :param Tmap: two dimensional array describing the surface thermal map [K]
    :param theta: list of magnetic latitude [radians] where Tmap is provided
    :param phi: list of magnetic longuitudets [radians] where Tmap is provided
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list energies where spectra should be computed [keV]
    :param nphot: total number of photons to be generated

    :returns: :sp: number of photons per energy bin
              :visible_map: two-dimensional array which contains only the temperature distribution on visible hemisphere
    

    """

    map_of_visible = np.zeros ((len(phi), len(theta)))

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    Bp = 1e14         ## G
    c   = 2.998e+10   ## speed of light in cm/s

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c
    g14 = G*Mns*Msol / R / R / sqrt(1.0 - xg ) / 1e14 ## Gudmundsson et al. (1983), eq. (2-3)

    #print ('xg = ', xg, ' log_10 of  mean non-redshifted Ts = ', log10(np.mean(Tmap)))

    Ts_inf = Tmap * sqrt(1 - xg)

    ## Here we prepare variables for integration over the visible hemisphere

    sp_red = np.zeros(len(eph))

    dtheta = theta[1] - theta[0]
    dphi   = phi[1] - phi[0]

    for i in range (0, len(phi)):
        for j in range (0, len(theta)):
            new_theta = sin(theta[j])*cos(phi[i])
            al = alpha (new_theta, Rns, Mns)



            #print (phi[i], theta[j], al)
            if al < pi / 2.0:
                Df = D_factor (new_theta, Rns, Mns)
                sp_red = sp_red +  Df * 15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(theta[j]) * np.cos(al) * np.power(eph, 2) / (np.exp(eph / kB / Ts_inf[i,j]) - 1.0) * dtheta * dphi
                map_of_visible[i,j] = Ts_inf[i, j]

    coeff = nphot / np.sum(sp_red)

    sp_red_n = np.asarray(sp_red * coeff, dtype=int)

    return [sp_red_n, map_of_visible]


## ++++
def spectra_any (Tmap, Rns, Mns, phase, chi, inc, eph):
    """
    |


    Calculate redshifted thermal spectra [erg s^-1 cm^-2 keV^-1] emitted by neutron star observed from arbitrary orientation.

    :math:`H_E^\\infty = \\frac{15 \\sigma_{SB}}{\\pi^5 k_B^4} \\int_\\mathrm{viz} \\frac{E^3 \\cos \\alpha \\; \\mathcal{D} \\sin \\theta d\\theta d\\phi}{\\exp(E/(k_B T_s^\\infty)) - 1}`

    :param Tmap: member of Tmap class containing surface temperature map.
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param phase: rotational phase [radian] 
    :param chi:  magnetic obliquity angle (angle between orientation of original dipolar magnetic field - top of the surface thermal map)
    :param inc: inclination of the observer with respect to the rotational axis.
    :param eph: list containing energy mesh :math:`E` [keV]

    :returns: :sp: redshifted spectra [erg s^-1 cm^-2 keV^-1]
              :visible_map: two-dimensional array which contains only the temperature distribution on visible hemisphere

    """

    map_of_visible = np.zeros ((len(Tmap.phi), len(Tmap.theta)))

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    Bp = 1e14         ## G
    c   = 2.998e+10   ## speed of light in cm/s

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c
    g14 = G*Mns*Msol / R / R / sqrt(1.0 - xg ) / 1e14 ## Gudmundsson et al. (1983), eq. (2-3)

    Ts_inf = Tmap.Ts * sqrt(1 - xg)

    ## Here we prepare variables for integration over the visible hemisphere

    dtheta = Tmap.theta[1] - Tmap.theta[0]
    dphi   = Tmap.phi[1] - Tmap.phi[0]

    theta1, phi1 = np.meshgrid (Tmap.theta, Tmap.phi)

    factor_int =   15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(theta1) * dtheta * dphi

    Dcosalpha = precompute_Dcos_alpha (Rns, Mns, chi, inc, phase, phi1, theta1)

    sp_red_en = np.zeros(len(eph))

    for i in range (0, len(eph)):

        sp_red_en[i] = np.sum (factor_int * Dcosalpha * np.power(eph[i], 3) / (np.exp(eph[i] / kB / Ts_inf) - 1.0))


    return [sp_red_en, Dcosalpha * Ts_inf]








def get_redshifted_spectra_any_photons (theta, phi, Tmap, Rns, Mns, phase, chi, inc, eph, nphot):
    """
    |

    Calculate thermal spectra [integer number of photons per energy bin] emitted by neutron star.

    :param Tmap: two dimensional array describing the surface thermal map [K]
    :param theta: list of magnetic latitude [radians] where Tmap is provided
    :param phi: list of magnetic longuitudets [radians] where Tmap is provided
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param phase: rotational phase [radian] 
    :param chi:  magnetic obliquity angle (angle between orientation of original dipolar magnetic field - top of the surface thermal map)
    :param inc: inclination of the observer with respect to the rotational axis.
    :param eph: list energies where spectra should be computed [keV]
    :param nphot: total number of photons to be generated

    :returns: :sp: number of photons per energy bin
              :visible_map: two-dimensional array which contains only the temperature distribution on visible hemisphere
    

    """

    map_of_visible = np.zeros ((len(phi), len(theta)))

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    Bp = 1e14         ## G
    c   = 2.998e+10   ## speed of light in cm/s

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c
    g14 = G*Mns*Msol / R / R / sqrt(1.0 - xg ) / 1e14 ## Gudmundsson et al. (1983), eq. (2-3)

    #print ('xg = ', xg, ' log_10 of  mean non-redshifted Ts = ', log10(np.mean(Tmap)))

    Ts_inf = Tmap * sqrt(1 - xg)

    ## Here we prepare variables for integration over the visible hemisphere

    sp_red = np.zeros(len(eph))

    dtheta = theta[1] - theta[0]
    dphi   = phi[1] - phi[0]

    theta1, phi1 = np.meshgrid (theta, phi)


## ==============================================
#    factor_int = sigma_SB * R * R * np.power(Ts_inf, 4.0) * np.sin(theta1) / pi * dtheta * dphi
#
#    res_int = []
#
#    for i in range (0, len(phases)):
#        Dcosalpha = precompute_Dcos2_alpha (Rns, Mns, chi, inc, phases[i], phi1, theta1)
#        res_int.append (np.sum(factor_int * Dcosalpha))

## ==============================================

    factor_int =   15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(theta1) * dtheta * dphi

    #np.power(eph, 2) / (np.exp(eph / kB / Ts_inf) - 1.0)

    #print ('Shape of factor_int is ', factor_int.shape)

    Dcosalpha = precompute_Dcos2_alpha (Rns, Mns, chi, inc, phase, phi1, theta1)

    #res_int = np.sum(factor_int * Dcosalpha)

    #res_int = Dcosalpha * 

    sp_red_en = np.zeros(len(eph))

    for i in range (0, len(eph)):

        sp_red_en[i] = np.sum (factor_int * Dcosalpha * np.power(eph[i], 2) / (np.exp(eph[i] / kB / Ts_inf) - 1.0))



#    for i in range (0, len(phi)):
#        for j in range (0, len(theta)):
#            new_theta = sin(theta[j])*cos(phi[i])
#            al = alpha (new_theta, Rns, Mns)



            #print (phi[i], theta[j], al)
#            if al < pi / 2.0:
#                Df = D_factor (new_theta, Rns, Mns)
#                sp_red = sp_red +  Df * 15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(theta[j]) * np.cos(al) * np.power(eph, 2) / (np.exp(eph / kB / Ts_inf[i,j]) - 1.0) * dtheta * dphi
#                map_of_visible[i,j] = Ts_inf[i, j]

    coeff = nphot / np.sum(sp_red_en)

    sp_red_n = np.asarray(sp_red_en * coeff, dtype=int)

    return [sp_red_n, Dcosalpha * Ts_inf]

def get_redshifted_phase_resolved_spectroscopy_photons (theta, phi, Tmap, Rns, Mns, phases, chi, inc, eph, nphot):
    """
    |

    Calculate thermal spectra [integer number of photons per energy bin] emitted by neutron star.

    :param Tmap: two dimensional array describing the surface thermal map [K]
    :param theta: list of magnetic latitude [radians] where Tmap is provided
    :param phi: list of magnetic longuitudets [radians] where Tmap is provided
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param phases: rotational phases [radian] 
    :param chi:  magnetic obliquity angle (angle between orientation of original dipolar magnetic field - top of the surface thermal map)
    :param inc: inclination of the observer with respect to the rotational axis.
    :param eph: list energies where spectra should be computed [keV]
    :param nphot: total number of photons to be generated

    :returns: :sp: number of photons per energy bin
    

    """

    sp_ph = []

    for i in range (0, len(phases)):
        res, vis = get_redshifted_spectra_any_photons (theta, phi, Tmap, Rns, Mns, phases[i], chi, inc, eph, nphot)
        sp_ph.append (res)

    sp_ph = np.asarray(sp_ph)

    return sp_ph





def two_BB (param, Teff, Rns, Mns):

    sc, sh, pc, ph = param

    sc = abs(sc)
    sh = abs(sh)
    pc = abs(pc)
    ph = abs(ph)

    eph, sp_Teff1 = single_BB (pc*Teff, Rns, Mns)
    eph, sp_Teff2 = single_BB (ph*Teff, Rns, Mns)

    spec = sc * sp_Teff1 + sh * sp_Teff2

    return [eph, spec]


#def chi2_2BB (param, Teff, Rns, Mns, spec):

#    eph, spec_synth = two_BB (param, Teff, Rns, Mns)

#    res = 0.0

#    spec_n  = np.asarray(spec) / np.max(spec)
#    spec_sn = np.asarray(spec_synth) / np.max(spec) 

#    for i in range (0, len(spec)):
#        if (spec[i] > np.max(spec) / 1e3):
#        if (spec[i] > 1):
#            res = res + (spec[i] - spec_synth[i])**2.0 / abs(spec[i])**2
#            res = res + (log(spec[i]) - log(spec_synth[i]))**2.0 / abs(log(spec_synth[i])) - better complete shape fit
#            res = res + (spec[i] - spec_synth[i])**2 / abs(spec[i])
#            res = res + (spec[i] - spec_synth[i])**2 / abs(spec_synth[i])
#            res = res + (spec_n[i] - spec_sn[i])**2 / abs(spec_sn[i])
#            res = res + (spec[i] - spec_synth[i])**2 / (spec[i] + spec_synth[i])           

#    res = res / 142.0
#    res = sqrt(res)

#    return res

def chi2_1BB (sc, Teff, Rns, Mns, spec, eph, nphot, L):
    """
    |

    Calculate chi^2 for synthetic spectra [counts per energy bin] composed of single blackbody.

    :param param: array of two parameters [sc, pc]
    :sc: fraction of total surface area covered with first hot spot
    :pc: temperature of first hot spot as a fraction of Teff 
    :param Teff: effective temperature of neutron star
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list of energies where spectra is to be calculated [keV]
    :param nphot: total number of received photons

    :returns: chi^2
    

    """

    spec_synth = single_BB_photons (sc, Teff, Rns, Mns, eph, nphot, L)

    res = 0.0

    for i in range (0, len(spec)):
        if (spec[i] >= 1):
            res = res + (spec[i] - spec_synth[i])**2 / (spec[i])


    return res

def Cstat_1BB (param, Teff, Rns, Mns, spec, eph, nphot, L):
    """
    |

    Calculate C-statistics for synthetic spectra [counts per energy bin] composed of single blackbody.

    :param sc: fraction of total surface area covered with hot spot
    :param Teff: effective temperature of hot spot
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list of energies where spectra is to be calculated [keV]
    :param nphot: total number of received photons
    :param L: luminosity of neutron star which we try to model

    :returns: C-stat
    

    """

#    sc, pc = param

#    R = Rns * 1e5

#    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.

#    spec_synth = single_BB_photons (pc*Teff, Rns, Mns, eph, nphot, integ=False)

#    Lcomp = sigma_SB * R*R * 4.0 * pi * sc * pow(pc*Teff, 4)

#    spec_synth = spec_synth * Lcomp / L
    
    #print ('Lcomp = ', Lcomp)

    spec_synth = examine_spectral_fit_1BB_photons (param, Teff, Rns, Mns, eph, nphot, L, integ=False)

    res = 0.0

    for i in range (0, len(spec)):
        if (spec[i] >= 1) and (spec_synth[i] > 0):
            res = res + 2.0 * ( spec_synth[i] - spec[i] + spec[i] * (log(spec[i]) - log(spec_synth[i])))
        elif (spec_synth[i] <= 0):
            res = res + np.inf


    return res





def chi2_2BB (param, Teff, Rns, Mns, spec, eph, nphot):
    """
    |

    Calculate chi^2 for synthetic spectra [number of photons per energy bin] composed of two blackbody.

    :param param: array of two parameters [sc, sh, pc, ph]
    :sc: fraction of total surface area covered with first hot spot
    :pc: temperature of first hot spot as a fraction of Teff 
    :param Teff: effective temperature of neutron star
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list of energies where spectra is to be calculated [keV]
    :param nphot: total number of received photons

    :returns: chi^2
    

    """


    spec_synth = two_BB_photons (param, Teff, Rns, Mns, eph, nphot)

    res = 0.0

    for i in range (0, len(spec)):
        if (spec[i] >= 1):
            res = res + (spec[i] - spec_synth[i])**2 / (spec[i])          

    return res

def Cstat_2BB (param, Teff, Rns, Mns, spec, eph, nphot, L):
    """
    |

    Calculate C-statistics for synthetic spectra [number of photons per energy bin] composed of two blackbody.

    :param param: array of two parameters [sc, sh, pc, ph]
    :sc: fraction of total surface area covered with first hot spot
    :pc: temperature of first hot spot as a fraction of Teff 
    :param Teff: effective temperature of neutron star
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param eph: list of energies where spectra is to be calculated [keV]
    :param nphot: total number of received photons
    :param L: luminosity of neutron star which we try to fit

    :returns: C-stat
    

    """

    #R = Rns * 1e5

    #sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.

    #spec_synth = two_BB_photons (param, Teff, Rns, Mns, eph, nphot, integ=False)

    #sc, sh, pc, ph = param

    #Lcomp = sigma_SB * R*R * 4.0 * pi * (sc * pow(pc*Teff, 4) + sh * pow(ph*Teff, 4))

    #Lcomp = compute_L_param (param, Teff, Rns, Mns)

    #print ('Lcomp = ', Lcomp)

    #spec_synth = spec_synth * Lcomp / L

    spec_synth = examine_spectral_fit_2BB_photons (param, Teff, Rns, Mns, eph, nphot, L, integ=False)

    res = 0.0

    for i in range (0, len(spec)):
        if (spec[i] >= 1) and (spec_synth[i] > 0):
            res = res + 2.0 * ( spec_synth[i] - spec[i] + spec[i] * (log(spec[i]) - log(spec_synth[i])))
        elif (spec[i] >= 1) and (spec_synth[i] <= 0):
            res = res + np.inf


    return res



def fit_spectral_model_chi2 (Teff, Rns, Mns, spec, eph, nphot):
    """
    |

    Fit spectra [spec] with two models: (1) single blackbody model and (2) sum of two blackbody model and choose the best model

    :param Teff: effective temperature of neutron star
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param spec: spectra which we try to fit
    :param eph: list of energies where spectra was calculated [keV]
    :param nphot: total number of received photons
    :param L: luminosity of neutron star which we try to fit

    :returns: C-stat
    

    """



    x2 = [0.2, 0.3, 0.9, 1.4]

    x1 = [1.3]

    res_1BB = minimize (chi2_1BB,     x1, method = 'Nelder-Mead',args=(Teff, Rns, Mns, spec, eph, nphot))
    res_2BB = minimize (chi2_2BB, x2, method = 'Nelder-Mead',args=(Teff, Rns, Mns, spec, eph, nphot))


   
    if res_2BB.fun < (res_1BB.fun - 2.7): ## 2BB model is significantly better

        return [res_2BB.x[0], res_2BB.x[1], res_2BB.x[2], res_2BB.x[3], res_1BB.fun, res_2BB.fun]

    else: ## 1BB model is as good

        return [0.5, 0.5, res_1BB.x[0], res_1BB.x[0], res_1BB.fun, res_2BB.fun]


def fit_spectral_model_Cstat (Teff, Rns, Mns, spec, eph, nphot, L):
    """
    |

    Fit spectra [spec] with two models: (1) single blackbody model and (2) sum of two blackbody model and choose the best model using C-statistics

    :param Teff: effective temperature of neutron star
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param spec: spectra which we try to fit
    :param eph: list of energies where spectra was calculated [keV]
    :param nphot: total number of received photons
    :param L: luminosity of neutron star which we try to fit

    :returns: list [s1, s2, p1, p2, Cstat BB1, Cstat BB2]
    

    """


    x2 = [0.2, 0.8, 0.9, 1.3]

    x1 = [0.7, 1.3]

    res_1BB = minimize (Cstat_1BB, x1, method = 'Nelder-Mead',args=(Teff, Rns, Mns, spec, eph, nphot, L))
    res_2BB = minimize (Cstat_2BB, x2, method = 'Nelder-Mead',args=(Teff, Rns, Mns, spec, eph, nphot, L))


   
    if res_2BB.fun < (res_1BB.fun - 2.7): ## 2BB model is significantly better

        return [res_2BB.x[0], res_2BB.x[1], res_2BB.x[2], res_2BB.x[3], res_1BB.fun, res_2BB.fun]

    else: ## 1BB model is as good

        return [res_1BB.x[0]/2, res_1BB.x[0]/2, res_1BB.x[1], res_1BB.x[1], res_1BB.fun, res_2BB.fun]






## Auxiliary function to compute the multiplication factors for lightcurve efficiently

def precompute_Dcos_alpha (Rns, Mns, chi, inc, phase, phi1, theta1):

    """
    |

    Calculate :math:`D \\cos\\alpha` factor. 

    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param chi:  magnetic obliquity angle (angle between orientation of original dipolar magnetic field - top of the surface thermal map)
    :param inc: inclination of the observer with respect to the rotational axis.
    :param phase: rotational phase [radian] 
    :param phi1: list of magnetic longuitudets [radians] where Tmap is provided
    :param theta1: list magnetic latitude [radians] where Tmap is provided

    :returns: list of D cos alpha factors
    

    """

    cospsi = (np.sin(chi)*np.cos(theta1) + np.sin(theta1)*np.cos(chi)*np.cos(phi1))*np.sin(inc)*np.cos(phase) + (-np.sin(chi)*np.sin(theta1)*np.cos(phi1) + np.cos(chi)*np.cos(theta1))*np.cos(inc) + np.sin(inc)*np.sin(phase)*np.sin(phi1)*np.sin(theta1)

    #print ('Comparison: cospsi = ', cospsi, ' cos theta = ', np.cos(theta1))



    Dfact = D_factor (cospsi, Rns, Mns)
    Dfact = np.nan_to_num (Dfact, copy = False)
    cos_al = cos_alpha (cospsi, Rns, Mns)
    cos_al = np.nan_to_num(cos_al, copy = False)
    b = (cos_al > 0).astype(int) ## masking invisible part of neutron star 
    b = np.nan_to_num (b, copy = False)

    return b*Dfact * cos_al 

## Auxiliary function to compute the multiplication factors for lightcurve efficiently taking into account cos^2 beaming

def precompute_Dcos2_alpha (Rns, Mns, chi, inc, phase, phi1, theta1):
    
    cospsi = (np.sin(chi)*np.cos(theta1) + np.sin(theta1)*np.cos(chi)*np.cos(phi1))*np.sin(inc)*np.cos(phase) + (-np.sin(chi)*np.sin(theta1)*np.cos(phi1) + np.cos(chi)*np.cos(theta1))*np.cos(inc) + np.sin(inc)*np.sin(phase)*np.sin(phi1)*np.sin(theta1)
    Dfact = D_factor (cospsi, Rns, Mns)
    Dfact = np.nan_to_num (Dfact, copy = False)
    cos_al = cos_alpha (cospsi, Rns, Mns)
    cos_al = np.nan_to_num(cos_al, copy = False)
    b = (cos_al > 0).astype(int)
    b = np.nan_to_num (b, copy = False)
    return b*Dfact * cos_al * cos_al * cos_al


## Efficient calculations of lightcurve - no beaming

def lightcurve (theta, phi, Tmap, Rns, Mns, phases, chi, inc):
    """
    |

    Calculate soft X-ray lightcurve

    :param Tmap: Tmap is the surface thermal map [K]
    :param theta: list magnetic latitude [radians] where Tmap is provided
    :param phi: list of magnetic longuitudets [radians] where Tmap is provided
    :param Rns: radius of neutron star [km]
    :param Mns: mass of neutron star [Solar mass]
    :param phases: list of phases [radian] 
    :param chi:  magnetic obliquity angle (angle between orientation of original dipolar magnetic field - top of the surface thermal map)
    :param inc: inclination of the observer with respect to the rotational axis.

    :returns: Soft X-ray lightcurve in units of intensity for each rotational phase [erg/s]
    

    """

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    c   = 2.998e+10   ## speed of light in cm/s

    theta1, phi1 = np.meshgrid (theta, phi)

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c

    Ts_inf = Tmap * sqrt(1 - xg)

    dtheta = theta[10] - theta[9]
    dphi   = phi[10] - phi[9]

    factor_int = sigma_SB * R * R * np.power(Ts_inf, 4.0) * np.sin(theta1) / pi * dtheta * dphi

    res_int = []

    for i in range (0, len(phases)):
        Dcosalpha = precompute_Dcos_alpha (Rns, Mns, chi, inc, phases[i], phi1, theta1)
        res_int.append (np.sum(factor_int * Dcosalpha))

    return res_int

def lightcurve_cos2 (theta, phi, Tmap, Rns, Mns, phases, chi, inc):

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    c   = 2.998e+10   ## speed of light in cm/s

    theta1, phi1 = np.meshgrid (theta, phi)

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c

    Ts_inf = Tmap * sqrt(1 - xg)

    dtheta = theta[10] - theta[9]
    dphi   = phi[10] - phi[9]

    factor_int = sigma_SB * R * R * np.power(Ts_inf, 4.0) * np.sin(theta1) / pi * dtheta * dphi

    res_int = []

    for i in range (0, len(phases)):
        Dcosalpha = precompute_Dcos2_alpha (Rns, Mns, chi, inc, phases[i], phi1, theta1)
        res_int.append (np.sum(factor_int * Dcosalpha))

    return res_int




