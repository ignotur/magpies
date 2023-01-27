from math import *
import numpy as np


def g14 (Rns, Mns):
    """
    NAME: g14

    PURPOSE: calculate the free fall acceleration at the surface of neutron star in units of g/1e14 cm/s/s

    INPUT: Rns - neutron star radius in units of km, Mns - neutron star mass in units of solar mass

    OUTPUT: free fall acceleration at the surface of neutron star in units of g/1e14 cm/s/s
    

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

## Function to compute the effective temperature

def compute_Teff (theta, phi, Rns, Tmap):
    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.

    Lns = compute_L (theta, phi, Rns, Tmap)
    return (pow(Lns / (4.0 * pi * sigma_SB * pow(Rns*1e5, 2)) , 1.0 / 4.0))

## Lensing factor following the article Poutanen (2020) A&A 640, A24 (2020)

def D_factor (cospsi, Rns, Mns):
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

def get_redshifted_spectra_pole_3D (theta, phi, Tmap, Rns, Mns):
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


            if al < pi / 2.0:
                Df = D_factor (cos(theta[j]), Rns, Mns)
                sp_red = sp_red +  Df * 15.0 * sigma_SB / ( pow(pi, 5) * pow(kB, 4)) * np.sin(theta[j]) * np.cos(al)  * np.power(Eph, 3) / (np.exp(Eph / kB / Ts_inf[i,j]) - 1.0) * dtheta * dphi
                map_of_visible[i,j] = Ts_inf[i, j]

    #Teff_inf = pow(10, 5.6) * sqrt(1 - xg)

    #res_1b =  en_red * 15.0*sigma_SB / (4.0 * pow(pi, 4) * pow(kB, 4)) * np.power(Eph, 3) / (np.exp(Eph / kB / Teff_inf) - 1.0)


    max_sp = np.max(sp_red)

    for i in range (0, len(sp_red)):
        if sp_red[i] < 1.0: #max_sp / 1e8:
            sp_red[i] = 0.0

    return [Eph, sp_red, map_of_visible]

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

def two_BB_diff (param, Teff, Rns, Mns, spec):

    eph, spec_synth = two_BB (param, Teff, Rns, Mns)

    res = 0.0

    for i in range (0, len(spec)):
        if (spec[i] > 1) and (spec_synth[i] > 1):
#            res = res + (spec[i] - spec_synth[i])**2.0 / abs(spec[i])**2
            res = res + (log(spec[i]) - log(spec_synth[i]))**2.0 / abs(log(spec_synth[i]))

    return res

## Auxiliary function to compute the multiplication factors for lightcurve efficiently

def precompute_Dcos_alpha (Rns, Mns, chi, inc, phase, phi1, theta1):
    
    cospsi = (np.sin(chi)*np.cos(theta1) + np.sin(theta1)*np.cos(chi)*np.cos(phi1))*np.sin(inc)*np.cos(phase) + (-np.sin(chi)*np.sin(theta1)*np.cos(phi1) + np.cos(chi)*np.cos(theta1))*np.cos(inc) + np.sin(inc)*np.sin(phase)*np.sin(phi1)*np.sin(theta1)
    Dfact = D_factor (cospsi, Rns, Mns)
    cos_al = cos_alpha (cospsi, Rns, Mns)
    b = (cos_al > 0).astype(int)
    return b*Dfact * cos_al

## Efficient calculations of lightcurve - no beaming

def lightcurve (theta, phi, Tmap, Rns, Mns, phases, chi, inc):

    sigma_SB = 5.670e-5 ## erg⋅cm^{−2}⋅s^{−1}⋅K^{−4}.
    kB = 8.617e-8       ## keV / K
    G = 6.67430e-8    ## cgs
    Msol = 2e33       ## Solar mass in gramms
    c   = 2.998e+10   ## speed of light in cm/s

    theta1, phi1 = np.meshgrid (theta, phi)

    R = Rns * 1e5

    xg = 2.0 * G * Mns*Msol / R / c / c

    Ts_inf = Tmap * sqrt(1 - xg)

    dtheta = theta[1] - theta[0]
    dphi   = phi[1] - phi[0]

    factor_int = sigma_SB * R * R * np.power(Ts_inf, 4.0) * np.sin(theta1) / pi * dtheta * dphi

    res_int = []

    for i in range (0, len(phases)):
        Dcosalpha = precompute_Dcos_alpha (Rns, Mns, chi, inc, phases[i], phi1, theta1)
        res_int.append (np.sum(factor_int * Dcosalpha))

    return res_int




