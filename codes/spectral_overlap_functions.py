import numpy as np
from scipy.interpolate import CubicSpline
import scipy.constants as const
from pathlib import Path
from datetime import datetime

def read_files(star, atmos, pigm):
    """Helper function to load spectra from
    specified files.
    
    Input:
    ------
    star:       path of file with spectral flux density spectrum 
    atmos:      path of file with atmosphere transmitance spectrum
    pigm:       path of file with extinction coefficient spectrum
    
    Returns:
    -------
    flux:     2D array with stellar spectral flux density spectrum 
                  [Angstrom; erg cm-2 s-1 A-1]
    trans:    2D array with atmosphere transmitance spectrum 
                  [microns; None (0-1)]
    epsilon:  2D array with pigment extinction coefficient spectrum
                  [nm; M-1 cm-1]
                  
    """
    
    flux = np.loadtxt(star,usecols=(0,1)) # [Angstrom, erg/cm2/s/A]
    trans = np.loadtxt(atmos,usecols=(0,1)) # [mum, None]
    epsilon = np.loadtxt(pigm,usecols=(0,1)) # [nm, M-1 cm-1]
    
    return flux, trans, epsilon

def get_exo_params(atmos, df_exo):
    """Helper function to extract spectral type (Stype) and habitable
    zone (hz) region from the atmosphere file name. With these parameters, 
    find the corresponding solar radius (Rs), exoplanet orbital 
    semi-major axis (sma), and Total Stellar Irradiance (S) from
    the auxiliary dataframe df_exo.
    
    Input:
    ------
    atmos:      File name with atmosphere transmitance spectrum
    df_exo:     Dataframe with host star/exoplanet parameters
    
    Returns:
    -------
    Stype:      String with stellar Spectral Type
    Rs:         Solar radius (in km)
    sma:        Exoplanet orbital semi-major axis (in km)
    S:          Total Stellar Irradiance (in W m-2)
    
    """
    
    # find indices corresponding to spectral type name
    i_first = atmos.index('_') + 1
    i_last = atmos.index('_hz')
    Stype = atmos[i_first:i_last]

    # find indices corresponding to habitable zone
    i_first = i_last + 1
    i_last = atmos.index('_os')
    hz = atmos[i_first:i_last]
    
    # get solar radius and exoplanet orbit (in km), and Stellar Irradiance
    Rs = 696340*df_exo.loc[np.where(df_exo.ST==Stype)[0][0], 'Rs']
    sma = 1.495978707e+8*df_exo.loc[np.where(df_exo.ST==Stype)[0][0], hz]
    
    S = df_exo.loc[np.where(df_exo.ST==Stype)[0][0], 'S_' + hz]
    
    return Stype, Rs, sma, S

def save_results(star, atmos, pigm, Rs, sma, S, spectra, rates):
    """Helper function to save results in an adequate format with
    all information regarding the files and exoplanet properties
    used in spectral absorption rate computation.
    
    Input:
    ------
    star:               file name with spectral flux density spectrum 
    atmos:              file name with atmosphere transmitance spectrum
    pigm:               file name with extinction coefficient spectrum
    Rs:                 stellar radius (in km or whatever units)
    sma:                Exoplanet orbital semi-major axis (same units as Rs)
    S:                  Total Stellar Irradiance (in W m-2)
    spectra             2D array with wavelength values and flux, atmospher, absorption,
                        and spectral absorption rate spectra
    rates:              1D array with :
                            total absorption rate [s-1] (gamma_t)
                            absorption rate at B band [s-1] (gamma_B)
                            absorption rate at Q band [s-1] (gamma_Q)
                            photosynthetic photon flux density [micromol photon m-2 s-1] (ppfd)                        
            
    """
    
    # datetime object containing current date and time
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    
    # define header
    header = ['##### EXOPHOT SPECTRAL OVERLAP FILE #####'  + '\n',
              '# Date: ' + dt_string + '\n',
              '# SED file: ' + star + '\n',
              '# atmosphere file: ' + atmos + '\n',
              '# photopigment file: ' + pigm + '\n',
              '# Stellar radius (Rs in km): ' + '{:.6e}'.format(Rs) + '\n',
              '# Exoplanet orbit semi-major axis (sma in km): ' + '{:.6e}'.format(sma) + '\n',
              '# Total absorption rate (gamma_t in s-1): ' + '{:.5e}'.format(rates[0]) + '\n',
              '# B band absorption rate (gamma_B in s-1): ' + '{:.5e}'.format(rates[1]) + '\n',
              '# Q band absorption rate (gamma_Q in s-1): ' + '{:.5e}'.format(rates[2]) + '\n',
              '# B-Q bands cut-off (in nm): ' + '{:.0f}'.format(500) + '\n',
              '# photosynthetic photon flux density [ppfd in micromol photon m-2 s-1]: ' + '{:.5e}'.format(rates[3]) + '\n',
              '# Total Stellar Irradiance (S in W m-2): ' + '{:.5e}'.format(S) + '\n',
              '# -------------------------------------------------' + '\n',
              '# Col 1: Wavelength (wl in Angstrom) ' + '\n',
              '# Col 2: Spectral Flux Density (F_l in J cm-2 s-1 A-1) ' + '\n',
              '# Col 3: Atmosphere transmitance (T, no units) ' + '\n',
              '# Col 4: Pigment absorption cross section (sigma_abs in cm-2) ' + '\n',
              '# Col 5: Spectral Absorption Rate (Gamma_lambda in A-1 s-1)' + '\n'
             ]

    # create output folder if doesn't exist
    output_folder = Path().absolute() / 'output'
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # set file name and path
    file_name = atmos[:-4] + '_' + pigm[:-4] + '.txt'
    path_name = output_folder / file_name

    # open and write header and spectral overlap
    with path_name.open("w") as f:
        f.writelines(header)
        for line in spectra:
            f.write('{:.2f}'.format(line[0]) + ' ' 
                    + '{:.5e}'.format(line[1]) + ' ' 
                    + '{:.5e}'.format(line[2]) + ' ' 
                    + '{:.5e}'.format(line[3]) + ' ' 
                    + '{:.5e}'.format(line[4]) + ' ' + '\n')
        f.close() 

        
#####################################################################################
#### FUNCTIONS TO PREPROCESS SPECTRA AND COMPUTE TOTAL AND SPECTRAL ABSORPTION RATE
#####################################################################################

def unit_converter(flux, trans, epsilon, Rs, sma):
    """Helper function to convert wavelength units to Angstroms, flux
    to J cm-2 s-1 A-1, and extinction coefficient to cross section. In addition, calculates
    flux at planet's top atmosphere (Multiply flux by (Rs/sma)^2).
    
    Inputs:
    -------
    flux:     2D array with stellar spectral flux density spectrum 
                  [Angstrom; erg cm-2 s-1 A-1]
    trans:    2D array with atmosphere transmitance spectrum 
                  [microns; None (0-1)]
    epsilon:  2D array with pigment extinction coefficient spectrum
                  [nm; M-1 cm-1]
    Rs:       stellar radius (in km or whatever units)
    sma:      Exoplanet orbital semi-major axis (same units as Rs)
    
    Returns:
    -------
    flux_u:   2D array with stellar spectral flux density spectrum 
                  [Angstrom; J cm-2 s-1 A-1]
    trans_u:  2D array with atmosphere transmitance spectrum 
                  [Angstrom; None (0-1)]
    sigma:    2D array with pigment absorption cross section spectrum
                  [Angstrom; cm-2]
                  
    """
    
    flux_u = flux.copy()
    trans_u = trans.copy()
    sigma = epsilon.copy()
        
    # from mum/nm to Angstrom
    trans_u[:,0] = 1e4*trans[:,0]
    sigma[:,0] = 10*epsilon[:,0]
    
    # flux conversion to J/cm2/s/A and multiplication by solid angle subtended by star from exoplanet, (Rs/sma)^2
    # https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas   
    flux_u[:,1] = flux[:,1]*1e-7*(Rs/sma)**2
    
    # from epsilon to cross section (cm2)
    sigma[:,1]=1000*np.log(10)*epsilon[:,1]/const.N_A
    
    return flux_u, trans_u, sigma


def wavelength_interpolate(flux, trans, sigma):
    """Helper function to interpolate wavelenghts values to 
    the combination of those of fluxes, transmitance and absorbance.
    Wavelength values truncated to those of the pigment range.
    
    Inputs:
    -------
    flux:     2D array with stellar spectral flux density spectrum 
                  [Angstrom; J cm-2 s-1 A-1]
    trans:    2D array with atmosphere transmitance spectrum 
                  [Angstrom; None (0-1)]
    sigma:    2D array with pigment absorption cross section spectrum
                  [Angstrom; cm-2]
    
    Returns:
    -------
    flux_int:   2D array with stellar spectral flux density spectrum (Interpolated) 
                  [Angstrom; J cm-2 s-1 A-1]
    trans_int:  2D array with atmosphere transmitance spectrum (Interpolated)
                  [Angstrom; None (0-1)]
    sigma_int:    2D array with pigment absorption cross section spectrum (Interpolated)
                  [Angstrom; cm-2]
                  
    """
    
    # truncate wavelength values to those of the pigment wavelength
    # range and remove repeated wavelength values (May happen due to
    # rounding errors, specially in M stars)
    # truncate
    wl_min, wl_max = sigma[0,0], sigma[-1,0]
    index_trun = (flux[:,0]>=wl_min) & (flux[:,0]<=wl_max)
    flux = flux[index_trun,:]
    index_trun = (trans[:,0]>=wl_min) & (trans[:,0]<=wl_max)
    trans = trans[index_trun,:]
    # remove repeated    
    index_unique = np.unique(flux[:,0], return_index = True) [1]
    flux = flux[index_unique,:]
    
    # store wavelength (wl) values
    wl_t = trans[:,0]
    wl_s = sigma[:,0]
    wl_f = flux[:,0]
    
    # create cubic spline functions of spectra
    cs_f = CubicSpline(wl_f, flux[:,1])
    cs_t = CubicSpline(wl_t, trans[:,1])
    cs_s = CubicSpline(wl_s, sigma[:,1])
    
    # concatenate all wavelength values and remove those repeated
    wl_int = np.unique(np.concatenate((wl_f,wl_t,wl_s)))
        
    # find spectral values at new wavelength values    
    f_int = cs_f(wl_int)
    t_int = cs_t(wl_int)
    s_int = cs_s(wl_int)
    
    # concatenate interpolated wavelengths and spectra
    flux_int = np.column_stack((wl_int,f_int))
    trans_int = np.column_stack((wl_int,t_int))
    sigma_int = np.column_stack((wl_int,s_int))
    
    return flux_int, trans_int, sigma_int

def ppfd_calculator(flux, trans):
    """Helper function to interpolate wavelenghts values to 
    the combination of those of fluxes and transmitance 
    to compute photosynthetic photon flux density at the Photosynthetic
    Active Region (between 400 and 700 nm).
    
    Inputs:
    -------
    flux:     2D array with stellar spectral flux density spectrum 
                  [Angstrom; J cm-2 s-1 A-1]
    trans:    2D array with atmosphere transmitance spectrum 
                  [Angstrom; None (0-1)]
    
    Returns:
    -------
    ppfd:   photosynthetic photon flux density [micromol photon m-2 s-1]
                  
    """
    
    # truncate wavelength values to those between 400 nm and 700 nm
    #  and remove repeated wavelength values (May happen due to
    # rounding errors, specially in M stars)
    # truncate
    wl_min, wl_max = 4000, 7000
    index_trun = (flux[:,0]>=wl_min) & (flux[:,0]<=wl_max)
    flux = flux[index_trun,:]
    index_trun = (trans[:,0]>=wl_min) & (trans[:,0]<=wl_max)
    trans = trans[index_trun,:]
    # remove repeated    
    index_unique = np.unique(flux[:,0], return_index = True) [1]
    flux = flux[index_unique,:]
    
    # store wavelength (wl) values
    wl_t = trans[:,0]
    wl_f = flux[:,0]
    
    # create cubic spline functions of spectra
    cs_f = CubicSpline(wl_f, flux[:,1])
    cs_t = CubicSpline(wl_t, trans[:,1])
    
    # concatenate all wavelength values and remove those repeated
    wl_int = np.unique(np.concatenate((wl_f,wl_t)))
    nu = const.c/(wl_int*1e-10) # photon energy nu=c/wl. convert units of wl to m!!
        
    # find spectral values at new wavelength values    
    f_int = cs_f(wl_int)
    t_int = cs_t(wl_int)
    
    # photosynthetic photon flux density (ppfd)
    phi_par = t_int*f_int/(const.h*nu)*1e10/const.N_A
    ppfd = np.trapz(phi_par,wl_int)
        
    return ppfd

def spectral_overlap(flux, trans, epsilon, Rs, sma):
    """Helper function to compute spectral absorption rate (Gamma_lambda),
    total absorption rate (rho_t), absorption rates at the Q and B band 
    (rho_Q, rho_B), and PAR flux density from stellar flux at exoplanet orbit, 
    atmosphere transmitance, and pigment absorption.
    
    Inputs:
    -------
    star:               file name with spectral flux density spectrum 
    atmos:              file name with atmosphere transmitance spectrum
    pigm:               file name with extinction coefficient spectrum
    Rs:                 stellar radius (in km or whatever units)
    sma:                Exoplanet orbital semi-major axis (same units as Rs)
    
    Returns:
    -------
    spectra:             2D array with:
                                Wavelength values (Interpolated) [Angstrom]
                                Stellar spectral flux density spectrum (Interpolated) [J cm-2 s-1 A-1]
                                Atmosphere transmitance spectrum (Interpolated) [None (0-1)]
                                Pigment absorption cross section spectrum (Interpolated) [cm-2]
                                Spectral absorption rate spectrum [s-1 A-1]
    gamma_t:             total absorption rate [s-1]
    gamma_B:             absorption rate at B band [s-1]
    gamma_Q:             absorption rate at Q band [s-1]
    ppfd:                photosynthetic photon flux density [micromol photon m-2 s-1]
    
    """
        
    # Convert units  
    flux_u, trans_u, sigma = unit_converter(flux, trans, epsilon, Rs, sma)
    
    # Interpolate spectra
    flux_int, trans_int, sigma_int = wavelength_interpolate(flux_u,trans_u,sigma)
    
    # calculate integrand and integrate    
    wl = flux_int[:,0]
    nu = const.c/(wl*1e-10) # photon energy nu=c/wl. convert units of wl to m!!
    spec_abs_rate = sigma_int[:,1]*trans_int[:,1]*flux_int[:,1]/(const.h*nu)
    
    cut_off = 5000 # 500 nm
    gamma_B = np.trapz(spec_abs_rate[wl<cut_off],wl[wl<cut_off])
    gamma_Q = np.trapz(spec_abs_rate[wl>=cut_off],wl[wl>=cut_off])
    gamma_t = gamma_B + gamma_Q
    
    # Calculate photosynthetic photon flux density (ppfd)
    ppfd = ppfd_calculator(flux_u,trans_u)
    
    # stack wavelength values and all relevant spectra to save
    spectra = np.column_stack((wl, flux_int[:,1], trans_int[:,1], sigma_int[:,1], spec_abs_rate))
    
    return spectra, [gamma_t, gamma_B, gamma_Q, ppfd]
