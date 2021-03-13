import numpy as np
from scipy.interpolate import CubicSpline
import scipy.constants as const
from pathlib import Path
from datetime import datetime


def unit_converter(flux, trans, epsilon, Rs, sma):
    """Helper function to convert wavelength units to Angstroms, flux
    to J cm-2 s-1 A-1, and extinction coefficient to cross section. In addition, calculates
    flux at planet (Multiply flux by (Rs/sma)^2)
    
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
    
    # flux conversion to J/cm2/s/A and multiplication by solid angle subtended by star from exoplanet, (Rs/sma)^22
    # https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas   
    flux_u[:,1] = flux[:,1]*1e-7*(Rs/sma)**2
    
    # from epsilon to cross section (cm2)
    sigma[:,1]=1000*np.log(10)*epsilon[:,1]/const.N_A
    
    return flux_u, trans_u, sigma


def wavelength_interpolate(flux,trans,sigma):
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
    
    # store wavelength (wl) values
    wl_f = flux[:,0]
    wl_t = trans[:,0]
    wl_s = sigma[:,0]
    
    # create cubic spline functions of spectra
    cs_f = CubicSpline(wl_f, flux[:,1])
    cs_t = CubicSpline(wl_t, trans[:,1])
    cs_s = CubicSpline(wl_s, sigma[:,1])
    
    # concatenate all wavelength values and remove those repeated
    wl_int = np.unique(np.concatenate((wl_f,wl_t,wl_s)))
    
    #truncate interpolated values to the pigment wavelength range
    wl_min, wl_max = wl_s[0], wl_s[-1]       
    wl_int = wl_int[(wl_int>=wl_min) & (wl_int<=wl_max)]
    
    # find spectral values at new wavelength values    
    f_int = cs_f(wl_int)
    t_int = cs_t(wl_int)
    s_int = cs_s(wl_int)
    
    # concatenate interpolated wavelengths and spectra
    flux_int = np.column_stack((wl_int,f_int))
    trans_int = np.column_stack((wl_int,t_int))
    sigma_int = np.column_stack((wl_int,s_int))
    
    return flux_int, trans_int, sigma_int


def spectral_overlap(flux, trans, epsilon, Rs, sma):
    """Helper function to compute spectral absorption rate and
    total absorption rate from stellar flux at exoplanet orbit, 
    atmosphere transmitance, and pigment absorption.
    
    Inputs:
    -------
    flux:     2D array with stellar spectral flux density spectrum (Interpolated)
                  [Angstrom; J cm-2 s-1 A-1]
    trans:    2D array with atmosphere transmitance spectrum (Interpolated) 
                  [Angstrom; None (0-1)]
    sigma:    2D array with pigment absorption cross section spectrum (Interpolated)
                  [Angstrom; cm-2]
    
    Returns:
    -------
    flux_int:            same as Input flux
    trans_int:           same as Input trans
    sigma_int:           same as Input sigma
    spec_abs_rate_int:   2D array with spectral absorption rate spectrum
                               [Angstrom; s-1 A-1]
    abs_rate:            total absorption rate [s-1]    
    
    """
        
    # Convert units  
    flux_u, trans_u, sigma = unit_converter(flux, trans, epsilon, Rs, sma)
    
    # Interpolate spectra
    flux_int, trans_int, sigma_int = wavelength_interpolate(flux_u,trans_u,sigma)
    
    # calculate integrand and integrate
    wl = flux_int[:,0]
    nu = const.c/(wl*1e-10) # photon energy nu=c/wl. convert units of wl to m!!
    spec_abs_rate = sigma_int[:,1]*trans_int[:,1]*flux_int[:,1]/(const.h*nu)
    abs_rate = np.trapz(spec_abs_rate,wl)
    
    # add wavelength values to overlapped spectrum
    spec_abs_rate_int = np.column_stack((wl,spec_abs_rate))
    
    return flux_int, trans_int, sigma_int, spec_abs_rate_int, abs_rate


def save_results(star, atmos, pigm, Rs, sma, spec_abs_rate_int, abs_rate):
    """Helper function to save results in an adequate format with
    all information regarding the files and exoplanet properties
    used in spectral absorption rate computation
    
    Input:
    ------
    star:               file name spectral flux density spectrum 
    atmos:              file name with atmosphere transmitance spectrum
    pigm:               file name with extinction coefficient spectrum
    Rs:                 stellar radius (in km or whatever units)
    sma:                Exoplanet orbital semi-major axis (same units as Rs)
    spec_abs_rate_int:  2D array with spectral absorption rate spectrum
                               [Angstrom; s-1 A-1]
    abs_rate:           total absorption rate [s-1]
            
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
             '# Stellar radius (Rs in km): ' + str(Rs) + '\n',
             '# Exoplanet semi-major axis (sma in km): ' + str(sma) + '\n',
             '# Absorption rate (s^-1): ' + '{:.5e}'.format(abs_rate) + '\n',
             '## Wavelength (Angstrom)  Spectral Absorption Rate (A^-1 s^-1)' + '\n']

    # create output folder if doesn't exist
    output_folder = Path().absolute() / 'output'
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # set file name and path
    file_name = atmos[:-4] + '_' + pigm[:-4] + '.txt'
    path_name = output_folder / file_name

    # open and write header and spectral overlap
    with path_name.open("w") as f:
        f.writelines(header)
        for line in spec_abs_rate_int:
            f.write('{:.2f}'.format(line[0]) + ' ' + '{:.5e}'.format(line[1]) + '\n')
        f.close()  

