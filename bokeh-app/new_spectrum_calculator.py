from os.path import join, dirname

from scipy import constants as scon
from scipy import interpolate 
import numpy as np

from kinematics import * 
from cosmo import *

# bubble spectrum
def h2_omega_phi(f, log10_T, log10_H_on_beta, log10_alpha, log10_eta, mod='Semi-analytic'):

    T = 10**log10_T
    H_on_beta = 10**log10_H_on_beta
    alpha = 10**log10_alpha
    eta = 10**log10_eta

    if mod=='Envelope':
        a = 3.
        b = 0.94
        c = 1.5
    
    elif mod=='Semi-analytic':
        a = 1.
        b = 2.61
        c = 1.5

    else:
        a = 3.
        b = 1.51
        c = 2.18

    # dilution coeff.
    R = 1.235 * 10**-5 * g_eq**(4/3) * g_star(log10_T)**(-1/3)
    
    # velocity factor  
    delta = (0.48 * v_w(eta, alpha)) / (1 + 5.3 * v_w(eta, alpha)**2 + 5 * v_w(eta, alpha)**4)

    # peak frequency today
    f_on_beta = (0.35) / (1+ 0.07 * v_w(eta, alpha) + 0.69 * v_w(eta, alpha)**4)
    f_p = 4.88 * 10**-8 * f_on_beta * H_on_beta**-1 * T * g_eq**(1/3) * g_star(log10_T)**(1/6)

    # spectral function
    S = (a + b)**c / (b * (f/f_p)**(-a/c) + a * (f/f_p)**(b/c))**c
        
    return R * delta * (k_phi(eta, alpha) * alpha / (1 + alpha))**2 * H_on_beta**2 * S

# sound wave spectrum
def h2_omega_sw(f, log10_T, log10_H_on_beta, log10_alpha, log10_eta):
    T = 10**log10_T
    H_on_beta = 10**log10_H_on_beta
    alpha = 10**log10_alpha
    eta = 10**log10_eta

    # dilution coeff.
    R = 8.5 * 10**-6 * (100/g_star(log10_T))**(1/3)
    
    # velocity factor  
    delta = v_w(eta, alpha)

    # peak frequency today
    f_p = 8.9 * 10**-8 * v_w(eta, alpha)**-1 * H_on_beta**-1 * T * (g_star(log10_T)/100) **(1/6)

    # spectral function 
    S = (f/f_p)**3 * (7 / (4 + 3* (f/f_p)**2))**(7/2)

    return R * delta * (k_sw(eta, alpha) * alpha / (1 + alpha))**2 * H_on_beta * S


# turbolence spectrum
def h2_omega_turb(f, log10_T, log10_H_on_beta, log10_alpha, log10_eta):
    T = 10**log10_T
    H_on_beta = 10**log10_H_on_beta
    alpha = 10**log10_alpha
    eta = 10**log10_eta
    
    # velocity factor  
    delta = 2.01 * 10**-1 * v_w(eta, alpha)

    # peak frequency today
    f_p = 27 * 10**-8 * v_w(eta, alpha)**-1 * H_on_beta**-1 * T * g_star(log10_T)**(1/6)

    # spectral function 
    h_star = (16.5*(1E-6)) * (T/100.) * ((g_star(log10_T)/100.)**(1/6))
    S = (f/f_p)**3 / ((1 + (f/f_p))**(11/3) * (1 + 8*np.pi*f/h_star))

    return 3.35 * 10**-4 * (100/g_star(log10_T))**1/3 * delta * (k_turb(eta, alpha) * alpha / (1 + alpha))**(3/2) * H_on_beta**1 * S
