from os.path import join, dirname

from scipy import constants as scon
from scipy import interpolate 
import numpy as np

from kinematics import * 
from cosmo import *

# bubble spectrum
def h2_omega(f, log10_T, log10_H_on_beta, log10_alpha, vw, contr='bubble', mod='Semi-analytic'):

    [T, H_on_beta, alpha, vw] = 10**np.array([log10_T, log10_H_on_beta, log10_alpha, np.log10(vw)])

    gs = g_star(log10_T)

    ############
    ## bubble ##
    ############fid
    if contr == 'bubble':
        if mod=='Envelope':
            a = 3.
            b = 0.94
            c = 1.5
            # peak frequency at emission 
            f_on_beta = 0.35/(1+0.07 * vw + 0.69 * vw**2)
        
        elif mod=='Semi-analytic':
            a = 1.
            b = 2.61
            c = 1.5
            # peak frequency at emission 
            f_on_beta = 0.1

        else:
            a = 0.7
            b = 2.3
            c = 1
            f_on_beta = 0.2

    
        # velocity factor  
        delta = 0.48 * vw**3 / (1 + 5.3 * vw**2 + 5 * vw**4)

        # efficiency factor
        k = 1

        [p, q] = [2 ,2]

        # spectral shape 
        def S(x):
            return (a + b)**c / (b * x**(-a/c) + a * x**(b/c))**c


    ############
    ## sound  ##
    ############
    if contr == 'sound':

        k = k_sw(vw, alpha)

        cs = 1 / np.sqrt(3)
        ups = 1 - (1 + 4 * (8 * np.pi)**1/3 * max(vw, cs) * H_on_beta * (1 + alpha)**1/2 / (3 * k * alpha)**1/2)**-1/2

        delta = 5.13 * 10**-1 * vw * ups


        [p, q] = [2 ,1]
    
        def S(x):
            return x**3 * (7/(4 + 3 * x**2))**(7/2)

        f_on_beta = 5.36 * 10**-1 / vw


    # dilution coeff.
    R = 7.69 * 10**-5 * gs**(-1/3)
    
    # peak frequency today in Hz
    f0 = 1.13 * 10**-7 * f_on_beta * H_on_beta**-1 * T * (gs/10)**(1/6)
        
    return R * delta * (k * alpha / (1 + alpha))**p * H_on_beta**q * S(f/f0)

