import numpy as np 
from numpy import genfromtxt
from scipy import interpolate
from os.path import join, dirname


data_folder = join(dirname(__file__),'./data/')


# define the sound-wave efficiency
def k_sw(v, alpha):

    kappa_A = (v**(6/5)) * (6.9*alpha) / (1.36 - 0.037*np.sqrt(alpha) + alpha)
    kappa_B = (alpha**(2/5)) / (0.017 + (0.997+alpha)**(2/5))
    kappa_C = np.sqrt(alpha) / (0.135 + np.sqrt(0.98 + alpha))
    kappa_D = alpha / (0.73 + 0.083*np.sqrt(alpha) + alpha)
    delta_kappa = - 0.9 * np.log(np.sqrt(alpha) / (1+np.sqrt(alpha)))
    xi_J = (np.sqrt(2*alpha/3 + alpha**2) + np.sqrt(1/3)) / (1 + alpha)
    c_s = 1 / np.sqrt(3)

    if v < c_s:
        return (c_s**(11/5) * kappa_A * kappa_B) / ((c_s**(11/5) - v**(11/5))*kappa_B + v*(c_s**(6/5))*kappa_A)

    elif v > xi_J:
        numerator = ((xi_J-1)**3) * (xi_J**(5/2)) * (v**(-5/2)) * kappa_C * kappa_D
        denominator = ((xi_J-1)**3 - (v-1)**3) * (xi_J**(5/2)) * kappa_C + ((v-1)**3) * kappa_D

        return numerator / denominator

    else:
        return kappa_B + (v - c_s)*delta_kappa + (((v-c_s)/(xi_J-c_s))**3) * (kappa_C - kappa_B - (xi_J-c_s)*delta_kappa)