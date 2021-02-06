import numpy as np 
from scipy import constants as scon
from os.path import join, dirname

# T span
T_span = 12.5*scon.year # s

# first 5 frequency bins

freq_bins = np.loadtxt(join(dirname(__file__), 'data/nano12_freq.txt'))

# median of log10(residuals (s))
log10_res_median = np.loadtxt(join(dirname(__file__), 'data/nano12_median.txt'))

# upper limit of log10(residuals (s))
log10_res_up = np.loadtxt(join(dirname(__file__), 'data/nano12_up.txt'))

# lower limit of log10(resuduals (s))
log10_res_down = np.loadtxt(join(dirname(__file__), 'data/nano12_down.txt'))

# NanoGrav
h0 = 2.27E-18 # Hz
f_yr = 1/scon.year # Hz

# PSD(f)
psd_median = ((10**log10_res_median)**2) * T_span # s^3
psd_up = ((10**log10_res_up)**2) * T_span # s^3
psd_down = ((10**log10_res_down)**2) * T_span # s^3

# h_c(f)^2
hcf2_median = 12*np.pi*np.pi*psd_median*(freq_bins**3) # dimensionless
hcf2_up = 12*np.pi*np.pi*psd_up*(freq_bins**3) # dimensionless
hcf2_down = 12*np.pi*np.pi*psd_down*(freq_bins**3) # dimensionless

# h^2*omega(f)
h100 = 3.24E-18 # Hz
h2omega_median = 2*np.pi*np.pi/(3*np.square(h100)) * (freq_bins**2) * hcf2_median # dimensionless
h2omega_up = 2*np.pi*np.pi/(3*np.square(h100)) * (freq_bins**2) * hcf2_up # dimensionless
h2omega_down = 2*np.pi*np.pi/(3*np.square(h100)) * (freq_bins**2) * hcf2_down # dimensionless

yerr_up = h2omega_up - h2omega_median
yerr_down = h2omega_median - h2omega_down
yerr = [yerr_down[:5], yerr_up[:5]]