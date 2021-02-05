import numpy as np
from scipy import constants as scon
from numpy import genfromtxt
from scipy import interpolate 

#import new_runaway
#from new_runaway_h2_omega import pt_h2_omega

# get deflagration velocity
#@function
def get_def_xi(eta, alpha, def_eta_tab='./data/def_eta_tab.csv', def_alpha_tab='./data/def_alpha_tab.csv', def_xi_tab='./data/def_xi_tab.csv'):
    
    eta_tab = genfromtxt(def_eta_tab, delimiter=',')
    alpha_tab = genfromtxt(def_alpha_tab, delimiter=',')
    xi_tab = genfromtxt(def_xi_tab, delimiter=',')

    xi_interp = interpolate.interp2d(eta_tab, alpha_tab, xi_tab, kind='cubic')

    return xi_interp(eta, alpha)

# get detonation velocity
#@function
def get_det_xi(eta, alpha, det_eta_tab='./data/det_eta_tab.csv', det_alpha_tab='./data/det_alpha_tab.csv', det_xi_tab='./data/det_xi_tab.csv'):
    
    eta_tab = genfromtxt(det_eta_tab, delimiter=',')
    alpha_tab = genfromtxt(det_alpha_tab, delimiter=',')
    xi_tab = genfromtxt(det_xi_tab, delimiter=',')

    xi_interp = interpolate.interp2d(eta_tab, alpha_tab, xi_tab, kind='cubic')

    return xi_interp(eta, alpha)
   
# deflagration if alpha < get_boundary(eta)
# detonation if get_boundary(eta) < alpha < alpha_inf
# runaway if alpha > alpha_inf
#@function 
def get_boundary(eta, table='./data/det2def_trans.csv'):
    
    # deflagration if alpha < get_boundary(eta)
    # detonation if get_boundary(eta) < alpha < alpha_inf
    # runaway if alpha > alpha_inf
    
    det2def_tab = genfromtxt(table, delimiter=',')
    
    det2def_interp = interpolate.interp1d(det2def_tab[::,0], det2def_tab[::,1], kind='cubic')
    
    return det2def_interp(eta)

# get g_star
#@function
def get_g_star(log10_T, g_star_path):
    
    # return the relativistic degree of freedom g_* for a given temperature T
    # log10_T: log10(T in GeV)
    # g_star_path: path to the table containing a column of T (in MeV) and a column of g_*
    # return: g_*
    
    # read the g_star table
    log10_T_table = np.log10(np.loadtxt(g_star_path)[:,0] / (1E3))
    log10_g_star_table = np.log10(np.loadtxt(g_star_path)[:,1])
    
    # compute g_star at T by interpolation
    log10_g_star = np.interp(log10_T, log10_T_table, log10_g_star_table)
    
    return 10**log10_g_star

# get the efficiency coefficient  
#@function
def get_nonrunaway_kappa(v, alpha):

    # efficiency factors
    kappa_A = (v**(6/5)) * (6.9*alpha) / (1.36 - 0.037*np.sqrt(alpha) + alpha)
    kappa_B = (alpha**(2/5)) / (0.017 + (0.997+alpha)**(2/5))
    kappa_C = np.sqrt(alpha) / (0.135 + np.sqrt(0.98 + alpha))
    kappa_D = alpha / (0.73 + 0.083*np.sqrt(alpha) + alpha)
    delta_kappa = - 0.9 * np.log(np.sqrt(alpha) / (1+np.sqrt(alpha)))
    xi_J = (np.sqrt(2*alpha/3 + alpha**2) + np.sqrt(1/3)) / (1 + alpha)
    c_s = 1 / np.sqrt(3)

    if v < c_s:

        kappa = (c_s**(11/5) * kappa_A * kappa_B) / ((c_s**(11/5) - v**(11/5))*kappa_B + v*(c_s**(6/5))*kappa_A)

    elif v > xi_J:

        numerator = ((xi_J-1)**3) * (xi_J**(5/2)) * (v**(-5/2)) * kappa_C * kappa_D
        denominator = ((xi_J-1)**3 - (v-1)**3) * (xi_J**(5/2)) * kappa_C + ((v-1)**3) * kappa_D
        kappa = numerator / denominator
    
    else:

        kappa = kappa_B + (v - c_s)*delta_kappa + (((v-c_s)/(xi_J-c_s))**3) * (kappa_C - kappa_B - (xi_J-c_s)*delta_kappa)

    return kappa

# get the type of phase transition
#@function
def pt_type(eta, alpha):
            
    # alpha_inf
    alpha_inf = eta + (1/3) * (1 - 0.85)
    
    if alpha > alpha_inf:
        
        return 'runaway'
    
    elif alpha < get_boundary(eta):
        
        return 'deflagration'
    
    else:
        
        return 'detonation'

    
    # main function to get the PSD of phase transition      
#@function
def pt_h2_omega(f, log10_T_star=None, log10_H_star_on_beta=None, log10_alpha_star=None, log10_eta=None, g_star_path='./data/g_star.txt', h=0.674, zp=6.9, gamma=4/3, epsilon_turb=0.1, components=2):
    
    # PSD of gravitational wave background from runaway and nonrunaway phase transition
    # the signal shapes come from Weir (2018)
    # only sound wave and turbulence contributions are considered
    # no spatial correlation
    
    # f: frequency (Hz)
    # log10_T_star: log10(T_*) where T_* is the temperature during PT (GeV)
    # log10_H_star_on_beta: log10(H_*/beta), where beta is the nucleation rate and H_* is the Hubble rate during PT. H_*/beta should be between 0.01 to 1
    # log10_alpha_star: log10(alpha_*), where alpha_* is the phase transition strength
    # log10_v: log10(v_w), where v_w is the wall velocity
    # g_star_path: path for the table with a column of T (in MeV) and a column of g_star
    # h: 100*h km/s/Mpc is the Hubble constant today, CURRENTLY NOT SUPPORTED
    # zp: parameter fitted to numerical simulations
    # gamma: adiabatic index for sound wave contribution
    # epsilon_turb: efficiency coefficient for turbulence
    # alpha_inf: threshould between runaway and nonrunaway phase transition
    
    # useful formulas:
        # h_c(f) = \sqrt(3/(2*pi^2)) * ((100 km/s/Mpc)/f) * \sqrt(h^2*\Omega(f))
        # S(f) = (1/(12*pi^2)) * h_c^2(f) / f^3

    df = np.diff(np.concatenate((np.array([0]), f[::components])))
    
    # compute the parameters
    eta = 10**log10_eta
    alpha = 10**log10_alpha_star
    T_star = 10**log10_T_star
    H_star_on_beta = 10**log10_H_star_on_beta

    # g_* during PT
    g_star = get_g_star(log10_T_star, g_star_path)

    # reference frequency (Hz)
    f_yr = 1/scon.Julian_year
    
    # Hubble rate today / h (100 km/s/Mpc)
    h0 = 3.24078E-18 # Hz

    # alpha_inf
    alpha_inf = eta + (1/3) * (1 - 0.85)
    
    ########################################################################################
    
    # type of phase transition
    phase_transition_type = pt_type(eta, alpha)

    # runaway
    if phase_transition_type == 'runaway':

        # velocity
        v = 1.

        # efficiency factor
        kappa_inf = alpha_inf / (0.73 + 0.083*np.sqrt(alpha_inf) + alpha_inf)
        kappa_phi = 1 - alpha_inf / alpha
        kappa_f = kappa_inf * alpha_inf / alpha

    # nonrunaway
    else:
        
        kappa_phi = 0.
        
        # deflagration
        if phase_transition_type == 'deflagration':
            
            v = get_def_xi(eta, alpha)[0]
            kappa_f = get_nonrunaway_kappa(v, alpha)
        
        # detonation
        elif phase_transition_type == 'detonation':
            
            v = get_det_xi(eta, alpha)[0]
            kappa_f = get_nonrunaway_kappa(v, alpha)
            
        else:
            
            raise ValueError('Phase transition type must be runaway, deflagration or detonation')
    

    ##########################################################################################
    
    # scalar field
        
    # parameters from U(1) in Table 1 in Lewicki and Vaskonen (July 2020)
    big_A = 3.10 / 100
    omega_bar_on_beta = 0.64
    a = 1.00
    b = 2.61
    c = 1.5
    
    # Hubble parameter during phase transition (Hz) (h = 0.674)
    h_star = (1.65e-5) * (T_star/100.) * ((g_star/100.)**(1/6))
    
    # peak freqeuncy (Hz)
    omega_0 = h_star * (1/H_star_on_beta) * omega_bar_on_beta
    
    # angular frequency
    omega = 2 * np.pi * f 
    
    # S_fit(omega)
    s_fit = (big_A * ((a+b)**c)) / (b*((omega/omega_0)**(-a/c)) + a*((omega/omega_0)**(b/c)))**c
    
    # h^2 * \Omega(f)
    h2_omega_env = (1.67E-5) * ((H_star_on_beta)**2) * ((kappa_phi*alpha/(1+alpha))**2) * ((100/g_star)**(1/3)) * s_fit

    ######################################################################################################
    
    # sound wave
        
    # rms fluid velocity ^ 2
    u_bar_f_2 = (3/4) * kappa_f * alpha/(1+alpha)
    
    # peak frequency (Hz)
    f_sw = (8.9*(1E-6)) * (1/v) * (1/H_star_on_beta) * (zp/10.) * (T_star/100.) * ((g_star/100.)**(1/6))
    
    # spectral form
    spectral_shape_sw = ((f/f_sw)**3) * ((7/(4+3*np.square(f/f_sw)))**(7/2))
    
    # h^2 * \Omega(f)
    h2_omega_sw = (8.5E-6) * ((100/g_star)**(1/3)) * (gamma**2) * (u_bar_f_2**2) * (H_star_on_beta) * v * spectral_shape_sw

    #####################################################################################################################
    
    # turbulence
    
    # kappa
    kappa_turb = kappa_f * epsilon_turb
    
    # peak frequency (Hz)
    f_turb = (27*(1E-6)) * (1/v) * (1/H_star_on_beta) * (T_star/100.) * ((g_star/100)**(1/6)) 
    
    # h_*
    h_star = (16.5*(1E-6)) * (T_star/100.) * ((g_star/100.)**(1/6))
    
    # spectral shape
    spectral_shape_turb = ((f/f_turb)**3) / (((1+f/f_turb)**(11/3)) * (1+8*np.pi*f/h_star))
    
    # h^2 * \Omega(f)
    h2_omega_turb = (3.35E-4) * (H_star_on_beta) * ((kappa_turb*alpha/(1+alpha))**(3/2)) * ((100./g_star)**(1/3)) * v * spectral_shape_turb
    
    ##################################################################################
    
    # total relic density
    h2_omega_sum = h2_omega_env + h2_omega_sw + h2_omega_turb
    
    return h2_omega_env, h2_omega_sw, h2_omega_turb, h2_omega_sum, kappa_f, kappa_phi, v


# spectrum from SMBHB

def smbhb_h2_omega(f, log10_A, gamma):
    
    h0 = 3.24078e-18 # Hz
    
    alpha = (3 - gamma) / 2
    A = 10 ** log10_A
    f_yr = 1 / scon.Julian_year # Hz
    
    hcf = A * ((f / f_yr) ** alpha)
    h2_omega = (2*np.pi*np.pi/3) * ((f/h0)**2) * (hcf**2)
    
    return h2_omega

# check if the PT is in the runaway regime 

def runaway_Q(log10_alpha, log10_eta):
    
    alpha = 10 ** log10_alpha
    eta = 10 ** log10_eta
    
    # alpha infinity
    alpha_inf = eta + (1/3) * (1-0.85)
    
    # the phase transition is runaway if alpha > alpha_inf
    if alpha > alpha_inf:
        return True
    
    else:
        return False



# frequency (Hz)
f = np.logspace(-11., np.log10(3e-8))

# h^2 * \Omega(f)
#pt = [pt_h2_omega(f, log10_T_list[i], log10_H_star_on_beta_list[i], log10_alpha_list[i], log10_eta_list[i])[-1] for i in range(len(log10_T_list))]

# T span
T_span = 12.5*scon.year # s

# first 5 frequency bins
freq_bins = np.loadtxt('data/nano12_freq.txt')

# median of log10(residuals (s))
log10_res_median = np.loadtxt('./data/nano12_median.txt')

# upper limit of log10(residuals (s))
log10_res_up = np.loadtxt('./data/nano12_up.txt')

# lower limit of log10(resuduals (s))
log10_res_down = np.loadtxt('./data/nano12_down.txt')

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