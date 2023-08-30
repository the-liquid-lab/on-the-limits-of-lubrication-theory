# -*- coding: utf-8 -*-
"""
@author: Clément & Arnaud
"""

#This file contain basic functions and define the equations nedded for 
#analytical resolution

############################## Import #########################################
import numpy as np
from mpmath import mp, cosh, sinh, tanh, exp, sqrt, findroot, j
import time
from scipy.optimize import curve_fit

#The package must be installed through "conda install gwr_inversion"
from gwr_inversion import gwr

######################## Basic signal evolutions ##############################
def decaying_sinusoid(t, om_dec, om_osc, amp = 1., phi = 0.):
    return amp*np.exp(- om_dec * t)*np.cos(om_osc * t + phi)

def decaying_sinusoid_zero_initial_velocity(t, om_dec, om_osc, amp = 1., phi = 0.):
    return amp*np.exp(- om_dec * t)*(np.cos(om_osc * t + phi) + np.sin(om_osc * t + phi) * om_dec / om_osc)

################# Basic length of Cortelezzi resolution #######################
def ker_sy (s, Oh, Bo, k, lbda):
    return 2*Oh/s*k*(k-lbda*tanh(k)) - Oh/s*(4*lbda*k*sinh(k)*(k*exp(-lbda)
            *(k*cosh(k)+lbda*sinh(k))-(k**2+lbda**2))+(k**2+lbda**2)**2
            *sinh(lbda))/(2*k*cosh(k)*(k*cosh(k)*sinh(lbda)-lbda*sinh(k)*cosh(lbda)))
            
def eta_sy (s, Oh, k, omega2, Kern):
    return 1/s*(1-omega2/(s**2+4*Oh*k**2*s+omega2+2*Oh*k**2*s*Kern))

######################## Reduction of the expression ##########################
#Reduce the expressions as functions of s and of the parameters Oh, Bo and k
def freeSurfaceLaplace(s, Oh, Bo, k):
    lbda = sqrt(k**2 + s/Oh)
    omega2 = (Bo+k**2)*k*tanh(k)
    ker = ker_sy (s, Oh, Bo, k, lbda)
    return eta_sy(s, Oh, k, omega2, ker)

def denom (s, Oh, Bo, k):
    lbda = sqrt(k**2 + s/Oh)
    omega2 = (Bo+k**2)*k*tanh(k)
    ker = ker_sy (s, Oh, Bo, k, lbda)
    return (s**2+4*Oh*k**2*s+omega2+2*Oh*k**2*s*ker)

############################## Laplace Inversion ##############################
#Inverse the Laplace transfrom and return the values of eta as a function 
#of a range of t and the parameters Oh, Bo and k
#This resolution uses the gwr packageand need mp package for variable definitions.
#The M_value is an important parameter for precision.
def freeSurface(t_all, Ohnumb, Bonumb, knumb, M_value = 32):
    store = time.time()
    Oh = mp.mpmathify(Ohnumb)
    Bo = mp.mpmathify(Bonumb)
    k = mp.mpmathify(knumb) 
    f = lambda s: freeSurfaceLaplace(s, Oh, Bo, k)
    a = [float(gwr(f, t, M_value)) for t in t_all]
    print (time.time()-store)
    return a


############ Expression of the growth rates and pulsations ####################
###Purely analytical expressions
#Lubrication growth rate
def om_lub(Oh, Bo, k):
    return (k**2*Bo+k**4)/(3*Oh)

#Inertial pulsation
def pulsation(Bo, k):
    return np.sqrt(np.abs(Bo + k**2)*k*np.tanh(k))

#Asymptotic solutions obtained from the normal mode in Cortelezzi's derivation
def om_normal_mode_viscous(Oh, Bo, k):
    return -pulsation(Bo, k)**2/(k**2*Oh*np.tanh(k))*(k-np.cosh(k)*np.sinh(k))/(1+2*k**2+np.cosh(2*k))
    
def puls_normal_mode_inertial(Oh, Bo, k):
    return pulsation(Bo, k) - (1/np.sinh(2*k)*np.sqrt(pulsation(Bo, k) * k**2*Oh/2)
            - pow(k**2*Oh,3./2.)/np.sqrt(2*pulsation(Bo, k))
            *(3-8*np.cosh(2*k)-14*np.cosh(4*k)+4*np.cosh(6*k))/(8*np.sinh(2*k)**3)) 

def om_normal_mode_inertial(Oh, Bo, k):
    if Bo > 0:
        return (1/np.sinh(2*k)*np.sqrt(pulsation(Bo, k) * k**2*Oh/2) +
                2*k**2*Oh * (np.cosh(4*k)+np.cosh(2*k)-1) / (np.cosh(4*k) -1)
                - pow(k**2*Oh,3./2.)/np.sqrt(2*pulsation(Bo, k))
                *(3-8*np.cosh(2*k)-14*np.cosh(4*k)+4*np.cosh(6*k))/(8*np.sinh(2*k)**3))
    else:
        return pulsation(Bo, k) - (1/np.sinh(2*k)*np.sqrt(pulsation(Bo, k) * k**2*Oh/2)
                + 2*k**2*Oh * (np.cosh(4*k)+np.cosh(2*k)-1) / (np.cosh(4*k) -1)
                + pow(k**2*Oh,3./2.)/np.sqrt(2*pulsation(Bo, k))
                *(3-8*np.cosh(2*k)-14*np.cosh(4*k)+4*np.cosh(6*k))/(8*np.sinh(2*k)**3))

###Resolution by normal mode analysis (zero of the denominator)
def om_analytic(Oh, Bo, k, guess = False):
    if guess:
        root_denom = findroot(lambda s: denom (s, Oh, Bo, k), guess)
    else: 
        try:
            root_denom = findroot(lambda s: denom (s, Oh, Bo, k), om_lub(Oh, Bo, k))
        except ValueError:
            root_denom = findroot(lambda s: denom (s, Oh, Bo, k), j*pulsation(Bo, k))
    return root_denom

#Growth rate and pulsations obtained by fit of the numerical solution.
def om_numerical(Oh, Bo, k, guess_value):
    M = 64
    om_relax = guess_value[0]
    om_0 = guess_value[1]
    logspan = np.array([1e-4,2e-4,5e-4,1e-3,2e-3,5e-3])
    linspan = np.linspace(0.01, 1., 50)
    loglinspan = np.concatenate((logspan,linspan))
    t_all = loglinspan * 10./max(om_0, abs(om_relax))
    sampled_eta = freeSurface(t_all, Oh, Bo, k, M)
    guess_value = list(guess_value)
    guess_value[3] = 7.*np.pi/4.
    guess_value = tuple(guess_value)
    popt = curve_fit(decaying_sinusoid, t_all, sampled_eta, p0=guess_value, bounds=([0.,0.,1.,0.],[np.inf, 5.*om_0, 2., 2.*np.pi]), sigma=(1.+10.*np.exp(-om_relax*t_all)))[0]
    return popt, t_all, sampled_eta