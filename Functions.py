# -*- coding: utf-8 -*-
"""
@author: Clément & Arnaud
"""
import numpy as np
from mpmath import mp, cosh, sinh, tanh, exp, sqrt, findroot, j
import time

#The package must be installed through "conda install gwr_inversion"
from gwr_inversion import gwr

## Functions and expressions declarations
def decaying_sinusoid(t, om_dec, om_osc):
    return np.exp(- om_dec * t)*np.cos(om_osc * t)

def better_sinusoid(t, om_dec, om_osc, amp, phi):
    return amp*np.exp(- om_dec * t)*np.cos(om_osc * t + phi)

def my_exp(t, om_dec):
      return np.exp(- om_dec * t)
  
#Declare the expressions of the kernel and eta
def ker_sy (s, Oh, Bo, k, lbda):
    return 2*Oh/s*k*(k-lbda*tanh(k)) - Oh/s*(4*lbda*k*sinh(k)*(k*exp(-lbda)
            *(k*cosh(k)+lbda*sinh(k))-(k**2+lbda**2))+(k**2+lbda**2)**2
            *sinh(lbda))/(2*k*cosh(k)*(k*cosh(k)*sinh(lbda)-lbda*sinh(k)*cosh(lbda)))
            
def eta_sy (s, Oh, k, omega2, Kern):
    return 1/s*(1-omega2/(s**2+4*Oh*k**2*s+omega2+2*Oh*k**2*s*Kern))

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

#Inverse the Laplace transfrom and return the values of eta as a function 
#of a range of t and the parameters Oh, Bo and k
def freeSurface(t_all, Ohnumb, Bonumb, knumb, M_value = 32):
    store = time.time()
    Oh = mp.mpmathify(Ohnumb)
    Bo = mp.mpmathify(Bonumb)
    k = mp.mpmathify(knumb) 
    f = lambda s: freeSurfaceLaplace(s, Oh, Bo, k)
    a = [float(gwr(f, t, M_value)) for t in t_all]
    print (time.time()-store)
    return a

#Calculation of the different growth rates and pulsations
def om_lub(Oh, Bo, k):
    return (k**2*Bo+k**4)/(3*Oh)

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
    return (1/np.sinh(2*k)*np.sqrt(pulsation(Bo, k) * k**2*Oh/2) +
            2*k**2*Oh * (np.cosh(4*k)+np.cosh(2*k)-1) / (np.cosh(4*k) -1)
            - pow(k**2*Oh,3./2.)/np.sqrt(2*pulsation(Bo, k))
            *(3-8*np.cosh(2*k)-14*np.cosh(4*k)+4*np.cosh(6*k))/(8*np.sinh(2*k)**3)) 

def om_analytic(Oh, Bo, k):
    try:
        root_denom = findroot(lambda s: denom (s, Oh, Bo, k), om_lub(Oh, Bo, k))
    except ValueError:
        root_denom = findroot(lambda s: denom (s, Oh, Bo, k), j*pulsation(Bo, k))
    return root_denom