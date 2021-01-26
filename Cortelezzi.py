# -*- coding: utf-8 -*-

##Import 

import numpy as np

from mpmath import mp
from mpmath import cosh, sinh, tanh, exp, sqrt

from gwr_inversion import gwr
M_value = 32

import matplotlib.pyplot as plt

## Function and expression declarations
#Declare the expressions of the kernel, A and eta
def ker_sy (s, Oh, Bo, k, lbda):
    return 2*Oh*k*(k-lbda*tanh(k))-Oh*(4*lbda*k*sinh(k)*(k*exp(-lbda*k)*(k*cosh(k)+lbda*sinh(k))-(k**2+lbda**2))+(k**2+lbda**2)**2*sinh(lbda))/(2*k*cosh(k)*(k*cosh(k)*sinh(lbda)-lbda*sinh(k)*cosh(lbda)))
def eta_sy (s, Oh, k, omega2, Kern):
    return 1/s*(1-omega2/(s**2+4*Oh*k**2*s+omega2+2*Oh*k**2*s*Kern))
def ALaplace_sy (s, Oh, k, lbda, omega2, z, Kern):
    return 1/sinh(lbda)*(-(-2*k*lbda*sinh(k)+(k**2+lbda**2)*sinh(lbda))*sinh(lbda*z)/(k*(-lbda*cosh(lbda)*sinh(k)+k*cosh(k)*sinh(lbda)))+2*sinh(lbda*(1+z)))*(-omega2/(s**2+4*Oh*k**2*s+omega2+2*Oh*k**2*s*Kern))

#Reduce the expressions as functions of s and of the parameters Oh, Bo and k
def freeSurfaceLaplace(s, Oh, Bo, k):
    lbda = sqrt(k**2 + s/Oh)
    omega2 = (Bo+k**2)*k*tanh(k)
    ker = ker_sy (s, Oh, Bo, k, lbda)
    return eta_sy(s, Oh, k, omega2, ker)

def omegaLaplace(s, z, Oh, Bo, k):
    lbda = sqrt(k**2 + s/Oh)
    omega2 = (Bo+k**2)*k*tanh(k)
    ker = ker_sy (s, Oh, Bo, k, lbda)
    return ALaplace_sy (s, Oh, k, lbda, omega2, z, ker)

#Inverse the Laplace transfrom and return functions of t and of the parameters Oh, Bo and k
def freeSurface(t, Oh, Bo, k):
    return mp.invertlaplace(lambda s: freeSurfaceLaplace(s, Oh, Bo, k), t)

def omega(t, z, Oh, Bo, k):
    return mp.invertlaplace(lambda s: omegaLaplace(s, z, Oh, Bo, k), t)

##Main function
##################################################################"
#def plotComparison(Oh Bo, k, nbOfrelax):
nbOfrelax = 1
Oh = mp.mpmathify('10')
Bo = mp.fraction(1,1000)
k = mp.fraction(1,10) 

#Resolution functions with numerical Laplace Transform GWR. Oh_, Bo_ and k_ must have been defined before
def eta(t):
    return freeSurface(t, Oh, Bo, k)

def vorti(t, z_):    
    return omega(t, z_, Oh, Bo, k)


#determine the maximal value of vorticity
t_all = np.linspace(1, 99, 100)*(mp.fraction(3,100)/Oh)
z_top = mp.mpmathify('-99/100')
vortmax = [vorti(t, z_top) for t in t_all]
maxv = max(np.array(vortmax))

#Spatial and time resolutions
tau_relax = 3*(Oh/(k**2*Bo+k**4))
timesOfInterest = np.array([1, 2, 5, 10])*mp.fraction(nbOfrelax*tau_relax, 10)
z_all = np.linspace(-99, -1, 99)*mp.fraction(1,100)
t_all_eta = np.linspace(1, 100)* mp.fraction(nbOfrelax*tau_relax, 100)

#solve the equation on omega with Cortelezzi model and lubrication
sampled_omega = np.array([[[2*i+vorti(timesOfInterest[i], z)/maxv, z] for i in range(4)] for z in z_all])
sampled_omega_lub = np.array([[[2*i+k**3/(Oh+k*Bo/Oh)*(-z)*np.exp(np.float(-timesOfInterest[i]/tau_relax))/maxv, z] for i in range(4)] for z in z_all])
sampled_t = t_all_eta/tau_relax
sampled_eta = [eta(t) for t in t_all_eta]
sampled_eta_lub = [np.exp(np.float(-t/tau_relax)) for t in t_all_eta]