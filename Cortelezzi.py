# -*- coding: utf-8 -*-
"""
@author: Clément & Arnaud
"""

## Import
import numpy as np

from mpmath import mp, findroot, j
from mpmath import cosh, sinh, tanh, exp, sqrt

from gwr_inversion import gwr #The package must be installed through "conda install gwr_inversion"
M_value = 32

import matplotlib.pyplot as plt
## Function and expression declarations

def decaying_sinusoid(t, om_dec, om_osc):
    return np.exp(- om_dec * t)*np.cos(om_osc * t)

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

#Inverse the Laplace transfrom and return functions of t and of the parameters Oh, Bo and k
def freeSurface(t, Ohnumb, Bonumb, knumb):
    Oh = mp.mpmathify(Ohnumb)
    Bo = mp.mpmathify(Bonumb)
    k = mp.mpmathify(knumb) 
    return gwr(lambda s: freeSurfaceLaplace(s, Oh, Bo, k), t, M_value)

def om_lub(Oh, Bo, k):
    return float((k**2*Bo+k**4)/(3*Oh))

def pulsation(Bo, k):
    return np.sqrt(np.abs(Bo + k**2)*k*np.tanh(k))
    
## Parameters figures
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[squaren,Gray]{SIunits} \usepackage{nicefrac}'

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'cm'
#font size
plt.rc('font', size=10)  # general font size
plt.rc('axes', labelsize=10, titlesize=10)
plt.rc('lines', markersize=8, markeredgewidth=0., linewidth=0.4)
plt.rc('legend', frameon=False, fancybox=False, numpoints=1, markerscale=1, 
       fontsize=10, handlelength=0.6, handletextpad=0.6, labelspacing=0.3)
plt.rc('xtick',  labelsize=8, direction='in', bottom='true', top='true')
plt.rc('ytick',  labelsize=8, direction='in', left='true', right='true')
plt.rc('savefig', bbox='tight', transparent=True, dpi=300) 

#%%
##Figure 1
def plotHeight(Oh, Bo, k, ax):
    
    nbOfrelax = 4.*max(1.,k**2/Oh)
    om_relax = om_lub(Oh, Bo, k)
    om_0 = pulsation(Bo, k)
    
    t_all = np.linspace(0.0001, 1., 300) * nbOfrelax * abs(1./om_relax)
    sampled_t = abs(t_all*om_relax)
    sampled_eta = [float(freeSurface(t, Oh, Bo, k)) for t in t_all]
    sampled_eta_lub = [np.exp(-t*om_relax) for t in t_all]
    
    try:
        root_denom = findroot(lambda s: denom (s, Oh, Bo, k), om_relax)
    except ValueError:
        root_denom = findroot(lambda s: denom (s, Oh, Bo, k), j*om_0)
    
    ax.set_title("Oh = " + str(Oh) + ", k = " + str(k))
    ax.plot(sampled_t[::8],np.abs(sampled_eta[::8]), '.b', ms = 6., label = r'Numerical resolution')
    ax.plot(sampled_t, np.abs(decaying_sinusoid(sampled_t, float(-mp.re(root_denom/om_relax)), float(mp.im(root_denom/om_relax)))), 'red', label = 'Analytical resolution')
    ax.plot(sampled_t,sampled_eta_lub, 'green', label = 'Lubrication theory')
    ax.set_xlabel('Time (in $\tau_{relax}$)')
    
fig, ax = plt.subplots(ncols = 2, figsize=(8, 4))
plotHeight(10., 0.001, 0.1, ax[0])
plotHeight(0.01, 0.001, 0.5, ax[1])

lines, labels = ax[-1].get_legend_handles_labels()
fig.legend(lines, labels, loc = 'lower center', borderaxespad=0.1, ncol=3)
ax[0].set_ylabel('Relative amplitude')
plt.tight_layout(pad=2.)

#%%
##Figure 2


#%%
##Figure 3

#%%
##Figure 4
from scipy import stats

def growth_rate(Oh, Bo, k):
    t_all = np.linspace(0.001, 25., 50)/k
    sampled_eta = [float(freeSurface(t, Oh, Bo, k)) for t in t_all]
    
    reg = stats.linregress(t_all[20:], np.log(sampled_eta[20:]))
    if (reg[2]<0.999):
        print(Oh, k, reg[2])
        plt.figure()
        plt.xlabel(r'Time (in $\tau_{relax}$ units)')
        plt.ylabel("Relative wave amplitude") 
        plt.semilogy(t_all*abs(om_lub(Oh, Bo, k)), sampled_eta, 'black', label = r'Cortelezzi \& Prosperetti')
        plt.semilogy(t_all*abs(om_lub(Oh, Bo, k)), np.exp(reg[1] + t_all*reg[0]), 'gray', label = 'Regression')
    return reg[0]

Bo = -0.5
Oh_list = [0.01, 1.]
k_list = np.linspace(0.005, 0.999, 100) * np.sqrt(-Bo)
k_list2 = np.linspace(0.005, 1., 100) * np.sqrt(-Bo)

om_gwr_Oh = []
om_lub_Oh = []
for Oh in Oh_list:
    om_gwr_Oh.append([growth_rate(Oh, Bo, k) for k in k_list])
    om_lub_Oh.append([np.abs(om_lub(Oh, Bo, k)) for k in k_list2])
om_potential = [pulsation(Bo, k) for k in k_list]
Colors = ['orange', 'green', 'black']

plt.figure()
plt.xlabel(r'k')
plt.ylabel(r'$\omega$')
plt.loglog(k_list, om_potential, lw=1.0, alpha = 0.4, color = Colors[-1], label = r'Potential')
for Oh, om_gwr, om_lub, c in zip(Oh_list, om_gwr_Oh, om_lub_Oh, Colors):
    plt.plot(k_list, np.abs(om_gwr), '--', lw=1.0, color = c, alpha = 0.8, label = r'Cortelezzi resolution, Oh = ' + str(Oh))
    plt.plot(k_list2, om_lub, '-', lw=1.0, alpha = 0.4, color = c, label = 'Lubrication, Oh = ' + str(Oh))
plt.legend()
plt.tight_layout(pad=0.)





















#%%


def ALaplace_sy (s, Oh, k, lbda, omega2, z, Kern):
    return 1/sinh(lbda)*(-(-2*k*lbda*sinh(k)+(k**2+lbda**2)*sinh(lbda))
             *sinh(lbda*z)/(k*(-lbda*cosh(lbda)*sinh(k)+k*cosh(k)*sinh(lbda)))
             +2*sinh(lbda*(1+z)))*(-omega2/(s**2+4*Oh*k**2*s+omega2+2*Oh*k**2*s*Kern))

def vortiLaplace(s, z, Oh, Bo, k):
    lbda = sqrt(k**2 + s/Oh)
    omega2 = (Bo+k**2)*k*tanh(k)
    ker = ker_sy (s, Oh, Bo, k, lbda)
    return k*ALaplace_sy (s, Oh, k, lbda, omega2, z, ker)


def vorti(t, z, Oh, Bo, k):
    return gwr(lambda s: vortiLaplace(s, z, Oh, Bo, k), t, M_value)

def solveOmega(Ohnumb, Bonumb, knumb, nbOfrelax, toI, z_all):
    Oh = mp.mpmathify(Ohnumb)
    Bo = mp.mpmathify(Bonumb)
    k = mp.mpmathify(knumb) 
    
    def omega(t, z_):    
        return float(vorti(t, z_, Oh, Bo, k))
       
    #Spatial and time resolutions
    tau_relax = abs(float(3*(Oh/(k**2*Bo+k**4))))
    t_all = np.linspace(0.0001, 1., 300) * nbOfrelax * tau_relax
    timesOfInterest = np.array(toI)*tau_relax

    #determine the maximal value of vorticity
    vortmax = [abs(omega(t, mp.mpmathify(z_all[0]))) for t in t_all]
    maxv = max(np.array(vortmax))
    
    #solve the equation on omega with Cortelezzi model and lubrication
    sampled_omega = np.array([[2*i+omega(timesOfInterest[i], mp.mpmathify(z))/maxv
        for i in range(len(timesOfInterest))] for z in z_all])
    sampled_omega_lub = np.array([[2*i+float(k**3/Oh+k*Bo/Oh)*(-z)
        *np.exp(-timesOfInterest[i]/tau_relax)/maxv
        for i in range(len(timesOfInterest))] for z in z_all])
    
    return sampled_omega, sampled_omega_lub

## Relaxation time

from scipy.optimize import curve_fit
def decaying_sinusoid(t, om, k):
    return np.exp(-om * t)* np.cos(k * t)

from scipy.signal import hilbert

def plotHeight(Ohnumb, Bonumb, knumb, ax):
    sampled_t, sampled_eta, sampled_eta_lub = solveEta(Ohnumb, Bonumb, knumb, 5.*max(1.,knumb**2/Ohnumb))
    try:
        popt, pcov = curve_fit(decaying_sinusoid, sampled_t,sampled_eta, p0=(min(Ohnumb/knumb**2, 0.1), 0.), bounds=(0,[1,np.inf]))
        om_cort = popt[0]
    except RuntimeError :
        om_cort = -1
        
    dacaying = np.exp(-om_cort * sampled_t)
    enveloppe = np.abs(hilbert(sampled_eta))
    
    om_lub = float((knumb**2*Bonumb+knumb**4)/(3*Ohnumb))
    om_0 = np.sqrt(abs((Bonumb+knumb**2)*knumb*np.tanh(knumb)))
    om_diff = knumb**2/Ohnumb
    
    ax.set_title(om_cort)
    ax.plot(sampled_t,np.abs(sampled_eta), 'black', label = r'Cortelezzi \& Prosperetti')
    ax.plot(sampled_t,sampled_eta_lub, 'grey', label = 'Lubricaion theory')
    ax.plot(sampled_t,enveloppe, 'green', label = 'Enveloppe')
    ax.plot(sampled_t,dacaying, 'red', label = 'Decaying')
    ax.set_ylim([0,1])
    
Oh_list = np.logspace(-2, 1, 3)
k_list = np.logspace(0, -2, 3)

fig, ax = plt.subplots(ncols = 3, nrows = 3, figsize=(8, 8))
j = 0
for Oh in Oh_list:
    i = 0
    for k in k_list:
        plotHeight(Oh, 0.001, k, ax[i,j])
        i+=1
    j+=1


# In[]:


Bonumb = 0.001
def om_cort(knumb, Ohnumb):
    sampled_t, sampled_eta, sampled_eta_lub = solveEta(Ohnumb, Bonumb, knumb, 5.*max(1.,knumb**2/Ohnumb))
    try:
        popt, pcov = curve_fit(decaying_sinusoid, sampled_t,sampled_eta, p0=(min(Ohnumb/knumb**2, 0.1), 0.), bounds=(0,[1,np.inf]))
        return popt[0]
    except RuntimeError :
        return -1

Oh_list = np.logspace(-3, 1, 20)
k_list = np.logspace(-2, 0, 20)
om = np.array([[om_cort(k, Oh) for Oh in Oh_list] for k in k_list])

print(om)


# In[]:


plt.contourf(Oh_list, k_list, np.log10(1-om), cmap='RdGy');
plt.xscale('log')
plt.yscale('log')
plt.plot(Oh_list,np.sqrt(Oh_list), linewidth = 2.)
cbar = plt.colorbar();
ctks = [-2, -1, 0]
cbar.set_ticks(ctks)
ctkls = ["$10^{%d}$"%(v) for v in ctks]
cbar.set_ticklabels(ctkls)


# # With Vorticity

# In[]:


def plotComparison(Ohnumb, Bonumb, knumb, nbOfrelax, toI = [0.05, 0.1]):
    
    sampled_t, sampled_eta, sampled_eta_lub =         solveEta(Ohnumb, Bonumb, knumb, nbOfrelax)
    
    z_all = np.linspace(-1, 0, 200)
    sampled_omega, sampled_omega_lub =         solveOmega(Ohnumb, Bonumb, knumb, nbOfrelax, toI, z_all)
    
    ### Figures ###
    
    plt.rc('figure', figsize=[6/2.54, 6/2.54])
    plt.figure()
    plt.xlabel(r'Time (in $\tau_{relax}$ units)')
    plt.ylabel("Relative wave amplitude")
    plt.xlim([0,nbOfrelax])
    plt.ylim([0,1])
    plt.rc('legend', loc='best')
    
    plt.plot(sampled_t,sampled_eta, 'black', label = r'Cortelezzi \& Prosperetti')
    plt.plot(sampled_t,sampled_eta_lub, 'gray', label = 'Lubricaion theory')
    plt.legend()
    plt.savefig('Oh_' + str(Ohnumb) + '_k_' + str(knumb) + '_eta.png')
    plt.tight_layout(pad=0.)

    plt.figure(figsize=[6./2.54, 6./2.54])
    plt.xlabel('Scaled vorticity')
    plt.ylabel('depth z')
    plt.xlim([-0.5, 2*len(sampled_omega[0])-0.5])
    plt.ylim([-1,0])
    plt.rc('legend', loc='best')
    ax = plt.gca()
    ax.axes.xaxis.set_ticks([])
    plt.tight_layout(pad=0.)
    
    for i in range(len(sampled_omega[0])):
        plt.plot(sampled_omega[:,i], z_all, 'r')
        plt.plot(2*i*np.ones(len(z_all)), z_all, 'r')
        plt.fill_betweenx(z_all, sampled_omega[:,i], 2*i, color = 'red', alpha=0.20)
        plt.plot(sampled_omega_lub[:,i],z_all, c = 'gray', linewidth=1, linestyle='dashed')
        plt.text(2.*i+0.2, -0.1, str(toI[i]) + r'$\tau_{relax}$', size=8)

    plt.savefig('Oh_' + str(Ohnumb) + '_k_' + str(knumb) + '_omg.png')

# In[]:
# ### Oh = 10, Bo = 0.001, k = 0.1
plotComparison(10, 0.001, 0.1, 1)                   # oilly film 

# ### Oh = 0.01, Bo = 0.001, k = 0.1
plotComparison(0.01, 0.001, 0.1, 1)                 # waterborne coating

# ### Oh = 0.005, Bo = 0.001, k = 0.1
plotComparison(0.005, 0.001, 0.1, 6, [0.3, 0.6])    # oscillations

# # Feuilletage
plotComparison(0.005, 0.001, 2., 2400, [1800, 2400])
plotComparison(0.005, 0.001, 1.5, 1800, [1200, 1800])
plotComparison(1e-4, 0.001, 2., 1./6.25e-6, [1./6.25e-6]) #tau = (Oh/k**4)


# ## Feuilletage simplifié

# In[]:


def Feuilletage(Ohnumb, Bonumb, knumb, nbOfrelax, toI = [0.05, 0.1]):

    plt.rc('figure', figsize=[6/2.54, 6/2.54])   
    z_all = np.linspace(-1, -0.999995, 1000)
    sampled_omega, sampled_omega_lub =         solveOmega(Ohnumb, Bonumb, knumb, nbOfrelax, toI, z_all)
    
    ### Figures ###
    plt.xlabel('Scaled vorticity')
    plt.ylabel('depth z')
    plt.xlim([-0.5, 1.5])
    plt.ylim([-1,-0.999995])

    plt.plot(sampled_omega[:,0], z_all, 'r')
    plt.fill_betweenx(z_all, sampled_omega[:,0], 0., color = 'red', alpha=0.20)
    
Oh = 1e-10
k = 100.
tau = (Oh/k**4)
Feuilletage(Oh, 0.001, k, 5.e-3/tau, [5.25e-3/tau])
 
sampled_t, sampled_eta, sampled_eta_lub = solveEta(0.01, -0.03, 0.05, 3.)
plt.rc('figure')
plt.rc('legend', loc='best')

plt.figure()
plt.xlabel(r'Time (in $\tau_{relax}$ units)')
plt.ylabel("Relative wave amplitude")
 #   plt.xlim([0,nbOfrelax])
 #   plt.ylim([1,10.])
plt.semilogy(sampled_t,sampled_eta, 'black', label = r'Cortelezzi \& Prosperetti')
plt.semilogy(sampled_t,sampled_eta_lub, 'gray', label = 'Lubricaion theory')
plt.legend()
plt.tight_layout(pad=0.)