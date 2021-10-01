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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker
import vapeplot

USETEX = True
clrs = vapeplot.palette('vaporwave')
clrlub=clrs[2]
clrpole=clrs[6]

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
def freeSurface(t, Oh, Bo, k):
    return gwr(lambda s: freeSurfaceLaplace(s, Oh, Bo, k), t, M_value)

## Main functions
def solveEta(Ohnumb, Bonumb, knumb, nbOfrelax):
    Oh = mp.mpmathify(Ohnumb)
    Bo = mp.mpmathify(Bonumb)
    k = mp.mpmathify(knumb) 
       
    #Time resolutions
    tau_relax = float(3*(Oh/(k**2*Bo+k**4)))
    t_all = np.linspace(0.0001, 1., 300) * nbOfrelax * abs(tau_relax)
    
    #solve the equation on eta with Cortelezzi model and lubrication
    sampled_t = abs(t_all/tau_relax)
    sampled_eta = [float(freeSurface(t, Oh, Bo, k)) for t in t_all]
    sampled_eta_lub = [np.exp(-t/tau_relax) for t in t_all]
    
    return sampled_t, sampled_eta, sampled_eta_lub


## Parameters figures
if USETEX:
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amssymb} \usepackage[squaren,Gray]{SIunits} \usepackage{nicefrac}'
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'cm'
else:
    plt.rcParams['font.sans-serif'] = "Roboto"
    plt.rcParams['font.weight'] = "light"
    plt.rcParams['font.family'] = "sans-serif" # always use sans-serif fonts

#font size
plt.rc('font', size=10)  # general font size
plt.rc('axes', labelsize=11, titlesize=10, linewidth=2.)
plt.rc('lines', markersize=8, markeredgewidth=0., linewidth=0.4)
#plt.rc('legend', frameon=False, fancybox=False, numpoints=1, markerscale=1, 
#       fontsize=10, handlelength=0.6, handletextpad=0.6, labelspacing=0.3)
plt.rc('xtick',  labelsize=10, direction='in', bottom='true', top='true')
plt.rc('ytick',  labelsize=10, direction='in', left='true', right='true')
plt.rc('savefig', bbox='tight', transparent=True, dpi=300) 

#%%
##Figure 1
def plotHeight(Ohnumb, Bonumb, knumb, ax, labelinset):
    sampled_t, sampled_eta, sampled_eta_lub = solveEta(Ohnumb, Bonumb, knumb, 4.*max(1.,knumb**2/Ohnumb))

    om_lub = float((knumb**2*Bonumb+knumb**4)/(3*Ohnumb))
    om_0 = np.sqrt(abs((Bonumb+knumb**2)*knumb*np.tanh(knumb)))
    try:
        root_denom = findroot(lambda s: denom (s, Ohnumb, Bonumb, knumb), om_lub)
    except ValueError:
        root_denom = findroot(lambda s: denom (s, Ohnumb, Bonumb, knumb), j*om_0)
    
   # ax.set_title("Oh = " + str(Ohnumb) + ", k = " + str(knumb))
    ax.plot(sampled_t, np.abs(decaying_sinusoid(sampled_t, 
                                                float(-mp.re(root_denom/om_lub)), 
                                                float(mp.im(root_denom/om_lub)))), 
            color=clrpole, ls=":", dash_capstyle="round", 
            linewidth=2, label = 'Analytical resolution')
    ax.plot(sampled_t,sampled_eta_lub, color=clrlub, ls=(0,(0.01,2)), linewidth=2, dash_capstyle="round", label = 'Lubrication theory')
    ax.plot(sampled_t,np.abs(sampled_eta), color=clrs[8], dash_capstyle="round", linewidth=2)
    if USETEX:
        text = r'Time (in $\tau_\mathrm{relax}$ units)'
    else:
        text = r'Time (in $\mathregular{\tau_{relax}}$ units)'        
    ax.set_xlabel(text, family = "Roboto", weight="ultralight")
    textOh = "Oh = " + str(Ohnumb) 
    textk  = "k = " + str(knumb)
    axinset = inset_axes(ax, width="40%", height="40%", borderpad=0.5)
    axinset.linewidth=0.5

#    axinset.tick_params(axis=u'both', which=u'both',length=0)
    axinset.tick_params(axis=u'both', which=u'both',width=0.4)
    props = dict(facecolor='white', alpha=0.8, edgecolor="None")
    plt.setp(axinset.get_xticklabels(), bbox=props,
             family = "Roboto", size=6, weight="light")
    plt.setp(axinset.get_yticklabels(), bbox=props,
             family = "Roboto", size=6, weight="light")
    ax.text(0.09, 0.92, textOh, size="x-large", 
            ha="left", va="center", family = "Source Sans Pro", weight="ultralight",
            transform=ax.transAxes)
    ax.text(0.15, 0.85, textk, size="x-large", 
            ha="left", va="center", family = "Source Sans Pro", weight="ultralight",
            transform=ax.transAxes)
    error_lub=np.subtract(sampled_eta,sampled_eta_lub)
    error_pole=np.subtract(sampled_eta,decaying_sinusoid(sampled_t, float(-mp.re(root_denom/om_lub)), float(mp.im(root_denom/om_lub))))
    axinset.semilogy(sampled_t, np.abs(error_lub), color=clrlub, linewidth=1.5)
    axinset.semilogy(sampled_t, np.abs(error_pole), color=clrpole, linewidth=1.5)
        ## set y ticks
    y_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = 5)
    axinset.yaxis.set_major_locator(y_major)
    y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    axinset.yaxis.set_minor_locator(y_minor)
    axinset.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axinset.set_xticks(sampled_t[-1]*np.arange(5)/4)

 #   axinset.set_xlim(0, 4), axinset.set_xticks(0.5*np.arange(1,8))
  #  axinset.set_ylim(0, 8), axinset.set_yticks(np.arange(1,8))
    axinset.grid(True, axis='both', which='both', linewidth=0.125, alpha=0.5) # which='major',
    axinset.set_axisbelow(True)
    for axis in ['top','bottom','left','right']:
        axinset.spines[axis].set_linewidth(0.5)
    axinset.patch.set_alpha(0.5)
    if USETEX:
        textxinset = r'\textbf{Time}'
        textyinset = r'\textbf{Abs. error}'
    else:
        textxinset = r'Time'
        textyinset = r'Abs. error'
    axinset.set_xlabel(text, family = "Roboto", weight="ultralight", fontsize=7)
    axinset.set_ylabel(text, family = "Roboto", weight="ultralight", fontsize=7)
    for i in axinset.xaxis.get_ticklines():
        # reach in and poke the guts 
        # USE AT YOUR OWN RISK
        i._marker._capstyle = 'round' 
        # this is not officially supported
    for i in axinset.yaxis.get_ticklines():
        # reach in and poke the guts 
        # USE AT YOUR OWN RISK
        i._marker._capstyle = 'round' 
        # this is not officially supported
    if labelinset:
        if USETEX:
            axinset.text(0.33, 0.25, r'\textbf{\textsf{discrete pole}}', size="x-small", 
                         ha="left", va="center", family = "sans-serif", weight="bold",
                         transform=axinset.transAxes, color=clrpole)
            axinset.text(0.15, 0.85, r'\textbf{\textsf{lubrication}}', size="x-small", 
                         ha="left", va="center", family = "sans-serif", weight="bold",
                         transform=axinset.transAxes, color=clrlub)
        else:
            axinset.text(0.33, 0.25, "discrete pole", size="x-small", 
                         ha="left", va="center", family = "Roboto Condensed", weight="bold",
                         transform=axinset.transAxes, color=clrpole)
            axinset.text(0.15, 0.85, "lubrication", size="x-small", 
                         ha="left", va="center", family = "Roboto Condensed", weight="bold",
                         transform=axinset.transAxes, color=clrlub)
  

fig, ax = plt.subplots(ncols = 2, figsize=(6+6/8, 3+3/8), sharey=True)
for i, axis in enumerate(ax):
    axis.set_box_aspect(1)
plotHeight(10, 0.001, 0.1, ax[0], True)
plotHeight(0.01, 0.001, 0.5, ax[1], False)

#ax.set_aspect('square')

#lines, labels = ax[-1].get_legend_handles_labels()
#fig.legend(lines, labels, loc = 'lower center', borderaxespad=0.1, ncol=3)
ax[0].set_ylabel('Relative amplitude', family = "Roboto", weight="ultralight")
fig.savefig("basic-plot.pdf")
plt.tight_layout(pad=2.)

