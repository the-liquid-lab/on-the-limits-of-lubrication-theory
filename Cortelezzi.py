# -*- coding: utf-8 -*-
"""
@author: Clément & Arnaud
"""
#Must install the fonts Roboto and Roboto condensed
#Must install the packages gwr_inversion and mpmath

#%% 
# Import
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import ticker
# Ensure working directory
import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

from Functions import freeSurface, om_lub, pulsation
from Functions import om_normal_mode_viscous
from Compute_datas import compute_datas_fig1, compute_datas_fig2, compute_datas_fig3
from Plots import fig_init, plot_fig2, plot_fig3
#Colors
import vapeplot
clrs = vapeplot.palette('vaporwave')
clrlub=clrs[2]
clrpole=clrs[6]

fig_init()

#%% Figure 1
#Comparison between lubrication, analytical and numerical results for 2 different situations : oscillations and relaxation
def plotHeight(Oh, Bo, k, ax, labelinset):

    t, eta, eta_lub, eta_ana = compute_datas_fig1(Oh, Bo, k)
    
    ax.plot(t, np.abs(eta_ana), color=clrpole, ls=":", dash_capstyle="round", 
            linewidth=2, label = 'Analytical resolution')
    ax.plot(t, eta_lub, color=clrlub, ls=(0,(0.01,2)), 
            linewidth=2, dash_capstyle="round", label = 'Lubrication theory')
    ax.plot(t, np.abs(eta), color=clrs[8], dash_capstyle="round", linewidth=2)
    text = r'Time (in $\mathregular{\tau_{relax}}$ units)'        
    ax.set_xlabel(text, family = "Roboto", weight="ultralight")
    textOh = "Oh = " + str(Oh) 
    textk  = "k = " + str(k)
    axinset = inset_axes(ax, width="40%", height="40%", borderpad=0.5)
    axinset.linewidth=0.5

    axinset.tick_params(axis=u'both', which=u'both',width=0.2)
    props = dict(facecolor='white', alpha=0.8, edgecolor="None")
    plt.setp(axinset.get_xticklabels(), bbox=props,
             family = "Roboto", size=6, weight="light")
    plt.setp(axinset.get_yticklabels(), bbox=props,
             family = "Roboto", size=6, weight="light")
    ax.text(0.09, 0.92, textOh, size="large", 
            ha="left", va="center", family = "Source Sans Pro", weight="ultralight",
            transform=ax.transAxes)
    ax.text(0.13, 0.85, textk, size="large", 
            ha="left", va="center", family = "Source Sans Pro", weight="ultralight",
            transform=ax.transAxes)
    error_lub=np.subtract(eta,eta_lub)
    error_pole=np.subtract(eta,eta_ana)
    axinset.semilogy(t, np.abs(error_lub), color=clrlub, linewidth=1.5)
    axinset.semilogy(t, np.abs(error_pole), color=clrpole, linewidth=1.5)
        ## set y ticks
    y_major = ticker.LogLocator(base = 10.0, numticks = 5)
    axinset.yaxis.set_major_locator(y_major)
    y_minor = ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    axinset.yaxis.set_minor_locator(y_minor)
    axinset.yaxis.set_minor_formatter(ticker.NullFormatter())
    axinset.set_xticks(t[-1]*np.arange(5)/4)

 #   axinset.set_xlim(0, 4), axinset.set_xticks(0.5*np.arange(1,8))
  #  axinset.set_ylim(0, 8), axinset.set_yticks(np.arange(1,8))
    axinset.grid(True, axis='both', which='both', linewidth=0.125, alpha=0.5) # which='major',
    axinset.set_axisbelow(True)
    for axis in ['top','bottom','left','right']:
        axinset.spines[axis].set_linewidth(0.5)
    axinset.patch.set_alpha(0.5)
#    textxinset = r'Time'
#    textyinset = r'Abs. error'
    axinset.set_xlabel(text, family = "Roboto", weight="ultralight", fontsize=6)
    axinset.set_ylabel(text, family = "Roboto", weight="ultralight", fontsize=6)
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
        axinset.text(0.33, 0.25, "discrete pole", size="x-small", 
                     ha="left", va="center", family = "Roboto Condensed", weight="bold",
                     transform=axinset.transAxes, color=clrpole)
        axinset.text(0.15, 0.85, "lubrication", size="x-small", 
                     ha="left", va="center", family = "Roboto Condensed", weight="bold",
                     transform=axinset.transAxes, color=clrlub)

            
    # ax.plot(t[::8],np.abs(eta), '.b', ms = 6., label = r'Numerical resolution')
    # ax.plot(t, np.abs(decaying_sinusoid(t, float(-mp.re(om_ana)), float(mp.im(om_ana)))), 'red', label = 'Analytical resolution')
    # ax.plot(t,eta_lub, 'green', label = 'Lubrication theory')
    # ax.set_xlabel('Time (in $\tau_{relax}$)')
    
fig, ax = plt.subplots(ncols = 2, figsize=(6+6/8, 3+3/8), sharey=True)
for i, axis in enumerate(ax):
    axis.set_box_aspect(1)
plotHeight(10, 0.001, 0.1, ax[0], True)
plotHeight(0.01, 0.001, 0.5, ax[1], False)

#ax.set_aspect('square')

#lines, labels = ax[-1].get_legend_handles_labels()
#fig.legend(lines, labels, loc = 'lower center', borderaxespad=0.1, ncol=3)
ax[0].set_ylabel('Relative amplitude', family = "Roboto", weight="ultralight")
plt.tight_layout(pad=1.)
fig.savefig("figure1.pdf")

#%% Figure 2

## Parameters
Bo = 0.001
k = 0.5
Oh_list = np.logspace(-3.5, 0.5, 600)

## Import datas
om_ana, om_0 = compute_datas_fig2(Oh_list, Bo, k)

## Plot datas
plot_fig2(Oh_list, Bo, k, om_ana, om_0)

#%% Figure 3

## Parameters
Oh_list = np.logspace(-3, 1, 50)
k_list = np.logspace(-2, 2, 50)
Bo = 1

## Import datas
#The data can be easily recompute but it takes about 1h.
#For time efficiency, numerical values are by default taken in the txt file. 
ComputeAllFig3 = False
err_lub, err_visc, err_in, splitline = compute_datas_fig3 (Oh_list, k_list, Bo, 'fig3_om_num.npy', ComputeAllFig3)

## Plot datas
plot_fig3(Oh_list, k_list, err_lub, err_visc, err_in, splitline)


#%% Figure 4 
#Rayleigh-Taylor
from scipy import stats
compute = False 

def growth_rate(Oh, Bo, k):
    t_all = np.linspace(0.001, 25., 50)/k
    sampled_eta = freeSurface(t_all, Oh, Bo, k)
    
    reg = stats.linregress(t_all[20:], np.log(sampled_eta[20:]))
    if (reg[2]<0.999):
        print(Oh, k, reg[2])
        plt.figure()
        plt.xlabel(r'Time (in $\tau_{relax}$ units)')
        plt.ylabel("Relative wave amplitude") 
        plt.semilogy(t_all*abs(om_lub(Oh, Bo, k)), sampled_eta, 'black', label = r'Cortelezzi \& Prosperetti')
        plt.semilogy(t_all*abs(om_lub(Oh, Bo, k)), np.exp(reg[1] + t_all*reg[0]), 'gray', label = 'Regression')
    return reg[0]

def om_normal_mode_inertial_Bo_neg(Oh, Bo, k):
    return pulsation(Bo, k) - (1/np.sinh(2*k)*np.sqrt(pulsation(Bo, k) * k**2*Oh/2)
            + 2*k**2*Oh * (np.cosh(4*k)+np.cosh(2*k)-1) / (np.cosh(4*k) -1)
            + pow(k**2*Oh,3./2.)/np.sqrt(2*pulsation(Bo, k))
            *(3-8*np.cosh(2*k)-14*np.cosh(4*k)+4*np.cosh(6*k))/(8*np.sinh(2*k)**3))

Bo = -0.5
Oh_list = [0.01, 1.]
k_list = np.linspace(0.005, 0.999, 100) * np.sqrt(-Bo)
k_list2 = np.linspace(0.005, 1., 100) * np.sqrt(-Bo)

if compute:
    om_gwr_Oh = []
    for Oh in Oh_list:
        om_gwr_Oh.append([growth_rate(Oh, Bo, k) for k in k_list])
    np.save("RayleighTaylor",om_gwr_Oh)
else:
    om_gwr_Oh = np.load("RayleighTaylor.npy")
   
om_potential = [pulsation(Bo, k) for k in k_list]


fig, ax = plt.subplots(1,2, figsize=(16,8))
    
Oh = Oh_list[0]
om_norm = [np.abs(om_normal_mode_inertial_Bo_neg(Oh, Bo, k)) for k in k_list2]
ax[1].plot(k_list, om_potential, lw=1.0, alpha = 0.4, color = 'black', label = r'Potential')
ax[1].plot(k_list2, om_norm, '-', lw=1.0, alpha = 0.4, color = 'red', label = 'Normal mode')

Oh = Oh_list[1]
ax[0].set_ylabel(r'$\omega$')
om_lub_list = [np.abs(om_lub(Oh, Bo, k)) for k in k_list2]
ax[0].plot(k_list2, om_lub_list, '-', lw=1.0, alpha = 0.4, color = 'blue', label = 'Lubrication')
om_norm = [np.abs(om_normal_mode_viscous(Oh, Bo, k)) for k in k_list2]
ax[0].plot(k_list2, om_norm, '-', lw=1.0, alpha = 0.4, color = 'red', label = 'Normal mode')

for Oh, axx, om_gwr in zip(Oh_list, [ax[1],ax[0]], om_gwr_Oh):
    axx.set_xlabel(r'k')
    axx.set_title('Oh = ' + str(Oh))
    axx.plot(k_list, np.abs(om_gwr), '--', lw=1.0, color = 'orange', alpha = 0.8, label = r'Cortelezzi resolution')
    axx.legend()

plt.tight_layout(pad=2.)

# Oh = 0.01
# k = 0.5
# t_all = np.linspace(0.001, 25., 50)/k
# sampled_eta = freeSurface(t_all, Oh, Bo, k)
# reg = stats.linregress(t_all[20:], np.log(sampled_eta[20:]))
# plt.figure()
# plt.xlabel(r'Time (in $\tau_{relax}$ units)')
# plt.ylabel("Relative wave amplitude") 
# plt.plot(t_all*abs(om_lub(Oh, Bo, k)), sampled_eta, 'black', label = r'Cortelezzi \& Prosperetti')
# plt.plot(t_all*abs(om_lub(Oh, Bo, k)), np.exp(reg[1] + t_all*reg[0]), 'gray', label = 'Regression')