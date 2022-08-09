# -*- coding: utf-8 -*-
"""
@author: Clément & Arnaud
"""

#%% Import
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker
import vapeplot
clrs = vapeplot.palette('vaporwave')
clrlub=clrs[2]
clrpole=clrs[6]

USETEX = True

from mpmath import mp, findroot, j 
from mpmath import cosh, sinh, tanh, exp, sqrt
from scipy.optimize import curve_fit
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

#%% Figure 1
#Comparison between lubrication, analytical and numerical results for 2 different situations : oscillations and relaxation
def om_analytic(Oh, Bo, k):
    try:
        root_denom = findroot(lambda s: denom (s, Oh, Bo, k), om_lub(Oh, Bo, k))
    except ValueError:
        root_denom = findroot(lambda s: denom (s, Oh, Bo, k), j*pulsation(Bo, k))
    return root_denom

def plotHeight(Oh, Bo, k, ax, labelinset):
    om_lub_relax = om_lub(Oh, Bo, k)
    om_0 = pulsation(Bo, k)
    om_ana = om_analytic(Oh, Bo, k)/om_lub_relax
    
    t_all = np.linspace(0.0001, 1., 300) * max(abs(5./om_lub_relax), 13./om_0)
    sampled_t = abs(t_all*om_lub_relax)
    sampled_eta = freeSurface(t_all, Oh, Bo, k)
    sampled_eta_lub = np.exp(-t_all*om_lub_relax)
    
    ax.plot(sampled_t, np.abs(decaying_sinusoid(sampled_t, 
                                                float(-mp.re(om_ana)), 
                                                float(mp.im(om_ana)))), 
            color=clrpole, ls=":", dash_capstyle="round", 
            linewidth=2, label = 'Analytical resolution')
    ax.plot(sampled_t,sampled_eta_lub, color=clrlub, ls=(0,(0.01,2)), linewidth=2, dash_capstyle="round", label = 'Lubrication theory')
    ax.plot(sampled_t,np.abs(sampled_eta), color=clrs[8], dash_capstyle="round", linewidth=2)
    if USETEX:
        text = r'Time (in $\tau_\mathrm{relax}$ units)'
    else:
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
    error_lub=np.subtract(sampled_eta,sampled_eta_lub)
    error_pole=np.subtract(sampled_eta,decaying_sinusoid(sampled_t, float(-mp.re(om_ana)), float(mp.im(om_ana))))
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

            
    # ax.plot(sampled_t[::8],np.abs(sampled_eta), '.b', ms = 6., label = r'Numerical resolution')
    # ax.plot(sampled_t, np.abs(decaying_sinusoid(sampled_t, float(-mp.re(om_ana)), float(mp.im(om_ana)))), 'red', label = 'Analytical resolution')
    # ax.plot(sampled_t,sampled_eta_lub, 'green', label = 'Lubrication theory')
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
Bo = 0.001
k = 0.5
Oh_list = np.logspace(-3.5, 0.5, 600)
om_ana = []
root_denom = j

for Oh in Oh_list:
    root_denom = findroot(lambda s: denom (s, Oh, Bo, k), root_denom)
    om_ana.append([float(mp.re(root_denom)), float(mp.im(root_denom))])
om_ana = np.array(om_ana)
om_0 = pulsation(Bo,k)

split = int(np.sum(om_ana[:,1]<0.018)/len(om_ana)*256)
# sample the colormaps that you want to use. Use 128 from each so we get 256
# colors in total
colors1 = plt.cm.Blues_r(np.linspace(0., 0.7, 256-split))
colors2 = plt.cm.Reds(np.linspace(0.3, 1, split))

# combine them and build a new colormap
colors = np.vstack((colors1, colors2))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

plt.figure(figsize=(5, 4))
p = [plt.scatter(0, 1, label = 'Natural pulsations', marker = 'P', s = 80, c = 'black'),
     plt.scatter(-0.93, 0, label = 'Split point', marker = '*', s = 60, c = 'black')]
plt.arrow(-0.91, 0.05, 0.05, 0.2, head_width = 0.02, color = 'black')
plt.arrow(-0.91, -0.05, 0.05, -0.2, head_width = 0.02, color = 'black')
plt.scatter(0, -1, marker = 'P', s = 80, c = 'black')
plt.scatter(om_ana[:,0]/om_0, om_ana[:,1]/om_0, s = 20, c = Oh_list, cmap=mymap, norm=matplotlib.colors.LogNorm())
plt.scatter(om_ana[:,0]/om_0, -om_ana[:,1]/om_0, s = 20, c = Oh_list, cmap=mymap, norm=matplotlib.colors.LogNorm())
plt.xlabel('$\omega_{relax}/\omega_0 = \Re(s/\omega_0)$', family = "Roboto", weight="ultralight")      
plt.ylabel('$\omega_{osc}/\omega_0 = \Im(s/\omega_0)$', family = "Roboto", weight="ultralight")
cbar = plt.colorbar(label = 'Oh')
plt.legend()
plt.tight_layout(pad=1.)
plt.savefig("figure2.pdf")

#%% Figure 3
# Relative error of different models compare to the numerical results.

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

def err_norm(relax, puls, om_num):
    relax_num = om_num[0] # 0 for decaying
    puls_num = om_num[1] # 1 for oscillation
    return np.sqrt((np.square(relax-relax_num) + 
                        np.square(puls-puls_num))/
                       (np.square(relax_num) + np.square(puls_num)))

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
    popt = curve_fit(better_sinusoid, t_all, sampled_eta, p0=guess_value, bounds=([0.,0.,1.,0.],[np.inf, 5.*om_0, 2., 2.*np.pi]), sigma=(1.+10.*np.exp(-om_relax*t_all)))[0]
    return popt, t_all, sampled_eta


#Compare the different models for a range of Oh and k.
def plotErrorOm (Oh_list, k_list, Bo, file_name, compute = False):
    #The data can be easily recompute but it takes about 1h.
    #For time efficiency, numerical values are by default taken in the txt file. 
    if compute:
        om_num = []
        value_k_Oh = (om_normal_mode_inertial(Oh_list[0], Bo, k_list[0]),
                       puls_normal_mode_inertial(Oh_list[0], Bo, k_list[0]),
                       1., 0.)
        for k in k_list:
            om_num_k = []
            for Oh in Oh_list:
                value_k_Oh = om_numerical(Oh, Bo, k, value_k_Oh)[0] 
                om_num_k.append(value_k_Oh)
            value_k_Oh = om_num_k[0]
            om_num.append(om_num_k)
        om_num = np.transpose(om_num,(2,0,1))
        np.save(file_name,om_num)
        
    #Numerical complex pulsation
    om_num = np.load(file_name)

    #Analytical error
    err_lub = err_norm(
        np.array([[om_lub(Oh, Bo, k) for Oh in Oh_list] for k in k_list]),
        np.array([[0 for Oh in Oh_list] for k in k_list]),
        om_num)
    err_visc = err_norm(
        np.array([[om_normal_mode_viscous(Oh, Bo, k) for Oh in Oh_list] for k in k_list]),
        np.array([[0 for Oh in Oh_list] for k in k_list]),
        om_num)
    err_in = err_norm(
        np.array([[om_normal_mode_inertial(Oh, Bo, k) for Oh in Oh_list] for k in k_list]),
        np.array([[puls_normal_mode_inertial(Oh, Bo, k) for Oh in Oh_list] for k in k_list]),
        om_num)
  
    #Figure parameter and contour's labels
    plt.figure(figsize=(5, 4))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$Oh$', family = "Roboto", weight="ultralight")      
    plt.ylabel('$k$', family = "Roboto", weight="ultralight")
    
    fmt = {}
    for l, s in zip([0.005, 0.05, 0.2], ['0.5 \%', '5 \%', '20 \%']):
        fmt[l] = s
        
    #Plot contour lines and fillings
    for err, c in zip([err_visc, err_lub, err_in],['red', 'grey', 'blue']):
        plt.contourf(Oh_list, k_list, err, levels = [-0.2, 0.2], colors = c, alpha = 0.2);
        cs = plt.contour(Oh_list, k_list, err, levels = [0.005, 0.05, 0.2], colors = c);
        plt.clabel(cs, fmt=fmt, fontsize=10)
    x = [pulsation(Bo, k)/(k**2/1.3115+1/0.732) for k in k_list]
    plt.plot(x, k_list, linewidth = 1.5, c = 'black')

Oh_list = np.logspace(-3, 1, 15)
k_list = np.logspace(-2, 2, 15)
Bo = 1
plotErrorOm (Oh_list, k_list, Bo, 'fig3_om_num.npy', True)
plt.tight_layout(pad=1.)
plt.savefig("figure3.pdf")

# #%% Visu_Figure 3
# # Not for the article : vue of the curve-fitting and comparison with models for different k, Oh.
# def plotGrowtRate(Oh, Bo, k, ax): 
#     om_num, t_all, sampled_eta = om_numerical(Oh, Bo, k)
#     if (Oh < pulsation(Bo, k)/(k**2/0.7+1/0.6)):
#         om_relax = om_normal_mode_inertial(Oh, Bo, k)
#     else: 
#         om_relax = om_normal_mode_viscous(Oh, Bo, k)
#     sampled_t = abs(t_all*om_relax)
    
#     ax.set_title(np.round(om_relax/om_num[0]-1,5))
#     ax.plot(sampled_t, np.abs(sampled_eta), 'black', label = r'Numerical resolution')
#     ax.plot(sampled_t, np.exp(- t_all * om_num[0]), 'red', label = 'Decaying')
#     ax.plot(sampled_t, np.abs(np.exp(- om_num[0] * t_all)*np.cos(om_num[1] * t_all)), 'blue', label = 'Decaying')
#     ax.set_ylim([0,1])
#     return om_num

# Bo = 1
# Oh = np.logspace(-3, 0, 4)
# k = np.logspace(-2, 2, 5)

# fig, ax = plt.subplots(ncols = len(Oh), nrows = len(k), figsize=(9, 9))

# om_num = [[0,pulsation(Bo, k0)] for k0 in k]
# for l in range(len(Oh)):
#     om_num = [plotGrowtRate(Oh[l], Bo, k[i], ax[len(k)-1-i,l]) for i in range(len(k))]
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
    om_gwr_Oh = np.load("RayleighTaylor")

om_lub_Oh = []
for Oh in Oh_list:
    om_lub_Oh.append([np.abs(om_normal_mode_viscous(Oh, Bo, k)) for k in k_list2])
om_potential = [pulsation(Bo, k) for k in k_list]

Colors = ['orange', 'green', 'black']

plt.figure()
plt.xlabel(r'k')
plt.ylabel(r'$\omega$')
plt.plot(k_list, om_potential, lw=1.0, alpha = 0.4, color = Colors[-1], label = r'Potential')
for Oh, om_gwr, om_lub, c in zip(Oh_list, om_gwr_Oh, om_lub_Oh, Colors):
    plt.plot(k_list, np.abs(om_gwr), '--', lw=1.0, color = c, alpha = 0.8, label = r'Cortelezzi resolution, Oh = ' + str(Oh))
    plt.plot(k_list2, om_lub, '-', lw=1.0, alpha = 0.4, color = c, label = 'Lubrication, Oh = ' + str(Oh))
plt.legend()
plt.tight_layout(pad=0.)