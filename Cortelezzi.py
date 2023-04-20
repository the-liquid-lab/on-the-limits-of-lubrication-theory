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
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patheffects as path_effects
from matplotlib import gridspec, ticker
from matplotlib.colors import LogNorm
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
from mpmath import mp, findroot, j 
from scipy.optimize import curve_fit

# Ensure working directory
import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

from Functions import denom, freeSurface, om_lub, pulsation, decaying_sinusoid, err_norm, om_numerical
from Functions import om_normal_mode_viscous, puls_normal_mode_inertial, om_normal_mode_inertial, om_analytic

#Colors
import Colors
import vapeplot
clrs = vapeplot.palette('vaporwave')
clrlub=clrs[2]
clrpole=clrs[6]



## Parameters figures
#The font Roboto and Roboto Condensed must be installed
p = plt.rcParams
p['text.usetex'] = False
p['font.family'] = 'sans-serif'
p["figure.figsize"] = 10.57, 8.3
p["font.sans-serif"] = ["Roboto Condensed"]
p["font.weight"] = "light"
p["ytick.minor.visible"] = True
p["xtick.minor.visible"] = True

#Figure parameter and contour's labels
plt.rc('font', size=12)  # general font size
plt.rc('axes', labelsize=11, titlesize=10, linewidth=1.)
plt.rc('lines', markersize=8, markeredgewidth=0., linewidth=0.4)
plt.rc('xtick',  labelsize=12, direction='in', bottom='true', top='true')
plt.rc('ytick',  labelsize=12, direction='in', left='true', right='true')
#plt.rc('legend', frameon=False, fancybox=False, numpoints=1, markerscale=1, 
#       fontsize=10, handlelength=0.6, handletextpad=0.6, labelspacing=0.3)
plt.rc('savefig', bbox='tight', transparent=True, dpi=300) 

#Old parameters for Latex
# plt.rcParams['text.usetex'] = True
# plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amssymb} \usepackage[squaren,Gray]{SIunits} \usepackage{nicefrac}'
# plt.rcParams['font.family'] = 'serif'
# plt.rcParams['font.serif'] = 'cm'

#%% Figure 1
#Comparison between lubrication, analytical and numerical results for 2 different situations : oscillations and relaxation

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
    y_major = ticker.LogLocator(base = 10.0, numticks = 5)
    axinset.yaxis.set_major_locator(y_major)
    y_minor = ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    axinset.yaxis.set_minor_locator(y_minor)
    axinset.yaxis.set_minor_formatter(ticker.NullFormatter())
    axinset.set_xticks(sampled_t[-1]*np.arange(5)/4)

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
plt.scatter(om_ana[:,0]/om_0, om_ana[:,1]/om_0, s = 20, c = Oh_list, cmap=mymap, norm=LogNorm())
plt.scatter(om_ana[:,0]/om_0, -om_ana[:,1]/om_0, s = 20, c = Oh_list, cmap=mymap, norm=LogNorm())
plt.xlabel('$\omega_{relax}/\omega_0 = \Re(s/\omega_0)$', family = "Roboto", weight="ultralight")      
plt.ylabel('$\omega_{osc}/\omega_0 = \Im(s/\omega_0)$', family = "Roboto", weight="ultralight")
cbar = plt.colorbar(label = 'Oh')
plt.legend()
plt.tight_layout(pad=1.)
plt.savefig("figure2.pdf")

#%% Figure 3

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
        
    Oh_full = np.array([[Oh for Oh in Oh_list] for k in k_list])
    k_full = np.array([[k for Oh in Oh_list] for k in k_list])

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

    
    def interpolate(X, Y, T):
        dR = (np.diff(X) ** 2 + np.diff(Y) ** 2) ** 0.5
        R = np.zeros_like(X)
        R[1:] = np.cumsum(dR)
        return np.interp(T, R, X), np.interp(T, R, Y), R[-1]

    fig = plt.figure(constrained_layout=False)
    nrows, ncols = 7, 5
    w0, w1, w2, w3, w4 = 83, .7, 20, 1., 1
    h1, h2, h3, h4, h5, h6, h7 = 20, 1, 20, 1, 20, 1, 20
    gspec = gridspec.GridSpec(ncols=ncols, nrows=nrows, figure=fig, width_ratios=[w0, w1, w2, w3, w4], height_ratios=[h1,h2,h3,h4,h5,h6,h7], hspace=0., wspace=0.)
    
    ax = plt.subplot(gspec[0:, 0],aspect=1)
    ax.set_xlim(0.001,10)
    ax.set_ylim(0.01,100)
    ax.set_xticks(np.logspace(-3, 1, 4 + 1))
    ax.set_xlabel("Oh")
    ax.set_yticks(np.logspace(-2, 2, 4 + 1))
    ax.set_ylabel("k")
    hb_main = plt.hexbin(Oh_full.flatten(), k_full.flatten(), C = err_in.flatten(), gridsize=(21,12), norm=LogNorm(vmax=1.,vmin=0.0001),cmap=Colors.cmap_amber,edgecolor='white',linewidths=1,xscale = 'log', yscale = 'log', reduce_C_function=np.mean)
    #trying to use https://stackoverflow.com/questions/15140072/how-to-map-number-to-color-using-matplotlibs-colormap
    ax.figure.canvas.draw()
    cols=hb_main.get_facecolors()
    hb_main_array = hb_main.get_array()
    #print("the shape of hb_main_array is", hb_main_array.shape)
    #print("the shape of cols is", cols.shape)
    Colors.cols_ind = np.clip((2*np.log10(hb_main_array)+8.).astype(int),0, 7)
    x = [pulsation(Bo, k)/(k**2/1.3115+1/0.732) for k in k_list]
    splitpoints=np.array([x,k_list]).T
    DC_to_FC = ax.transData.transform
    FC_to_DC = ax.transData.inverted().transform
    NDC_to_FC = ax.transAxes.transform
    FC_to_NDC = ax.transAxes.inverted().transform
    DC_to_NDC = lambda x: FC_to_NDC(DC_to_FC(x))
    NDC_to_DC = lambda x: FC_to_DC(NDC_to_FC(x))
    splitpointsNDC=DC_to_NDC(splitpoints)

    path = TextPath(
        (0, -0.75), "SPLIT LINE", prop=FontProperties(size=10, family="Roboto Condensed", weight="bold")
    )
    vert = path.vertices
    vert.flags.writeable = True
    xmin, xmax = vert[:, 0].min(), vert[:, 0].max()
    ymin, ymax = vert[:, 1].min(), vert[:, 1].max()
    vert -= (xmin + xmax) / 2, (ymin + ymax) / 2
    vert *= 0.003
    X0, Y0, D = interpolate(splitpointsNDC[:,0], splitpointsNDC[:,1], .3 + vert[:, 0])
    X1, Y1, _ = interpolate(splitpointsNDC[:,0], splitpointsNDC[:,1], .3 + vert[:, 0] + 0.1)

    # Transform text vertices
    dX, dY = X1 - X0, Y1 - Y0
    norm = np.sqrt(dX ** 2 + dY ** 2)
    dX, dY = dX / norm, dY / norm
    X0 += -vert[:, 1] * dY
    Y0 += +vert[:, 1] * dX
    vert[:, 0], vert[:, 1] = X0, Y0
    X, Y, _ = interpolate(splitpointsNDC[:,0], splitpointsNDC[:,1], np.linspace(0, vert[:, 0].min()-0.155, 10))
    splitpointsNDC_ = np.array([X,Y]).T
    splitpointsfirsthalf = NDC_to_DC(splitpointsNDC_)
    plt.plot(
        splitpointsfirsthalf[:,0],splitpointsfirsthalf[:,1], color=Colors.colors["blue grey"][5], linewidth=2, markersize=5, marker="o", markevery=[0, -1], linestyle="--", dash_capstyle="round"
    )
    X, Y, _ = interpolate(splitpointsNDC[:,0], splitpointsNDC[:,1], np.linspace(vert[:, 0].max() - .1, D , 200))
    splitpointsNDC_ = np.array([X,Y]).T
    splitpointssecondhalf = NDC_to_DC(splitpointsNDC_)
    plt.plot(
        splitpointssecondhalf[:,0],splitpointssecondhalf[:,1], color=Colors.colors["blue grey"][5], linewidth=2, markersize=5, marker="o", markevery=[0, -1], linestyle="--", dash_capstyle="round"
    )
    # Faint outline
    patch = PathPatch(
        path,
        facecolor="white",
        zorder=10,
        alpha=0.25,
        edgecolor="white",
        linewidth=1.25,
        transform = ax.transAxes,
    )
    ax.add_artist(patch)
    # Actual text
    patch = PathPatch(
        path, facecolor=Colors.colors["blue grey"][5], zorder=30, edgecolor=Colors.colors["blue grey"][5], linewidth=0.0,    transform = ax.transAxes
    )
    ax.add_artist(patch)
    #ax.plot(x, k_list, linewidth = 4, c = Colors.colors["orange"][3], linestyle="--", dash_capstyle="round")
    text = ax.text(
        0.1,
        0.95,
        "Oscillating",
        va="center",
        transform=ax.transAxes,
        size=18,
        ha="center",
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )
    text = ax.text(
        0.1,
        0.9,
        "modes",
        va="center",
        transform=ax.transAxes,
        size=18,
        ha="center",
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )
    text = ax.text(
        0.9,
        0.95,
        "Damped",
        va="center",
        transform=ax.transAxes,
        size=18,
        ha="center",
        color="black"
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )
    text = ax.text(
        0.9,
        0.9,
        "modes",
        va="center",
        transform=ax.transAxes,
        size=18,
        ha="center",
        color="black"
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )

    ax = plt.subplot(gspec[0, 2],aspect=1)
    ax.tick_params(which="both", direction="in")
    ax.tick_params(which="both", right=True)
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.set_axisbelow(True)

    ax.set_xlim(-7, 3)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.00))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1.00))
    ax.set_xticklabels([])

    ax.set_ylim(-8, 2)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1.00))
    ax.set_yticklabels([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.grid(color=".9", linestyle="--")

    Zx = [-2.]
    Zy = [-3.]
    ax.scatter(Zx, Zy, s=50, zorder=20, edgecolor="black", facecolor=Colors.colors["amber"][2], linewidth=0.5)
    Zx = [-6.]
    Zy = [-5.]
    ax.scatter(Zx, Zy, s=50, zorder=20, edgecolor="black", facecolor=Colors.colors["d.orange"][2], linewidth=0.5)
    # note, the points lie on the line 0.5 x - 2

    ax.text(
        -2.8,
        -4.3,
        r"$\Delta \omega$",
        va="top",
        #    transform=ax.transAxes,
        size=12,
        ha="right",
        usetex=True,
    )

    ax.text(
        -5.7,
        -5.5,
        r"$\omega$",
        va="top",
        #    transform=ax.transAxes,
        size=12,
        ha="right",
        usetex=True,
    )

    ax.text(
        -0.5,
        -1.5,
        r"$\omega_\mathrm{model}$",
        va="top",
        #    transform=ax.transAxes,
        size=12,
        ha="right",
        usetex=True,
    )

    ax.annotate('', xy=(-5.8, 0.5*(-5.8)-2.), xytext=(-2.2, 0.5*(-2.2)-2.),
                arrowprops=dict(facecolor='black', arrowstyle='<->'))
                     
    ax = plt.subplot(gspec[2, 2],aspect=1)
    ax.set_xlim(0.001,10)
    ax.set_ylim(0.01,100)
    hb_in = plt.hexbin(Oh_full.flatten(), k_full.flatten(), C = err_in.flatten(), gridsize=(21,12), norm=LogNorm(vmax=1.,vmin=0.0001),cmap=Colors.cmap_amber,edgecolor='white',linewidths=.25,xscale = 'log', yscale = 'log', reduce_C_function=np.mean)
    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction='in')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    hb_in_array = hb_in.get_array()

    text = ax.text(
        0.05,
        0.9,
        "Inertial model",
        va="center",
        transform=ax.transAxes,
        size=12,
        ha="left",
        #    usetex=True,
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )
    ax.text(
        0.85,
        0.95,
        r"$$\frac{\left\|\Delta\omega\right\|}{\left\|\omega\right\|}$$",
        va="top",
        transform=ax.transAxes,
        size=10,
        ha="center",
        usetex=True,
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )

    ax = plt.subplot(gspec[2, 4])
    cb = plt.colorbar(hb_in, cax=ax,aspect=20/1, ticks=[1.5e-4, 1e-2, 7e-1])
    cb.ax.set_yticklabels(['0.01 %', '1 %', '100 %'])
    ax.tick_params(axis=u'both', which=u'both',length=0)
    for t in cb.ax.get_yticklabels():
        t.set_horizontalalignment('left')
        t.set_fontsize('8')
    #cb.ax.yaxis.set_tick_params(pad=33)

    ax = plt.subplot(gspec[4, 2],aspect=1)
    ax.set_xlim(0.001,10)
    ax.set_ylim(0.01,100)
    hb_visc = plt.hexbin(Oh_full.flatten(), k_full.flatten(), C = err_visc.flatten(), gridsize=(21,12), norm=LogNorm(vmax=1.,vmin=0.0001),cmap=Colors.cmap_lblue,edgecolor='white',linewidths=.25,xscale = 'log', yscale = 'log', reduce_C_function=np.mean)
    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction='in')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    hb_visc_array = hb_visc.get_array()

    text = ax.text(
        0.05,
        0.9,
        "Viscous model",
        va="center",
        transform=ax.transAxes,
        size=12,
        ha="left",
        #    usetex=True,
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )
    ax.text(
        0.35,
        0.8,
        r"$$\frac{\left\|\Delta\omega\right\|}{\left\|\omega\right\|}$$",
        va="top",
        transform=ax.transAxes,
        size=10,
        ha="center",
        usetex=True,
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )

    ax = plt.subplot(gspec[4, 4])
    cb = plt.colorbar(hb_visc, cax=ax,aspect=20/1, ticks=[1.5e-4, 1e-2, 7e-1])
    cb.ax.set_yticklabels(['0.01 %', '1 %', '100 %'])
    ax.tick_params(axis=u'both', which=u'both',length=0)
    for t in cb.ax.get_yticklabels():
        t.set_horizontalalignment('left')
        t.set_fontsize('8')

    ax = plt.subplot(gspec[6, 2],aspect=1)
    ax.set_xlim(0.001,10)
    ax.set_ylim(0.01,100)
    hb_lub = plt.hexbin(Oh_full.flatten(), k_full.flatten(), C = err_lub.flatten(), gridsize=(21,12), norm=LogNorm(vmax=1.,vmin=0.0001),cmap=Colors.cmap_pink,edgecolor='white',linewidths=.25,xscale = 'log', yscale = 'log', reduce_C_function=np.mean)
    plt.tick_params(
        axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    direction='in')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    hb_lub_array = hb_lub.get_array()

    text = ax.text(
        0.05,
        0.9,
        "Lubrication theory",
        va="center",
        transform=ax.transAxes,
        size=12,
        ha="left",
        #    usetex=True,
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )
    ax.text(
        0.45,
        0.8,
        r"$$\frac{\left\|\Delta\omega\right\|}{\left\|\omega\right\|}$$",
        va="top",
        transform=ax.transAxes,
        size=10,
        ha="center",
        usetex=True,
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
    )

    ax = plt.subplot(gspec[6, 4])
    cb = plt.colorbar(hb_lub, cax=ax,aspect=20/1, ticks=[1.5e-4, 1e-2, 7e-1])
    cb.ax.set_yticklabels(['0.01 %', '1 %', '100 %'])
    ax.tick_params(axis=u'both', which=u'both',length=0)
    for t in cb.ax.get_yticklabels():
        t.set_horizontalalignment('left')
        t.set_fontsize('8')


    counter = 0
    for val_in, val_visc, val_lub in zip(hb_in_array, hb_visc_array, hb_lub_array):
        if ((val_in < 10.) & (val_in < val_visc) & (val_in < val_lub)): 
            ind = np.clip((2*np.log10(val_in)+8.).astype(int),0, 7)
            cols[counter] = Colors.cols_amber_rgba[ind]
        elif ((val_visc < 10.) & (val_visc < val_in) & (val_visc > 0.1*val_lub) & (val_visc < val_lub)):
            ind = np.clip((2*np.log10(val_visc)+8.).astype(int),0, 7)
            cols[counter] = Colors.cols_indigo_rgba[ind]
        elif ((val_visc < 10.) & (val_visc < val_in) & (val_visc < val_lub)):
            ind = np.clip((2*np.log10(val_visc)+8.).astype(int),0, 7)
            cols[counter] = Colors.cols_lblue_rgba[ind]
        elif ((val_lub < 10.) & (val_lub < val_in) & (val_lub > 0.1*val_visc) & (val_lub < val_visc)):
            ind = np.clip((2*np.log10(val_lub)+8.).astype(int),0, 7)
            cols[counter] = Colors.cols_dpurple_rgba[ind]
        else:
            cols[counter] = [0, 0, 0, 1]
        counter+=1
        
        hb_main.set(array=None,facecolors=cols)


Oh_list = np.logspace(-3, 1, 50)
k_list = np.logspace(-2, 2, 50)
Bo = 1
plotErrorOm (Oh_list, k_list, Bo, 'fig3_om_num.npy', False)
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