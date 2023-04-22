# -*- coding: utf-8 -*-
"""
@author: Clément & Arnaud
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import gridspec, ticker
from Colors import cmap_amber, cmap_lblue, cmap_pink, colors
from Colors import cols_amber_rgba, cols_dpurple_rgba, cols_lblue_rgba, cols_indigo_rgba

########################### Initialisation ####################################
def fig_init():
    ## Parameters figures
    #The fonts Roboto and Roboto Condensed must be installed
    p = plt.rcParams
    p['text.usetex'] = False
    p['font.family'] = 'sans-serif'
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
    plt.rc('savefig', bbox='tight', transparent=True, dpi=300) 

    #Old parameters for Latex
    # plt.rcParams['text.usetex'] = True
    # plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amssymb} \usepackage[squaren,Gray]{SIunits} \usepackage{nicefrac}'
    # plt.rcParams['font.family'] = 'serif'
    # plt.rcParams['font.serif'] = 'cm'


def plot_fig2(Oh_list, Bo, k, om_ana, om_0):
    ## Creation of the colormap.
    #Two colormap are sampled to be printed on each side of the split point
    split = int(np.sum(om_ana[:,1]<0.018)/len(om_ana)*256)
    colors1 = plt.cm.YlOrBr_r(np.linspace(0.1, 0.9, 256-split))
    colors2 = plt.cm.Blues(np.linspace(0.3, 1, split))
    # These are then combined and used to build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
    ## Ploting datas
    plt.figure(constrained_layout=False, figsize = (8,4))
    [plt.scatter(0, 1, label = 'Natural pulsations', marker = 'P', s = 80, c = 'black'),
         plt.scatter(-0.93, 0, label = 'Split point', marker = '*', s = 60, c = 'black')]
    plt.scatter(0, -1, marker = 'P', s = 80, c = 'black')
    for i in [1, -1]:
        plt.arrow(-0.91, i*0.05, 0.05, i*0.2, head_width = 0.02, color = 'black')
        plt.scatter(om_ana[:,0]/om_0, i*om_ana[:,1]/om_0, s = 20, c = Oh_list, 
                    cmap=mymap, norm=mcolors.LogNorm())
    
    ## Axes titles
    plt.xlabel('$\omega_{relax}/\omega_0 = \Re(s/\omega_0)$', usetex=True)      
    plt.ylabel('$\omega_{osc}/\omega_0 = \Im(s/\omega_0)$', usetex=True)
    plt.colorbar(label = 'Oh')
    plt.legend()
    plt.tight_layout(pad=1.)
    
    ##Save figure
    plt.savefig("figure2.pdf")
    
############################# Figure 3 ########################################
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects as path_effects

def interpolate(X, Y, T):
    dR = (np.diff(X) ** 2 + np.diff(Y) ** 2) ** 0.5
    R = np.zeros_like(X)
    R[1:] = np.cumsum(dR)
    return np.interp(T, R, X), np.interp(T, R, Y), R[-1]

def add_text(ax, x, y, string, va="center", ha="center", size = 18, usetex = False):
    text = ax.text(
    x,
    y,
    string,
    va=va,
    transform=ax.transAxes,
    size=size,
    ha=ha,
    color="black",
    usetex = usetex
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), 
         path_effects.Normal()]
    )

def simple_text(ax, x, y, string):
    ax.text(
    x,
    y,
    string,
    va="top",
    size=12,
    ha="right",
    usetex = True
    )
    
def draw_split_line(ax, splitpoints, color):
    DC_to_FC = ax.transData.transform
    FC_to_DC = ax.transData.inverted().transform
    NDC_to_FC = ax.transAxes.transform
    FC_to_NDC = ax.transAxes.inverted().transform
    DC_to_NDC = lambda x: FC_to_NDC(DC_to_FC(x))
    NDC_to_DC = lambda x: FC_to_DC(NDC_to_FC(x))
    splitpointsNDC=DC_to_NDC(splitpoints)
        
    path = TextPath(
        (0, -0.75), "SPLIT LINE", prop=FontProperties(size=10, weight="bold")
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
    X, Y, _ = interpolate(splitpointsNDC[:,0], splitpointsNDC[:,1], 
                          np.linspace(0, vert[:, 0].min()-0.155, 10))
    splitpointsNDC_ = np.array([X,Y]).T
    splitpointsfirsthalf = NDC_to_DC(splitpointsNDC_)
    plt.plot(
        splitpointsfirsthalf[:,0],splitpointsfirsthalf[:,1], color=color, 
        linewidth=2, markersize=5, marker="o", markevery=[0, -1], 
        linestyle="--", dash_capstyle="round"
    )
    X, Y, _ = interpolate(splitpointsNDC[:,0], splitpointsNDC[:,1], 
                          np.linspace(vert[:, 0].max() - .1, D , 200))
    splitpointsNDC_ = np.array([X,Y]).T
    splitpointssecondhalf = NDC_to_DC(splitpointsNDC_)
    plt.plot(
        splitpointssecondhalf[:,0],splitpointssecondhalf[:,1], color=color, 
        linewidth=2, markersize=5, marker="o", markevery=[0, -1], 
        linestyle="--", dash_capstyle="round"
    )
    # Faint outline
    patch = PathPatch(path, color="white", zorder=10,
        alpha=0.25, linewidth=1.25, transform = ax.transAxes)
    ax.add_artist(patch)
    # Actual text
    patch = PathPatch(path, color=color, 
        zorder=30, linewidth=0.0, transform = ax.transAxes)
    ax.add_artist(patch)

#Hexagonal plot
#trying to use https://stackoverflow.com/questions/15140072/how-to-map-number-to-color-using-matplotlibs-colormap
def hexagonal_plot(Oh_full, k_full, err, cmap, lw=1):
    return plt.hexbin(Oh_full.flatten(), k_full.flatten(), 
                         C = err.flatten(), gridsize=(21,12), 
                         norm=mcolors.LogNorm(vmax=1.,vmin=0.0001),
                         cmap=cmap,edgecolor='white',
                         linewidths=lw,xscale = 'log', yscale = 'log', 
                         reduce_C_function=np.mean)


def draw_scheme(ax):
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
    ax.scatter(Zx, Zy, s=50, zorder=20, edgecolor="black", facecolor=colors["amber"][2], linewidth=0.5)
    Zx = [-6.]
    Zy = [-5.]
    ax.scatter(Zx, Zy, s=50, zorder=20, edgecolor="black", facecolor=colors["d.orange"][2], linewidth=0.5)
    # note, the points lie on the line 0.5 x - 2

   
    simple_text(ax, -2.8, -4.3, r"$\Delta \omega$")
    simple_text(ax, -5.7, -5.5, r"$\omega$")
    simple_text(ax, -0.5, -1.5, r"$\omega_\mathrm{model}$")

    ax.annotate('', xy=(-5.8, 0.5*(-5.8)-2.), xytext=(-2.2, 0.5*(-2.2)-2.),
                arrowprops=dict(facecolor='black', arrowstyle='<->'))



def draw_subplot(ax, string, x, y):
    ax.set_xlim(0.001,10)
    ax.set_ylim(0.01,100)

    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction='in')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    add_text(ax, 0.05, 0.9, string, size=12, ha="left")

    add_text(ax, x,y, r"$$\frac{\left\|\Delta\omega\right\|}{\left\|\omega\right\|}$$",
        va="top", size=10, usetex=True)

def draw_colorbar(ax, hb):
    cb = plt.colorbar(hb, cax=ax,aspect=20/1, ticks=[1.5e-4, 1e-2, 7e-1])
    cb.ax.set_yticklabels(['0.01 %', '1 %', '100 %'])
    ax.tick_params(axis=u'both', which=u'both',length=0)
    for t in cb.ax.get_yticklabels():
        t.set_horizontalalignment('left')
        t.set_fontsize('8')
        
def draw_major_plot(ax, hb_main, hb_array, splitline):
    #Formatting
    ax.set_xlim(0.001,10)
    ax.set_ylim(0.01,100)
    ax.set_xticks(np.logspace(-3, 1, 4 + 1))
    ax.set_xlabel("Oh")
    ax.set_yticks(np.logspace(-2, 2, 4 + 1))
    ax.set_ylabel("k")
    
    #Add text
    add_text(ax, 0.9, 0.9, "modes")
    add_text(ax, 0.9, 0.95, "Damped")
    add_text(ax, 0.1, 0.9, "modes")
    add_text(ax, 0.1, 0.95, "Oscillating")
    
    #Add Split line
    draw_split_line(ax, splitline, colors["blue grey"][5])

    #Mege errors
    ax.figure.canvas.draw()
    cols=hb_main.get_facecolors()
    counter = 0
    for val_in, val_visc, val_lub in zip(hb_array[0], hb_array[1], hb_array[2]):
        if ((val_in < 10.) & (val_in < val_visc) & (val_in < val_lub)): 
            ind = np.clip((2*np.log10(val_in)+8.).astype(int),0, 7)
            cols[counter] = cols_amber_rgba[ind]
        elif ((val_visc < 10.) & (val_visc < val_in) & (val_visc > 0.1*val_lub) & (val_visc < val_lub)):
            ind = np.clip((2*np.log10(val_visc)+8.).astype(int),0, 7)
            cols[counter] = cols_indigo_rgba[ind]
        elif ((val_visc < 10.) & (val_visc < val_in) & (val_visc < val_lub)):
            ind = np.clip((2*np.log10(val_visc)+8.).astype(int),0, 7)
            cols[counter] = cols_lblue_rgba[ind]
        elif ((val_lub < 10.) & (val_lub < val_in) & (val_lub > 0.1*val_visc) & (val_lub < val_visc)):
            ind = np.clip((2*np.log10(val_lub)+8.).astype(int),0, 7)
            cols[counter] = cols_dpurple_rgba[ind]
        else:
            cols[counter] = [0, 0, 0, 1]
        counter+=1
        
        hb_main.set(array=None,facecolors=cols)
        
def plot_fig3(Oh_list, k_list, err_lub, err_visc, err_in, splitline):
    Oh_full = np.array([[Oh for Oh in Oh_list] for k in k_list])
    k_full = np.array([[k for Oh in Oh_list] for k in k_list])

    ## Initialise the figure
    fig = plt.figure(constrained_layout=False, figsize=(10.57, 8.3))
    gspec = gridspec.GridSpec(ncols=5, nrows=7, figure=fig, 
                              width_ratios=[83, .7, 20, 1., 1], 
                              height_ratios=[20, 1, 20, 1, 20, 1, 20], 
                              hspace=0., wspace=0.)

    ## Draw first subplot (scheme)
    ax = plt.subplot(gspec[0, 2],aspect=1)
    draw_scheme(ax)

    ## Draw subplots with hexagons for inertial, viscous and lubrication
    hb_array = []
    for i, err, cmap, string, x, y in zip([2, 4, 6], 
                                    [err_in, err_visc, err_lub], 
                                    [cmap_amber, cmap_lblue, cmap_pink],
                                    ["Inertial model", "Viscous model", "Lubrication theory"],
                                    [0.85, 0.35, 0.45],
                                    [0.95, 0.8, 0.8]):

        ax = plt.subplot(gspec[i,2],aspect=1)
        hb = hexagonal_plot(Oh_full, k_full, err, cmap, 0.25)
        hb_array.append(hb.get_array())
        draw_subplot(ax, string, x, y)
        ax = plt.subplot(gspec[i,4])
        draw_colorbar(ax, hb)

    ## Draw Major subplot
    ax = plt.subplot(gspec[0:, 0],aspect=1)           
    hb_main = hexagonal_plot(Oh_full, k_full, err_in, cmap_amber)
    draw_major_plot(ax, hb_main, hb_array, splitline)
    
    ## Save figure
    plt.savefig("figure3.pdf")
