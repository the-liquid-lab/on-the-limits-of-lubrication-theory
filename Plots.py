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
clrpole = colors["teal"][6]
clrlub = colors["pink"][6]
clrlaplace = colors["amber"][7]
clrinertia = colors["blue"][1]


########################### Initialisation ####################################
def fig_init():
    ## Parameters figures
    #The fonts Roboto and Roboto Condensed must be installed
    p = plt.rcParams
    p["font.sans-serif"] = ["Roboto Condensed"]
    p["ytick.minor.visible"] = True
    p["xtick.minor.visible"] = True

    #Figure parameter and contour's labels
    plt.rc('font', size=10.5, weight = "light", family = "sans-serif")  # general font size
    plt.rc('axes', labelsize=10.5, titlesize=10.5, linewidth=1.)
    plt.rc('lines', markersize=8, markeredgewidth=0., linewidth=0.4)
    plt.rc('scatter', edgecolors = "black")
    plt.rc('xtick',  labelsize=10.5, direction='in', bottom='true', top='true')
    plt.rc('ytick',  labelsize=10.5, direction='in', left='true', right='true')
    plt.rc('legend',  fontsize=8)
    plt.rc('savefig', bbox='tight', transparent=True, dpi=300) 
    plt.rc('text',  usetex=False)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')

############################# Figure 1 ########################################
##Function to plot the evolution of the amplitude (top figure)
def plot_amplitude(Oh, k, t, eta, eta_lub, eta_ana, ax, hide_axis):
    ax.set_ylim([-0.05,1.05])
    ax.plot(t, np.abs(eta_ana), color=clrpole, ls=":", dash_capstyle="round", 
            linewidth=2, label = 'Analytical resolution')
    ax.plot(t, eta_lub, color=clrlub, ls=(0,(0.01,2)), 
            linewidth=2, dash_capstyle="round", label = 'Lubrication theory')
    ax.plot(t, np.abs(eta), color=clrlaplace, solid_capstyle="round", linewidth=1)
    
    ax.text(0.09, 0.9, r"Oh = " + str(Oh), transform=ax.transAxes)
    ax.text(0.14, 0.81, r"k = " + str(k), transform=ax.transAxes)
    
    if hide_axis:
        plt.setp(ax.get_yticklabels(), visible=False)
    else:
        textxlabel = r'($\times\tau_\mathregular{relax}$)'
        ax.text(0.80, -0.08, textxlabel, transform=ax.transAxes)
        ax.set_ylabel('Amplitude')
        plot_inset(t, eta, eta_lub, eta_ana, ax)

##Function to plot the error on the amplitude (inset in tht top left figure)
def plot_inset(t, eta, eta_lub, eta_ana, ax):
    ax = ax.inset_axes([0.57,0.57,0.4,0.4])
    ax.linewidth=0.5
    ax.set_ylabel('Error', size="x-small", labelpad = 1.0)
    ax.set_xlabel(r't/$\tau_\mathregular{relax}$', size="x-small", labelpad = 1.0)
    error_lub=np.subtract(eta,eta_lub)
    error_pole=np.subtract(eta,eta_ana)
    ax.semilogy(t, np.abs(error_lub), color=clrlub, solid_capstyle='round', linewidth=1.5)
    ax.semilogy(t, np.abs(error_pole), color=clrpole, solid_capstyle='round', linewidth=1.5)
    y_major = ticker.LogLocator(base = 10.0, numticks = 5)
    ax.yaxis.set_major_locator(y_major)
    y_minor = ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    ax.yaxis.set_minor_locator(y_minor)
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xticks(np.arange(6))
    ax.grid(True, axis='both', which='both', linewidth=0.125, alpha=0.5) # which='major',
    ax.set_axisbelow(True)
    ax.spines[['top','bottom','left','right']].set_linewidth(0.5)
    ax.patch.set_alpha(0.5)
    ax.tick_params(axis='both', which='major', labelsize=6, width=0.5)
    ax.tick_params(axis='both', which='minor', labelsize=6, width=0.25)
    ax.text(0.25, 0.3, "discrete pole", transform=ax.transAxes, 
                size = "x-small", weight="regular", color=clrpole)
    ax.text(0.15, 0.77, "lubrication", transform=ax.transAxes, 
                size = "x-small", weight="regular", color=clrlub)


##In this section, the functions are made to plot the bottom graph which is 
#separated into two parts

#Initialisation of each part and plot the branch cut
def branch_cut(ax, xmin, xmax, y, ad_split, xtrsl = 0., dilat = 1.):   
    ax.set_xlim(xmin, xmax)
    ax.set_xticks([])
    ax.set_ylim(-y, y)
    ax.set_yticks([])
    
    ax.spines[['top', 'right', 'bottom']].set_visible(False)
    if xmax > 0:
        ax.spines['left'].set_position(('data', 0))
    else:
        ax.spines['left'].set_visible(False)
    
    tax = np.linspace(max(ad_split, xmin+2e-5),xmax-1e-5,3)
    if ad_split < xmax:
        ax.plot (tax, 0. * tax, solid_capstyle='round', linewidth = 1)

    if ad_split > xmin:
        x0 = min(ad_split, xmax)
        tcut = np.linspace(xmin+xtrsl,x0,200)
        zcut = -0.05*y*np.sin(2.*np.pi*9.*dilat/(xmax-xmin)*(tcut-x0))
        ax.plot (tcut, zcut, color=clrlaplace, solid_capstyle='round', linewidth = 1.5)
        if ad_split  < xmax:
            ax.scatter([ad_split], [0], s=10, zorder=20, facecolor=clrlaplace, linewidth=0.5)

    #Draw the diagonal for separation (X position depends on left or right side)
    d = .03
    if xmax > 0:
        X = (d/2, d+d/2)
    else:
        X = (1-2*d, 1.)
    ax.plot(X, (0.5-d, 0.5+d), transform=ax.transAxes,
                  solid_capstyle='round', linewidth = 1, color='k', clip_on=False)

#Draw the points for the different omega
def draw_point(ax, Zl, clr, marker='o'):
        ax.scatter(Zl[0], Zl[1], s=20, zorder=10, facecolor="None", linewidth=0.75)
        ax.scatter(Zl[0], Zl[1], s=20, zorder=20, edgecolor="None", facecolor=clr, alpha=0.75)
 
def plot_om(ax_left, ax_right, ad_om_lub, ad_om_ana_r, ad_om_ana_i, i):
    Za = [[ad_om_ana_r,ad_om_ana_r], [ad_om_ana_i,-ad_om_ana_i]]
    draw_point(ax_right, Za, clrpole)
    Zl = [-ad_om_lub, 0.]
    Zo = [[0., 0.] , [1., -1.]]
    if i==0:
        draw_point(ax_right, Zl, clrlub)
    if i==1:
        draw_point(ax_left ,Zl, clrlub)
        ax_right.scatter(Zo[0], Zo[1], s=80, zorder=10, marker="*", facecolor=clrinertia, linewidth=0.25)


#Add the annotations with arrows and text
def arrow(ax, x1, y1, x2, y2, rad = "0.3"):
    ax.annotate("", (x1,y1), size="x-small", xytext=(x2, y2), 
                textcoords="offset points",
                arrowprops=dict(arrowstyle="->", linewidth=0.5, 
                                connectionstyle="arc3,rad=" + rad))

def add_arr_text(ax_left, ax_right, ad_om_lub, ad_om_ana_r, ad_om_ana_i, ad_split, i):
    kw_left = dict(size="x-small", transform=ax_left.transAxes)
    kw_right = dict(size="x-small", transform=ax_right.transAxes)
    #Add arrows and text 
    if i == 0:       
        arrow(ax_left, ad_split - 1e-4, 0.,-8, 15, rad="-0.3")
        arrow(ax_left, ad_split, 0.,-8, -20)
        ax_left.text(0., 0.7, "branch \ncut", **kw_left)
        ax_left.text(0.15, 0.23, "branch \npoint", **kw_left)
            
        arrow(ax_right, -ad_om_lub-1e-5, -0.05,-8, -15, rad="-0.3")
        ax_right.text(-0.02, 0.23, r'$\omega_\mathrm{lub}$', **kw_right, usetex = True)
        ax_right.text(0.35, 0.7, "pole", **kw_right)
        
        arrow(ax_right, ad_om_ana_r+1e-5, 0.05,12, 15, rad="-0.3")
        
    else:
        arrow(ax_left,-ad_om_lub-0.01, -0.1,-8, -15,rad="-0.3")
        ax_left.text(0.25, 0.23, r'$\omega_\mathrm{lub}$', **kw_left, usetex = True)
        arrow(ax_right,ad_om_ana_r-0.02, ad_om_ana_i, -20, 5,rad="-0.3")
        ax_right.text(0.3, 0.85, "pole", **kw_right)
        arrow(ax_right,-0.01, -1.01, -20, -5,rad="0.3")
        ax_right.text(0.25, 0.08, "pulsation", **kw_right)
        
    ax_right.text(0.43, 0.98, r'$\mathfrak{I}(s/\omega_0)$', **kw_right, usetex = True)
    ax_right.text(0.76, 0.56, r'$\mathcal{R}(s/\omega_0)$', **kw_right, usetex = True) 
    
def plot_fig1(Oh_list, k_list, all_datas):
    ## Ploting datas, create different sections
    fig = plt.figure(constrained_layout=False, figsize = (5,3.9775))
    w0, w1, w2, w3, w4 = 0.82, 1.63, 0.1, 0.82, 1.63
    h1, h2, h3 = 2.45, 0.3, 1.225
    gspec = gridspec.GridSpec(ncols=5, nrows=3, figure=fig, 
                              width_ratios=[w0, w1, w2, w3, w4], 
                              height_ratios=[h1,h2,h3], hspace=0., wspace=0.)

    ## Two cases for the different (Oh, k) couples
    for Oh, k, datas, i in zip(Oh_list, k_list, all_datas, [0,1]):
        #Extract datas 
        t, eta, eta_lub, eta_ana, om_lub, om_0, om_ana_r, om_ana_i, split = datas
        ad_om_lub, ad_om_ana_r, ad_om_ana_i, ad_split = (om_lub, om_ana_r, om_ana_i, split)/om_0
    
        #Define and plot the top figure
        ax = fig.add_subplot(gspec[0,0+3*i:2+3*i])
        plot_amplitude(Oh, k, t, eta, eta_lub, eta_ana, ax, i)    

        #Define and plot the bottom figure       
        ax_left = plt.subplot(gspec[2, 0+3*i])
        ax_right = plt.subplot(gspec[2, 1+3*i])
        #Declare the limits of the bottom graph
        #Plot the branch cut, horizontal lines and separations
        y = 1.+i/2.
        if i == 0:
            x1,x2,x3,x4 = ad_split - 2e-4, ad_split + 1e-4, -4.5e-4, 1.5e-4
            branch_cut(ax_left, x1, x2, y, ad_split)
            branch_cut(ax_right, x3, x4, y, ad_split)
        else:
            x1,x2,x3,x4 = -9, -8.5, -0.75, 0.25
            branch_cut(ax_left, x1, x2, y, ad_split, xtrsl = 0.1)
            branch_cut(ax_right, x3, x4, y, ad_split, xtrsl = 0.03, dilat = 2.)

        plot_om(ax_left, ax_right, ad_om_lub, ad_om_ana_r, ad_om_ana_i, i)
        add_arr_text(ax_left, ax_right, ad_om_lub, ad_om_ana_r, ad_om_ana_i, ad_split, i)
    
    plt.tight_layout(pad=1.)
    plt.savefig("figure1.pdf")
    
############################# Figure 2 ########################################
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
    plt.figure(constrained_layout=False, figsize = (5,2.5))
    [plt.scatter(0, 1, label = 'Inviscid pulsations',  s=100, marker="*", c=clrinertia, linewidth=0.5),
         plt.scatter(-0.93, 0, label = 'Split point', marker = 'P', s = 80, c = clrlaplace, linewidth=0.5)]
    plt.scatter(0, -1,  s=100, marker="*", c=clrinertia, linewidth=0.5)
    for i in [1, -1]:
        plt.arrow(-0.91, i*0.05, 0.05, i*0.2, head_width = 0.02)
        plt.scatter(om_ana[:,0]/om_0, i*om_ana[:,1]/om_0, s = 10, c = Oh_list, 
                    edgecolor = 'None', cmap=mymap, norm=mcolors.LogNorm())
    
    ## Axes titles
    plt.xlabel('$\omega_{relax}/\omega_0 = \Re(s/\omega_0)$', usetex=True)      
    plt.ylabel('$\omega_{osc}/\omega_0 = \Im(s/\omega_0)$', usetex=True)
    plt.colorbar(label = 'Oh')
    plt.legend(loc = 3)
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

def add_text(ax, x, y, string, va="center", ha="center", size = 12, usetex = False):
    text = ax.text(
    x,
    y,
    string,
    va=va,
    transform=ax.transAxes,
    size=size,
    ha=ha,
    usetex = usetex
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2, foreground="white"), 
         path_effects.Normal()]
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
        (0, -0.75), "SPLIT LINE", prop=FontProperties(size = 11, weight="bold")
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
def hexagonal_plot(Oh_full, k_full, err, cmap, lw=0.5):
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
    ax.set_xticklabels([])

    ax.set_ylim(-8, 2)
    ax.set_yticklabels([])
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['bottom', 'left']].set_position(('data', 0))
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.grid(color=".9", linestyle="--")

    ax.scatter(-2, -3, s=20, zorder=10, facecolor=colors["amber"][2], linewidth=0.3)
    ax.scatter(-6, -5, s=20, zorder=10, facecolor=colors["d.orange"][2], linewidth=0.3)
    # note, the points lie on the line 0.5 x - 2

    ax.text(-4.2, -5.9, r"$\Delta \omega$", size=8, usetex = True)
    ax.text(-6, -7.1, r"$\omega$", size=8, usetex = True)
    ax.text(-5.2, -2.1, r"$\omega_\mathrm{model}$", size=8, usetex = True)

    ax.annotate('', xy=(-5.9, 0.5*(-5.9)-2.), xytext=(-2.1, 0.5*(-2.1)-2.),
                arrowprops=dict(facecolor='black', lw = 0.5, arrowstyle='<->'))

def draw_subplot(ax, string, x):
    ax.set_xlim(0.001,10)
    ax.set_ylim(0.01,100)

    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction='in')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    add_text(ax, 0.05, 0.85, string, size=8, ha="left")

    add_text(ax, x,0.7, r"$$\frac{\left\|\Delta\omega\right\|}{\left\|\omega\right\|}$$",
        va="top", size=8, usetex=True)

def draw_colorbar(ax, hb):
    cb = plt.colorbar(hb, cax=ax,aspect=15/1, ticks=[1.5e-4, 1e-2, 7e-1])
    cb.ax.set_yticklabels(['0.01 %', '1 %', '100 %'])
    ax.tick_params(axis=u'both', which=u'both',length=0)
    for t in cb.ax.get_yticklabels():
        t.set_horizontalalignment('left')
        t.set_fontsize('6')
        
def draw_major_plot(ax, hb_main, hb_array, splitline):
    #Formatting
    ax.set_xlim(0.001,10)
    ax.set_ylim(0.01,100)
    ax.set_xticks(np.logspace(-3, 1, 4 + 1))
    ax.set_xlabel(r"$Oh$", usetex = True)
    ax.set_yticks(np.logspace(-2, 2, 4 + 1))
    ax.set_ylabel(r"$k$", usetex = True)
    
    #Add text
    add_text(ax, 0.85, 0.9, "Damped\nmodes")
    add_text(ax, 0.15, 0.9, "Oscillating\nmodes")
    
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
    fig = plt.figure(constrained_layout=False, figsize=(5, 3.9))
    gspec = gridspec.GridSpec(ncols=5, nrows=7, figure=fig, 
                              width_ratios=[80, 2, 20, 0.5, 1], 
                              height_ratios=[18, 2, 18, 2, 18, 2, 18], 
                              hspace=0., wspace=0.)

    ## Draw first subplot (scheme)
    ax = plt.subplot(gspec[0, 2],aspect=1)
    draw_scheme(ax)

    ## Draw subplots with hexagons for inertial, viscous and lubrication
    hb_array = []
    for i, err, cmap, string, x in zip([2, 4, 6], 
                                    [err_in, err_visc, err_lub], 
                                    [cmap_amber, cmap_lblue, cmap_pink],
                                    ["Inertial model", "Viscous model", "Lubrication"],
                                    [0.75, 0.35, 0.35]):

        ax = plt.subplot(gspec[i,2],aspect=1)
        hb = hexagonal_plot(Oh_full, k_full, err, cmap, 0.1)
        hb_array.append(hb.get_array())
        draw_subplot(ax, string, x)
        ax = plt.subplot(gspec[i,4])
        draw_colorbar(ax, hb)

    ## Draw Major subplot
    ax = plt.subplot(gspec[0:, 0],aspect=1)           
    hb_main = hexagonal_plot(Oh_full, k_full, err_in, cmap_amber)
    draw_major_plot(ax, hb_main, hb_array, splitline)
    
    ## Save figure
    plt.savefig("figure3.pdf")

############################# Figure 4 ########################################
def plot_fig4(Oh_list, k_list, k_list2, om_gwr_Oh, om_potential, om_norm_in, om_lub_list, om_norm_visc):
    
    fig, ax = plt.subplots(1,2, figsize=(5, 4))
        
    ax[1].plot(k_list, om_potential, lw=1.0, alpha = 0.4, color = 'black', label = r'Potential')
    ax[1].plot(k_list, om_norm_in, '-', lw=1.0, alpha = 0.4, color = 'red', label = 'Normal mode')
    
    ax[0].set_ylabel(r'$\omega$', usetex = True)
    ax[0].plot(k_list2, om_lub_list, '-', lw=1.0, alpha = 0.4, color = 'blue', label = 'Lubrication')
    ax[0].plot(k_list2, om_norm_visc, '-', lw=1.0, alpha = 0.4, color = 'red', label = 'Normal mode')
    
    for Oh, axx, om_gwr in zip(Oh_list, [ax[1],ax[0]], om_gwr_Oh):
        axx.set_xlabel(r'$k$', usetex = True)
        axx.set_title('$Oh = ' + str(Oh) + '$', usetex = True)
        axx.plot(k_list, np.abs(om_gwr), '--', lw=1.0, color = 'orange', alpha = 0.8, label = r'Cortelezzi resolution')
        axx.legend(loc = 3)
    
    plt.tight_layout(pad=1.)
    ## Save figure
    plt.savefig("figure4.pdf")
