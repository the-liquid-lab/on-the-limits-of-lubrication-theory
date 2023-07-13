# -*- coding: utf-8 -*-
"""
@author: Clément & Arnaud
"""
#Must install the fonts Roboto and Roboto condensed
#Must install the packages gwr_inversion and mpmath

#%% 
# Import
import numpy as np
# Ensure working directory
import os
os.chdir(os.path.abspath(os.path.dirname(__file__)))

from Compute_datas import datas_fig1, datas_fig2, datas_fig3, datas_fig4
from Plots import fig_init, plot_fig1, plot_fig2, plot_fig3, plot_fig4
fig_init()

#%% Figure 1

## Parameters
Bo = 0.001
Oh_list = [10, 0.01]
k_list = [0.1, 0.5]
all_datas = []

#Compute datas 
#t, eta, eta_lub, eta_ana, om_lub, om_0, om_ana_r, om_ana_i, split
for Oh, k in zip(Oh_list, k_list):
    all_datas.append(datas_fig1(Oh, Bo, k))

plot_fig1(Oh_list, k_list, all_datas)
#%% Figure 2

## Parameters
Bo = 0.001
k = 0.5
Oh_list = np.logspace(-3.5, 0.5, 600)

## Import datas
om_ana, om_0 = datas_fig2(Oh_list, Bo, k)

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
err_lub, err_visc, err_in, splitline = datas_fig3 (Oh_list, k_list, Bo, 'fig3_om_num.npy', ComputeAllFig3)

## Plot datas
plot_fig3(Oh_list, k_list, err_lub, err_visc, err_in, splitline)


#%% Figure 4 
#Rayleigh-Taylor
## Parameters
Bo = -0.5
Oh_list = [0.01, 1.]
k_list = np.linspace(0.005, 0.999, 100) * np.sqrt(-Bo)
k_list2 = np.linspace(0.005, 1., 100) * np.sqrt(-Bo)

## Import datas
om_gwr_Oh, om_potential, om_norm_in, om_lub_list, om_norm_visc = datas_fig4 (Oh_list, k_list, k_list2, Bo, "RayleighTaylor.npy")

## Plot datas
plot_fig4(Oh_list, k_list, k_list2, om_gwr_Oh, om_potential, om_norm_in, om_lub_list, om_norm_visc)