# -*- coding: utf-8 -*-
"""
@author: Clément & Arnaud
"""

import numpy as np
from mpmath import mp, j
from Functions import freeSurface, om_lub, pulsation, om_analytic, om_numerical, decaying_sinusoid, err_norm
from Functions import om_normal_mode_viscous, puls_normal_mode_inertial, om_normal_mode_inertial


############################# Figure 1 ########################################
# This function returns height evolution with time with three different approaches : 
#    * Complete resolution using Cortelezzi's initial value resolution
#    * Analytical resolution by normal mode analysis
#    * Lubrication approximation
# These approaches are compared for two cases on Figure 1.

def compute_datas_fig1(Oh, Bo, k):
    om_lub_relax = om_lub(Oh, Bo, k)
    om_0 = pulsation(Bo, k)
    om_ana = om_analytic(Oh, Bo, k)/om_lub_relax
    
    t_all = np.linspace(0.0001, 1., 300) * max(abs(5./om_lub_relax), 13./om_0)
    sampled_t = abs(t_all*om_lub_relax)
    sampled_eta = freeSurface(t_all, Oh, Bo, k)
    sampled_eta_lub = np.exp(-t_all*om_lub_relax)
    sampled_eta_ana = decaying_sinusoid(sampled_t, float(-mp.re(om_ana)), 
                                        float(mp.im(om_ana)))
                                                
    return sampled_t, sampled_eta, sampled_eta_lub, sampled_eta_ana

############################# Figure 2 ########################################
# For a given couple (Bo, k), this function returns the analytical value of the 
# complexe growth rate (obtained by normal mode analysis) for a range of Oh values.
def compute_datas_fig2(Oh_list, Bo, k):
    om_ana = []
    root_denom = j
    
    for Oh in Oh_list:
        root_denom = om_analytic(Oh, Bo, k, root_denom)
        om_ana.append([float(mp.re(root_denom)), float(mp.im(root_denom))])
    om_ana = np.array(om_ana)
    om_0 = pulsation(Bo,k)
    
    return om_ana, om_0

############################# Figure 3 ########################################

def compute_om_num(Oh_list, k_list, Bo):
    om_num = []
    value_k_Oh = (om_normal_mode_inertial(Oh_list[0], Bo, k_list[0]),
                   puls_normal_mode_inertial(Oh_list[0], Bo, k_list[0]), 1., 0.)
    for k in k_list:
        om_num_k = []
        for Oh in Oh_list:
            value_k_Oh = om_numerical(Oh, Bo, k, value_k_Oh)[0] 
            om_num_k.append(value_k_Oh)
        value_k_Oh = om_num_k[0]
        om_num.append(om_num_k)
    om_num = np.transpose(om_num,(2,0,1))
    
    return om_num

#Compare the different models for a range of Oh and k.
def compute_datas_fig3 (Oh_list, k_list, Bo, file_name, compute = False):
    #The data can be easily recompute but it takes about 1h.
    #For time efficiency, numerical values are by default taken in the txt file. 
    if compute:
        om_num = compute_om_num(Oh_list, k_list, Bo)
        np.save(file_name,om_num)
    else:
        om_num = np.load(file_name)
    
    # Analytical growth rates and pulsations
    om_lub_list = [[om_lub(Oh, Bo, k) 
                      for Oh in Oh_list] for k in k_list]
    om_visc_list = np.array([[om_normal_mode_viscous(Oh, Bo, k) 
                      for Oh in Oh_list] for k in k_list])
    om_inert_list = np.array([[om_normal_mode_inertial(Oh, Bo, k) 
                      for Oh in Oh_list] for k in k_list])
    puls = np.array([[puls_normal_mode_inertial(Oh, Bo, k) 
                      for Oh in Oh_list] for k in k_list])

    # Calculation of the errors and splitline
    err_lub = err_norm(om_lub_list, np.zeros_like(om_lub_list), om_num)
    err_visc = err_norm(om_visc_list, np.zeros_like(om_visc_list), om_num)
    err_in = err_norm(om_inert_list, puls, om_num)
    splitline = np.array([[pulsation(Bo, k)/(k**2/1.3115+1/0.732),k] for k in k_list])    

    return err_lub, err_visc, err_in, splitline
