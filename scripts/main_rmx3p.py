'''
This file is the main file to run the code of Remix-3P
'''

# include the main library path (the parent folder) in the path environment variable
import os,sys
root_folder = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(root_folder)

# import the library as a package (defined in __init__.py) => function calls are done through the lpackage (eg om.solve_model)
import source as source

# import other packages
import numpy as np
import matplotlib.pyplot as plt

plt.close("all")

#%% Main

## Definitions of the plant 
# Parameters 
eps_is_exp = 0.50
eps_is_pp = 0.5
eta_elec = 0.9
fluid = 'water'
SH_expansion = 10 #K
RK_cycle = source.procF.Rankine_cycle(fluid, eps_is_exp = eps_is_exp, eps_is_pp = eps_is_pp, eta_e = eta_elec, SH_expansion = SH_expansion)

## Simulation of the Rankine 
T_cd = 100 + 273.15 #K 
T_out_boiler = 220 + 273.15 #K 
W_dot_exp_el = 10000 #W

RK_cycle.solve_Rankine(T_cd = T_cd, T_boiler = T_out_boiler, W_dot_exp_el = W_dot_exp_el)
RK_cycle.print_results_Rankine()
RK_cycle.plot_Ts_ph_diagrams()


## Volumetric machine simulation 

#` Airsquared P17H043D-AC-LC-SH
scroll_exp = source.procF.expander(r_v_in = 3.64, V_cp = 134.5*1e-6, N_max = 3450)
scroll_exp.speed_intervals(RK_cycle.V_dot_wf)




#%%
# Simulution with varying T_boiler

# T_boiler_vect = np.linspace(200,400, 500) + 273.15
# RK_cycle.simulation_T_boiler(T_boiler_vect, W_dot_exp_el = W_dot_exp_el, T_cd = T_cd )
# RK_cycle.plot_T_boiler()




