import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
# include the main library path (the parent folder) in the path environment variable
import os,sys
root_folder = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(root_folder)
# import the library as a package (defined in __init__.py) => function calls are done through the lpackage (eg om.solve_model)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times",
    "font.size": 14,
    "legend.fontsize": 12,              # Set legend font size
    "lines.linewidth": 1.5,             # Set default line width
    "figure.figsize": (4.5, 3.7),
    "figure.constrained_layout.use": True
})

#%% Functions and classes

class trigeneration_plant(): 
    """
    Class to create the trigeneration plant of the project Remix_3P
    """   
    def __init__(self):
        self.random_number = 5
        
class expander():
    """
    Class to create an expander for the Rankine cycle
    """    
    def __init__(self, r_v_in=3, V_cp=None, V_exp=None, N_max = 3000):
        
        
        self.eps_v = 0.9
        self.r_v_in = r_v_in
        self.N_max = N_max
        self.N_min = 1000
        # Exactly one must be provided
        if (V_cp is None) == (V_exp is None):
            raise ValueError("Provide exactly one of V_cp or V_exp.")
    
        if V_cp is not None:
            self.V_cp = V_cp
            self.V_exp = V_cp/r_v_in
        else:
            self.V_exp = V_exp
            self.V_cp = V_exp*r_v_in
    
    def speed_intervals(self, V_dot_in):
        N = 60 * V_dot_in * self.eps_v / self.V_exp
    
        print(f"{self.N_min} < {N:.1f} < {self.N_max}")
    
        if self.N_min < N < self.N_max:
            print("Speed is within the allowed range.")
        else:
            print("Speed is outside the allowed range.")
    
        return N
       

class Rankine_cycle():
    
    """
    Class to create the Rankine cycle of the project Remix_3P
    """
    def __init__(self, fluid, eps_is_exp = 0.7, eps_is_pp = 0.5, eta_e = 0.9, SH_expansion = 30):
        """
        Parameters
        ----------

        Returns
        -------
        """
        self.wf = fluid
        self.eps_is_exp = eps_is_exp
        self.eps_is_pp = eps_is_pp 
        self.eta_e = eta_e 
        
        self.SH_boiler = SH_expansion 
        self.SC_condenser = 5 # K
        self.HEX_pp = 5 # K
        
    def solve_Rankine(self, T_cd = 300, T_boiler = 470, W_dot_exp_el = 10000, T_air = 300): 
        
        # Condensing and evaporating temperatures and pressures
        # self.T_cd = T_air + self.HEX_pp + self.SC_condenser
        self.T_cd = T_cd
        self.T_ev = T_boiler - self.SH_boiler
        
        self.p_cd = PropsSI( 'P',  'T', self.T_cd, 'Q', 0.5, self.wf)
        self.p_ev = PropsSI( 'P',  'T', self.T_ev, 'Q', 0.5, self.wf)
        self.r_p = self.p_ev/self.p_cd
        
        ## Resolution
        
        # Expander
        self.W_dot_exp = W_dot_exp_el/self.eta_e
        self.h_su_exp = PropsSI( 'H',  'T', self.T_ev  + self.SH_boiler, 'P', self.p_ev, self.wf)
        self.s_su_exp = PropsSI( 'S',  'T', self.T_ev  + self.SH_boiler, 'P', self.p_ev, self.wf)
        self.h_ex_is_exp = PropsSI( 'H',  'S', self.s_su_exp, 'P', self.p_cd, self.wf)
        self.h_ex_exp = self.h_su_exp - self.eps_is_exp*(self.h_su_exp - self.h_ex_is_exp)
        self.m_dot_wf = self.W_dot_exp/(self.h_su_exp - self.h_ex_exp)
        self.V_dot_wf = self.m_dot_wf/PropsSI( 'D',  'T', self.T_ev  + self.SH_boiler, 'P', self.p_ev, self.wf)
        
        
        # Condenser
        self.h_su_pp = PropsSI( 'H',  'T', self.T_cd  - self.SC_condenser, 'P', self.p_cd, self.wf)
        self.Q_dot_cd = self.m_dot_wf*(self.h_ex_exp - self.h_su_pp)
        
        # Pump 
        self.s_su_pp = PropsSI( 'S',  'T', self.T_cd  - self.SC_condenser, 'P', self.p_cd, self.wf)
        self.h_ex_is_pp = PropsSI( 'H',  'S', self.s_su_pp, 'P', self.p_ev, self.wf)
        self.h_ex_pp = self.h_su_pp + (self.h_ex_is_pp - self.h_su_pp)/self.eps_is_pp
        self.W_dot_pp = self.m_dot_wf*(self.h_ex_pp - self.h_su_pp)
        
        # Evaporator 
        self.Q_dot_ev = self.m_dot_wf*(self.h_su_exp - self.h_ex_pp)
        
        # Energy balance 
        self.E_res = self.Q_dot_ev + self.W_dot_pp - self.Q_dot_cd - self.W_dot_exp
        self.eta = (self.W_dot_exp - self.W_dot_pp)/self.Q_dot_ev
        self.eta_C = 1 - T_air/T_boiler
        
        # Separation between the solar concentrator 
        self.h_ev_sat_liq = PropsSI( 'H',  'Q', 0, 'P', self.p_ev, self.wf)
        self.Q_dot_solar = self.m_dot_wf*(self.h_ev_sat_liq - self.h_ex_pp)
        self.Q_dot_boiler = self.m_dot_wf*(self.h_su_exp - self.h_ev_sat_liq)
        
    
    def property_calculation(self):
        # Data calculations 
        self.T_sat_low = PropsSI(  'T',  'P', self.p_cd, 'Q' , 0.5, self.wf)
        self.T_sat_high = PropsSI(  'T',  'P', self.p_ev, 'Q' , 0.5, self.wf)  
        self.s_liq_sat_low = PropsSI(  'S',  'P', self.p_cd, 'Q' , 0, self.wf)
        self.h_liq_sat_low = PropsSI( 'H',  'P', self.p_cd, 'Q' , 0, self.wf)
        self.s_vap_sat_low = PropsSI(  'S', 'P', self.p_cd, 'Q' , 1, self.wf)
        self.h_vap_sat_low = PropsSI(  'H',  'P', self.p_cd, 'Q' , 1, self.wf)
        self.s_liq_sat_high = PropsSI(  'S',  'P', self.p_ev, 'Q' , 0, self.wf)
        self.h_liq_sat_high = PropsSI(  'H',  'P', self.p_ev, 'Q' , 0, self.wf)
        self.s_vap_sat_high = PropsSI(  'S',  'P', self.p_ev, 'Q' , 1, self.wf)
        self.h_vap_sat_high = PropsSI(  'H',  'P', self.p_ev, 'Q' , 1, self.wf)
        
        
        # Point 1 su pp 
        self.T_su_pp = self.T_cd  - self.SC_condenser
        self.s_su_pp = PropsSI('S', 'P', self.p_cd, 'H', self.h_su_pp, self.wf)      
        self.p_su_pp = self.p_cd
        
        # Point 2 ex pp 
        self.p_ex_pp = self.p_ev
        self.s_ex_pp = PropsSI('S', 'P', self.p_ev, 'H', self.h_ex_pp, self.wf)     
        self.T_ex_pp = PropsSI('T', 'P', self.p_ev, 'H', self.h_ex_pp, self.wf)     
        
        # Point 3 ex solar plant 
        self.p_ex_solar = self.p_ev
        self.h_ex_solar = self.h_liq_sat_high
        self.s_ex_solar = self.s_liq_sat_high
        self.T_ex_solar = self.T_sat_high
        
        # Point 4 ex boiler 
        self.p_ex_boiler = self.p_ev
        self.h_ex_boiler = self.h_su_exp
        self.s_ex_boiler = self.s_su_exp
        self.T_ex_boiler = self.T_ev  + self.SH_boiler
        
        # Point 5 ex expander 
        self.h_ex_exp 
        self.p_ex_exp = self.p_cd
        self.s_ex_exp = PropsSI('S', 'P', self.p_cd, 'H', self.h_ex_exp, self.wf)     
        self.T_ex_exp = PropsSI('T', 'P', self.p_cd, 'H', self.h_ex_exp, self.wf)  
        self.Q_ex_exp = PropsSI('Q', 'P', self.p_cd, 'H', self.h_ex_exp, self.wf)     
        
    def print_results_Rankine(self):         
        print(
        f"p_ev               = {self.p_ev/1e5:8.4f} [bar]\n" 
        f"p_cd               = {self.p_cd/1e5:8.4f} [bar]\n" 
        f"r_p                = {self.p_ev/self.p_cd:8.2f} [-]\n" 
        f"m_dot              = {self.m_dot_wf*1000:8.2f} [g/s]\n"        
        f"eta (efficiency)   = {self.eta:8.2f} [-]\n"
        f"eta_C (Carnot)     = {self.eta_C:8.2f} [-]\n"
        f"W_dot_exp          = {self.W_dot_exp:8.0f} [W]\n"
        f"Q_dot_cd           = {self.Q_dot_cd:8.0f} [W]\n"
        f"W_dot_pp           = {self.W_dot_pp:8.0f} [W]\n"
        f"E_dot_res          = {self.E_res:8.2f} [W]\n"
        f"Q_dot_ev           = {self.Q_dot_ev:8.0f} [W]\n"
        f"Q_dot_solar        = {self.Q_dot_solar:8.0f} [W]\n"
        f"Q_dot_boiler       = {self.Q_dot_boiler:8.0f} [W]\n")
        
    def initialize_vectors(self, size):
        self.p_ev_vect = np.zeros(size)
        self.p_cd_vect = np.zeros(size)
        self.r_p_vect = np.zeros(size)
        self.m_dot_wf_vect = np.zeros(size)
        self.eta_vect = np.zeros(size)
        self.eta_C_vect = np.zeros(size)
        self.W_dot_exp_vect = np.zeros(size)
        self.W_dot_pp_vect = np.zeros(size)
        self.Q_dot_ev_vect = np.zeros(size)
        self.Q_dot_cd_vect = np.zeros(size)
        self.Q_dot_solar_vect = np.zeros(size)
        self.Q_dot_boiler_vect = np.zeros(size)

    def simulation_T_boiler(self, T_boiler_vect, W_dot_exp_el = 10000, T_cd = 300):
        N = len(T_boiler_vect)
        self.T_out_boiler_vect = T_boiler_vect
        self.initialize_vectors(N)

        for i in range(N): 
            self.solve_Rankine(T_cd = T_cd, T_boiler = T_boiler_vect[i], W_dot_exp_el = W_dot_exp_el)
            self.p_ev_vect[i] = self.p_ev
            self.p_cd_vect[i] = self.p_cd
            self.r_p_vect[i] = self.r_p
            self.m_dot_wf_vect[i] = self.m_dot_wf
            self.eta_vect[i] = self.eta
            self.eta_C_vect[i] = self.eta_C
            self.W_dot_exp_vect[i] = self.W_dot_exp
            self.W_dot_pp_vect[i] = self.W_dot_pp
            self.Q_dot_ev_vect[i] = self.Q_dot_ev
            self.Q_dot_cd_vect[i] = self.Q_dot_cd
            self.Q_dot_solar_vect[i] = self.Q_dot_solar
            self.Q_dot_boiler_vect[i] = self.Q_dot_boiler
            
    def plot_Ts_ph_diagrams(self):
        
        
        self.property_calculation()
    
        # Colors for the plots
        color1 = u'#1f77b4'
        color2 = u'#ff7f0e'
        color_blue = 'indigo'
        color_red = u'#d62728'
        color2 = 'darkorange'
        
        ############# T-s diagram 
        p_crit = PropsSI( 'PCRIT', self.wf)/1e5
        T_crit = PropsSI( 'TCRIT', self.wf)   # K
        
        
        # Generate the base of the diagram
        N_div = 300
        T_sat = np.linspace(20, (T_crit-0.1)-273.15 , N_div)
        T_sat_bell = np.linspace(0, (T_crit-0.1)-273.15 , N_div)
        s_sat_0 = np.zeros(len(T_sat))
        s_sat_1 = np.zeros(len(T_sat))
        for tt in range(len(T_sat_bell)):
            T_sat_iter = T_sat_bell[tt]
            s_sat_0[tt] = PropsSI(  'S', 'T' , T_sat_iter + 273.15, 'Q' , 0, self.wf)
            s_sat_1[tt] = PropsSI(  'S', 'T' , T_sat_iter + 273.15, 'Q' , 1, self.wf)  
        iso_P = np.array([0.1, 0.5, 2 ,5, 10, 20, 35,60])# bar
        s_range = np.linspace(500, 10000 , N_div) #J/kgK
        T_iso_P = np.zeros((len(iso_P), len(s_range)))
        T_sat_iso_P = np.zeros(len(iso_P))
        s_sat_liq = np.zeros(len(iso_P))
        s_sat_vap = np.zeros(len(iso_P))
        for pp in range(len(iso_P)):
            T_sat_iso_P[pp] = PropsSI(  'T', 'P' , iso_P[pp]*1e5, 'Q' , 0.5, self.wf)
            s_sat_liq[pp] = PropsSI(  'S', 'P' , iso_P[pp]*1e5, 'Q' , 0, self.wf)
            s_sat_vap[pp] = PropsSI(  'S', 'P' , iso_P[pp]*1e5, 'Q' , 1, self.wf)
            for ss in range(len(s_range)):
                if (s_range[ss] >= s_sat_liq[pp]) and (s_range[ss] <= s_sat_vap[pp]):
                    T_iso_P[pp][ss] = T_sat_iso_P[pp] - 273.15
                else:
                    T_iso_P[pp][ss] = PropsSI(  'T', 'P' , iso_P[pp]*1e5, 'S' , s_range[ss], self.wf) - 273.15
        iso_Q = np.linspace(0.2, 0.8, 4)# °C    
        s_iso_Q = np.zeros((len(iso_Q), len(T_sat)))
        for qq in range(len(iso_Q)):
            for tt in range(len(T_sat)):
                s_iso_Q[qq][tt] = PropsSI('S', 'Q' , iso_Q[qq], 'T', T_sat[tt] + 273.15, self.wf)
        T_sat_vect = np.concatenate((T_sat_bell, T_sat_bell[::-1]), axis=None)
        s_sat_vect = np.concatenate((s_sat_0, s_sat_1[::-1]), axis=None)
        
        ### Plot
        fig,ax = plt.subplots() 
        plt.xlim(300, 8800)
        plt.ylim(0, 300)
        for pp in range(len(   iso_P  )):
            plt.plot(s_range, T_iso_P[pp],'k--', linewidth = 1.2, zorder = 10)
        plt.plot(s_sat_vect,T_sat_vect,'k',linewidth = 2, zorder = 10)
        for qq in range(len(iso_Q)):
            plt.plot(s_iso_Q[qq], T_sat,'--', color = color1, linewidth = 0.8, zorder = 10)  
            plt.text(s_iso_Q[qq][0] -40, 3 , ''.join([str(round(iso_Q[qq],2))]), color = 'k', zorder = 5)  
        plt.xlabel('Entropy [J/(kg$\cdot$K)]')
        plt.ylabel('Temperature [°C]') 
        # legend = plt.legend(loc="upper left", framealpha=1) 
        # legend.get_frame().set_facecolor('white')
        # legend.get_frame().set_alpha(1)  
        # legend.set_zorder(50) 
        s_crit = PropsSI('S', 'T', T_crit, 'P', p_crit*1e5 ,self.wf) - 20
        for pp in range(len(iso_P)):
            plt.text(s_crit-120, T_sat_iso_P[pp] - 273.15 + 3, ''.join([str(iso_P[pp]),' bar']), color = 'k', fontsize = 14, zorder = 5)
        


        t1 = [self.T_su_pp , self.T_ex_pp, self.T_ex_solar, self.T_ev, self.T_ex_boiler, self.T_ex_exp, self.T_cd, self.T_su_pp] 
        s1 = [self.s_su_pp, self.s_ex_pp, self.s_ex_solar, self.s_vap_sat_high, self.s_ex_boiler, self.s_ex_exp,self.s_liq_sat_low, self.s_su_pp]
        t1 = [x -273.15 for x in t1]

        T_sat_vect = np.concatenate((T_sat, T_sat[::-1]), axis=None)
        s_sat_vect = np.concatenate((s_sat_0, s_sat_1[::-1]), axis=None)
         
        ##### Details on the plot 
        dotsize = 5
        fontsize = 16
        # Cycle line 
        plt.plot(s1,t1,'darkorange', linewidth = 1.5)
        # # Base points
        plt.plot(s1[0], t1[0], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(s1[0]+80, t1[0]-15,'1', color = 'black',fontsize = fontsize)
        plt.plot(s1[1], t1[1], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(s1[1]-120, t1[1]+5,'2', color = 'black',fontsize = fontsize)
        plt.plot(s1[2], t1[2], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(s1[2]-250, t1[2]+3,'3', color = 'black',fontsize = fontsize)
        plt.plot(s1[4], t1[4], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(s1[4], t1[4]+3,'4', color = 'black',fontsize = fontsize)
        plt.plot(s1[5], t1[5], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(s1[5], t1[5]+3,'5', color = 'black',fontsize = fontsize)


        # ########### P-h diagram 
        N_div = 5000
        P_sat = np.linspace(0.01, p_crit , N_div)
        h_sat_0 = np.zeros(len(P_sat))
        h_sat_1 = np.zeros(len(P_sat))
        for pp in range(len(P_sat)):
            P_sat_iter = P_sat[pp]
            h_sat_0[pp] = PropsSI('H', 'P', P_sat_iter*1e5, 'Q', 0, self.wf)
            h_sat_1[pp] = PropsSI('H', 'P', P_sat_iter*1e5, 'Q', 1, self.wf)

        iso_T = np.array([50,100,150,200,250])# °C
        P_range = np.linspace(0.03, 100, N_div)*1e2 # kPa
        h_iso_T = np.zeros((len(iso_T), len(P_range)))
        P_sat_iso_T = np.zeros(len(iso_T))
        h_sat_liq = np.zeros(len(iso_T))
        h_sat_vap = np.zeros(len(iso_T))
        for tt in range(len(iso_T)):
            #saturation pressure at the iso-T value
            P_sat_iso_T[tt] = PropsSI(  'P', 'T' , iso_T[tt] + 273.15, 'Q' , 0.5, self.wf)
            #useless but maybe useful in the future
            h_sat_liq[tt] = PropsSI(  'H', 'T' , iso_T[tt] + 273.15, 'Q' , 0, self.wf)
            h_sat_vap[tt] = PropsSI(  'H', 'T' , iso_T[tt] + 273.15, 'Q' , 1, self.wf)
            for pp in range(len(P_range)):
                if (P_range[pp] >= P_sat_iso_T[tt] - 0.1) and (P_range[pp] <= P_sat_iso_T[tt] + 0.1):
                    h_iso_T[tt][pp] = np.NaN
                else:
                    h_iso_T[tt][pp] = PropsSI(  'H', 'T' , iso_T[tt] + 273.15, 'P' , P_range[pp]*1e3, self.wf)
        h_iso_Q = np.zeros((len(iso_Q), len(P_range)))
        for qq in range(len(iso_Q)):
            for pp in range(len(P_range)):
                h_iso_Q[qq][pp] = PropsSI('H', 'Q' , iso_Q[qq], 'P' , P_range[pp]*1e3, self.wf)
        P_sat_vect = np.concatenate((P_sat, P_sat[::-1]), axis=None)
        h_sat_vect = np.concatenate((h_sat_0, h_sat_1[::-1]), axis=None)
        
        ###Plot
        plt.figure()
        plt.yscale('log')
        plt.ylim(0.01, 100)
        plt.xlim(50, 3000)
        # plt.xticks([200, 300, 400, 500])
        # # Force custom tick positions
        # plt.gca().yaxis.set_major_locator(FixedLocator([0.5, 1,2, 5, 10, 30]))
        # # Optional: make tick labels show as plain numbers instead of scientific notation
        # plt.gca().yaxis.set_major_formatter(ScalarFormatter())

        # Plot iso-T line
        for k in range(len(   iso_T  )):
            plt.plot(h_iso_T[k]/1000,P_range/100,'k--',linewidth = 1.2, zorder = 10)
            plt.text(1500, P_sat_iso_T[k]/1e5*1.07 , ''.join([str(iso_T[k]),' °C']), color = 'k', zorder = 5)
        # Plot iso-Q line  
        for qq in range(len(   iso_Q  )):
            plt.plot(h_iso_Q[qq]/1000,P_range/100,'--', color = color1, linewidth = 0.8, zorder = 10)    
            plt.text(h_iso_Q[qq][0]/1000 - 100, 0.018 , ''.join([str(round(iso_Q[qq],2))]), color = 'k', zorder = 5)
        # Plot saturation bell     
        plt.plot(h_sat_vect/1000,P_sat_vect,'k',linewidth = 2, zorder = 10, label = 'Pure ref. bell')
        plt.xlabel('Enthalpy [kJ/kg]')
        plt.ylabel('Pressure [bar]') 

        p1 = [self.p_su_pp/1e5, self.p_ex_pp/1e5, self.p_ex_solar/1e5, self.p_ex_boiler/1e5, self.p_ex_exp/1e5, self.p_su_pp/1e5] 
        h1 = [self.h_su_pp, self.h_ex_pp, self.h_ex_solar, self.h_ex_boiler, self.h_ex_exp, self.h_su_pp]     
        
        
        ##### Details on he plot 
        dotsize = 5
        fontsize = 16
        plt.plot(np.array(h1)/1000,p1,'darkorange', linewidth = 1.5)
         
        plt.plot(h1[0]/1000, p1[0], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(h1[0]/1000, p1[0]*1.5,'1', color = 'black',fontsize = fontsize)
       
        plt.plot(h1[1]/1000, p1[1], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(h1[1]/1000, p1[1]*1.1,'2', color = 'black',fontsize = fontsize)
              
        plt.plot(h1[2]/1000, p1[2], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(h1[2]/1000+100, p1[2]*1.1,'3', color = 'black',fontsize = fontsize)
       
        plt.plot(h1[3]/1000, p1[3], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(h1[3]/1000, p1[3]*1.1,'4', color = 'black',fontsize = fontsize)
        
        plt.plot(h1[4]/1000, p1[4], marker="o", markersize = dotsize, markeredgecolor = 'black', markerfacecolor="black")
        plt.text(h1[4]/1000, p1[4]*1.5,'5', color = 'black',fontsize = fontsize)
         
        
    def plot_T_boiler(self):
            
        color = [u'#00008B', u'#1f77b4', u'#17becf', u'#ff7f0e', u'#9DCA1C', u'#00004D'] # arc


        # Plot boiler power
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Boiler capacity [kW]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.Q_dot_boiler_vect/1000, color=u'#1f77b4')
        # plt.xticks([150, 175, 200,225, 250])
        
        # Pump
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Pump cons. [kW]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.W_dot_pp_vect/1000, color=u'#1f77b4')
        # plt.xticks([150, 175, 200,225, 250])
        
        # Plot boiler power
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Solar capacity [kW]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.Q_dot_solar_vect/1000, color=u'#1f77b4')
        # plt.xticks([150, 175, 200,225, 250])
        
        # Mass flow rate 
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('WF mass flow rate [g/s]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.m_dot_wf_vect*1000, color=u'#1f77b4')
        # plt.xticks([150, 175, 200,225, 250])
        
        # Cycle efficiency
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Cycle efficiency [\%]', color='black')   
        plt.plot(self.T_out_boiler_vect-273.15,self.eta_C_vect*100, color='darkorange', label = 'Carnot')
        plt.plot(self.T_out_boiler_vect-273.15,self.eta_vect*100, color=u'#1f77b4', label = 'Real')
        # plt.xticks([150, 175, 200,225, 250])
        plt.legend(loc='best')   
        
        # Condenser capacity
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Cond. capacity [kW]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.Q_dot_cd_vect/1000, color=u'#1f77b4')
        
        # Evaporation capacity
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Evap. capacity [kW]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.Q_dot_ev_vect/1000, color=u'#1f77b4')
        # plt.xticks([150, 175, 200,225, 250])
        
        # Pressure ratio
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Pressure ratio [-]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.r_p_vect, color=u'#1f77b4')
        # plt.xticks([150, 175, 200,225, 250])
        
        # Evap pressure
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Evap. pressure [bar]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.p_ev_vect/100000, color=u'#1f77b4')
        # plt.xticks([150, 175, 200,225, 250])

        # Condensation pressure
        plt.figure(figsize=(4.5/1.5, 3.7/1.5))
        plt.grid()
        plt.xlabel('Boiler outlet temperature [°C]')
        plt.ylabel('Cond.pressure [bar]', color='black')          
        plt.plot(self.T_out_boiler_vect-273.15, self.p_cd_vect/100000, color=u'#1f77b4')
        # plt.xticks([150, 175, 200,225, 250])
        
        

        
        
