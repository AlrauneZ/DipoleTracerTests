#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function collection on analytical solutions for flow in a dipole flow setting 
in log-normal layered hydraulic conductivity media following the paper 

"Revisitation of the dipole tracer test for heterogeneous porous formations"
A. Zech, C. D’Angelo, S. Attinger, A. Fiori, J. Hydrol, 2018
https://doi.org/10.1016/j.advwatres.2018.03.006


@author: A. Zech
"""

import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.special import erf

DEF_settings = dict(
    por  = 0.3, # aquifer porosity
    dist = 10,  # distance between pumping and injection well in [m]
    q_dis = 1,  # unit dischange (injection/pumping rate)
    )

def time_to_dimensionless(time, # real time in [time-unit]
                          por  = 0.3, # aquifer porosity
                          dist = 10,  # distance between pumping and injection well in [m]
                          q_dis = 1.,  # unit dischange (injection/pumping rate) [m^3/time-unit]
                          depth = 1.,  # well/aquifer depth in [m]
                ):

    """ Formula for the dimensionless time from real time as function of 
    setting parameters
    
    INPUT
    -----
        t (ndarray)  - array with real time values 
        por (float)  - constant porosity value of the aquifer
        dist (float) - the distance between injection and pumping well
        q_dis (float) - unit discharge (injection/pumping rate) [m^3/time-unit] 
        depth (float) - depth of well (depth of aquifer in fully penetrating well setting) in [m]
    """
    tau = q_dis * time/(4.*np.pi*por*dist*dist*depth)
    
    return tau

def setup_time(tmin = 1e-13, 
               tmax = 1e3, 
               timesteps=50,
               dimensionless = True,
               ):

    """
    Function to setup time array for BTC analysis
        
    Input
    -----
       
        dimensionless (Boolean, default True):
            - if True: time is given in dimensionless form
            - if False: time array given in dimensional form

    """    

    # tau  - time in dimensionless form given as array in log-scale

    tau = np.logspace(np.log10(tmin), np.log10(tmax), timesteps)
    
    return tau
    

def btc_inst_flux(tau,
            var,
            ):
    
    """
    Breakthrough curve (BTC) of tracer, i.e. time dependent concentration C(t)
        - in a dipole flow setting 
        - for instantaneouse pulse injection
        - flux proportional injection mode
        - observed at the extraction well (at distance a)
        - in a layered aquifer with hydraulic conductivity K_i of layers
            is log-normal distributed
        - for a given dimensionless time (as function of settings)

    Functional implementation of equation 11 of Zech et al., 2018
 
    In flux proportional injection mode, the mass of solute entering each layer
    is proportional to the local velocity at the injection well, which is 
    variable in the vertical and is proportional to the hydraulic conductivity 
    Ki of each layer. This settings takes the local heterogeneity into account.    
    This is the typical injection condition in applications.
    
    Input
    -----
       tau (ndarray)  -  array with time values for BTC in dimensionless form
       var (float) - log-conductivity variance (i.e. variance of underlying normal distribution)
    
    
    Output
    -----
        pdf  (ndarray, length of tau) -  array with concentration values 
                                         [BTC corresponds to travel time pdf]  
    """

    pdf = np.zeros(len(tau))

    for i,ti in enumerate(tau):
        integrand = lambda theta : (1 / np.sin(theta)**2 * (1 - theta / np.tan(theta))) * np.exp(-((0.5*var + np.log(1 / np.sin(theta)**2 * (1 - theta / np.tan(theta))) - np.log(ti))**2 / (2*var)))
        pdf[i] = quad(integrand, 0, np.pi)[0]
    pdf = pdf / (np.pi * tau**2 * np.sqrt(2 * np.pi * var))    

    return pdf

def btc_const_flux_integral(tau,
            var,
            ):
    
    """
    Breakthrough curve (BTC) of tracer, i.e. time dependent concentration C(t)
        - in a dipole flow setting 
        - for constant injection
        - flux proportional injection mode
        - observed at the extraction well (at distance a)
        - in a layered aquifer with hydraulic conductivity K_i of layers
            is log-normal distributed
        - for a given dimensionless time (as function of settings)

    In flux proportional injection mode, the mass of solute entering each layer
    is proportional to the local velocity at the injection well, which is 
    variable in the vertical and is proportional to the hydraulic conductivity 
    Ki of each layer. This settings takes the local heterogeneity into account.    
    This is the typical injection condition in applications.

    Implementation of solution (BTC = cdf) as numerical cumulative integral 
    of pdf for instantaneous injection.
    
    Input
    -----
        tau (ndarray)  -  array with time values for BTC in dimensionless form
        var (float) - log-conductivity variance (i.e. variance of underlying normal distribution)
   
    
    Output
    -----
        tau_cdf (ndarray, length of tau -1)  -  array with adapted dim.-less time values 
        cdf  (ndarray, length of tau -1) - array with cdf values  
    """

    pdf = btc_inst_flux(tau,var)

    tau_cdf = 0.5*(tau[1:] + tau[:-1]) 
    cdf = cumulative_trapezoid(y = pdf,x = tau)

    return tau_cdf,cdf
    
def btc_const_flux(tau,
               var,
               ):
     
    """
    Breakthrough curve (BTC) of tracer, i.e. time dependent concentration C(t)
        - in a dipole flow setting 
        - for constant injection
        - flux proportional injection mode
        - observed at the extraction well (at distance a)
        - in a layered aquifer with hydraulic conductivity K_i of layers
            is log-normal distributed
        - for a given dimensionless time (as function of settings)

    Functional implementation of equation 15 of Zech et al., 2019 for constant
    input (without end of injection)
 
    In flux proportional injection mode, the mass of solute entering each layer
    is proportional to the local velocity at the injection well, which is 
    variable in the vertical and is proportional to the hydraulic conductivity 
    Ki of each layer. This settings takes the local heterogeneity into account.    
    This is the typical injection condition in applications.

    
    Input
    -----
        tau (ndarray)  -  array with time values for BTC in dimensionless form
        var (float) - log-conductivity variance (i.e. variance of underlying normal distribution)
   
    
    Output
    -----
        tau_cdf (ndarray, length of tau -1)  -  array with adapted dim.-less time values 
        cdf  (ndarray, length of tau -1) - array with cdf values  
    """

    pdf = np.zeros(len(tau))
    
    for i,ti in enumerate(tau):
        integrand = lambda theta : 1 + erf((0.5*var - np.log(1./np.sin(theta)**2 * (1-theta/np.tan(theta))) + np.log(ti))/(np.sqrt(2 * var)))
        pdf[i] = quad(integrand, 0, np.pi)[0]

    pdf = 0.5* pdf / np.pi     

    
    return pdf

def btc_step_flux(tau,
               var,
               step = 1,
               ):
     
    """
    Breakthrough curve (BTC) of tracer, i.e. time dependent concentration C(t)
        - in a dipole flow setting 
        - for constant injection of a duration "step
        - flux proportional injection mode
        - observed at the extraction well (at distance a)
        - in a layered aquifer with hydraulic conductivity K_i of layers
            is log-normal distributed
        - for a given dimensionless time (as function of settings)

    Functional implementation of equation 15 of Zech et al., 2019 for constant
    input (without end of injection)
 
    In flux proportional injection mode, the mass of solute entering each layer
    is proportional to the local velocity at the injection well, which is 
    variable in the vertical and is proportional to the hydraulic conductivity 
    Ki of each layer. This settings takes the local heterogeneity into account.    
    This is the typical injection condition in applications.

    Finite duration of constant injection is given in dimensionless time. If
    duration is given in real time "t_dur", it has to be converted using the 
    aquifer and setting parameters via:
        step  = Q*t_dur/(4*pi*por*dist*dist*depth)
        --> see also function for converting time to dimensionless time:
         time_to_dimensionless()
    
    Input
    -----
        tau (ndarray)  -  array with time values for BTC in dimensionless form
        var (float) - log-conductivity variance (i.e. variance of underlying normal distribution)
        step (float) - duraction of (constant) injection in dimensionless time
        
   
    
    Output
    -----
        tau_cdf (ndarray, length of tau -1)  -  array with adapted dim.-less time values 
        cdf  (ndarray, length of tau -1) - array with cdf values  
    """

    ### PHI (tau)
    pdf1 = np.zeros(len(tau))   
    for i,ti in enumerate(tau):
        integrand = lambda theta : 1 + erf((0.5*var - np.log(1./np.sin(theta)**2 * (1-theta/np.tan(theta))) + np.log(ti))/(np.sqrt(2 * var)))
        pdf1[i] = quad(integrand, 0, np.pi)[0]
    pdf1 = 0.5* pdf1 / np.pi     

    ###
    
    ### only calculate stop when input duration time is shorter than observation time
    if tau[-1]>step: 
        arg = np.argmin(abs(tau-step))
        print(tau[arg],step)
        pdf_step = np.zeros_like(tau)
        pdf_step[arg:] = pdf1[:-arg]

        pdf  = pdf1 - pdf_step
    else:
        pdf = pdf1
            
    return pdf

def btc_inst_resident(tau,
            var,
            ):
    
    """
    Breakthrough curve (BTC) of tracer, i.e. time dependent concentration C(t)
        - in a dipole flow setting 
        - for instantaneouse pulse injection
        - in resident injection mode
        - observed at the extraction well (at distance a)
        - in a layered aquifer with hydraulic conductivity K_i of layers
            is log-normal distributed
        - for a given dimensionless time (as function of settings)

    Functional implementation of equation 8 of Zech et al., 2018
    
    Resident injection mode assumes that the mass of solute entering in each 
    layer from the injection well is constant for all layers - despite the 
    different flow velocities in each layer given the heterogeneous K_i 
    distribution.
    
    Input
    -----
       tau (ndarray)  -  array with time values for BTC in dimensionless form
       var (float) - log-conductivity variance (i.e. variance of underlying normal distribution)
    
    
    Output
    -----
        pdf  (ndarray, length of tau) -  array with concentration values 
                                         [BTC corresponds to travel time pdf]  
    """

    pdf = np.zeros(len(tau))

    for i,ti in enumerate(tau):
        integrand = lambda theta : np.exp(- (0.5*var + np.log(1./np.sin(theta)**2 * (1 - theta / np.tan(theta))) - np.log(ti))**2 / (2.*var))
        pdf[i] = quad(integrand, 0, np.pi)[0]
    pdf = pdf / (np.pi * tau * np.sqrt(2 * np.pi * var))    

    return pdf
          
def btc_const_resident(tau,
                       var,
                       ):
    
    """
    Breakthrough curve (BTC) of tracer, i.e. time dependent concentration C(t)
        - in a dipole flow setting 
        - for constant injection
        - in resident injection mode
        - observed at the extraction well (at distance a)
        - in a layered aquifer with hydraulic conductivity K_i of layers
            is log-normal distributed
        - for a given dimensionless time (as function of settings)

    Functional implementation of equation 8 of Zech et al., 2018
    
    Input
    -----
        tau (ndarray)  -  array with time values for BTC in dimensionless form
        var (float) - log-conductivity variance (i.e. variance of underlying normal distribution)
   
    
    Output
    -----
        tau_cdf (ndarray, length of tau -1)  -  array with adapted dim.-less time values 
        cdf  (ndarray, length of tau -1) - array with cdf values  
    """

    pdf = btc_inst_resident(tau,var)

    tau_cdf = 0.5*(tau[1:] + tau[:-1]) 
    cdf = cumulative_trapezoid(y = pdf,x = tau)

    return tau_cdf,cdf
    

