#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript reprocuding Figure 2 of Manuscript 

"Revisitation of the dipole tracer test for heterogeneous porous formations"
A. Zech, C. D’Angelo, S. Attinger, A. Fiori
https://doi.org/10.1016/j.advwatres.2018.03.006

@author: Alraune Zech
"""
import DipoleTracerTest as dtt
# import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

### ===========================================================================
### Key words to specify modus of script:
### ===========================================================================

textsize  = 12
var_list = [0.1, 1, 4, 8]
tmin  = -5
tmax  = 4

tau = dtt.setup_time(tmin = 10**tmin,tmax = 10**tmax, timesteps=500)

### ===========================================================================
### Prepape plot
### ===========================================================================

plt.figure(dpi = 300, figsize=(6,6))

for v in var_list:


    cdf = dtt.btc_const_flux(tau,v)
    plt.plot(tau, cdf,label = v)

    tau_cdf,cdf2 = dtt.btc_const_flux_integral(tau,v)
    plt.plot(tau_cdf, cdf2,ls = '--',c = 'k')

plt.axis([1e-5,1e4, 0,1])
plt.xscale("log")
plt.grid(True)
plt.tick_params(labelsize = textsize)
plt.legend(title = "$\sigma^2=$", fontsize = textsize, title_fontsize = textsize)
plt.xlabel(r"$\tau$", fontsize = textsize)
plt.ylabel("$C/C_0$", fontsize = textsize)
# plt.savefig("../results/Fig02_cdf.pdf")
