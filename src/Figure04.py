#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript reprocuding Figure 4 of Manuscript 

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
var_list = [0.1, 1]
colors = ['C3','C0']
tmin  = -3
tmax  = 1

tau = dtt.setup_time(tmin = 10**tmin,tmax = 10**tmax, timesteps=500)

plt.figure(dpi = 300, figsize=(6,6))

for i,v in enumerate(var_list):
    pdf = dtt.btc_inst_flux(tau,v)
    plt.plot(tau, pdf,c = colors[i],label = r"flux, $\sigma^2=${:.1f}".format(v))

for i,v in enumerate(var_list):
    pdf = dtt.btc_inst_resident(tau,v)
    plt.plot(tau, pdf,c = colors[i],ls = '--',label = r"resident, $\sigma^2=${:.1f}".format(v))

plt.axis([1e-3,1e1, 0,1.4])
plt.xscale("log")
plt.ylabel("PDF",fontsize = textsize)
plt.xlabel(r"$\tau$",fontsize = textsize)
plt.grid(True)
plt.legend(fontsize = textsize)
# plt.savefig("../results/Fig04_Res_Flux.pdf")