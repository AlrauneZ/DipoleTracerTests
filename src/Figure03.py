#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript reprocuding Figure 3 of Manuscript 

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
tmin  = -13
tmax  = 3

tau = dtt.setup_time(tmin = 10**tmin,tmax = 10**tmax, timesteps=500)

plt.figure(dpi = 300, figsize=(6,6))
for v in var_list:
    pdf = dtt.btc_inst_flux(tau,v)
    plt.plot(tau, pdf,label = v)

plt.axis([1e-13,1e3, 1e-4,1e3])
plt.xscale("log")
plt.yscale("log")
plt.grid(True)
plt.legend(title = r'$\sigma^2 = $',fontsize = textsize)
plt.xlabel(r"$\tau$",fontsize = textsize)
plt.ylabel("PDF",fontsize = textsize)
# plt.savefig("../results/Fig03_pdf.pdf")

