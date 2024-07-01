[![DOI](https://zenodo.org/badge/645850928.svg)](https://zenodo.org/doi/10.5281/zenodo.12608026)

# Dipole Tracer Tests

Function collection on analytical solutions for flow in a dipole flow setting 
in log-normal layered hydraulic conductivity media following the paper 

"Revisitation of the dipole tracer test for heterogeneous porous formations"
A. Zech, C. D’Angelo, S. Attinger, A. Fiori, J. Hydrol, 2018
https://doi.org/10.1016/j.advwatres.2018.03.006

## Background and Overview 

Analytical solutions describe the breakthrough curve (BTC) of a tracer in 
dipole flow (also named two well flow) setting  where water is injected at one
well and extracted at another well. Tracer is added to the injection well either
as a instantaneous pulse (infinitisimaly short duration) or constantly. 

The solutions consider advective transport only. For dipole tests, it is well 
know that advection is the most significant source of spreading due to the 
non-uniform flow configuration. Local dispersion mechanisms like hydrodynamic 
dispersion or molecular diffusion are less significant and can be neglected.

The aquifer under consideration is non-homogeneous/heterogeneous where the structure
is layered with a random distribution of hydraulic conductivity values $K_i$ of each
layer follow a stochastic distribution. In the manuscript Zech et al., 2018, 
general solutions for the BTC are expressed as integrals. Specific solutions
are derived considering the spatial distribution of $K$ to be log-normal to a 
mean $K_G$ and a variance $\sigma_K^2$. These specific solutions are implemented here.

The BTC solutions present the dependence on both, the dipole setting and the 
aquifer configuration. The analysis is carried out by considering the travel 
time of a generic solute particle, from the injection to the pumping well. 
The probability density function (pdf) of such travel time is identical 
to the BTC of a solute instantaneous pulse. The corresponding cumulative 
density function (cdf's) describes the BTC of the solute for constant input. 
                                            
## Structure

Please try to organize your example in the given Structure
- `data/` - here you should place your input data
- `src/` - here you should place your python scripts
- `results/` - here your computed results and plots should be stored
- `README.md` - please describe your example in the readme, potentially showing results
- `LICENSE` - the default license is MIT, you can use another one if wanted


## Python environment

To make the example reproducible, it would be a good practice to provide one of
the following files:
- `requirements.txt` - requirements for [pip](https://pip.pypa.io/en/stable/user_guide/#requirements-files) to install all needed packages


## Contact

You can contact us via a.zech@uu.nl.


## License

MIT © 2024
