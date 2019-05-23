# ExoplanetsSysSim
Welcome to the ExoplanetsSysSim package for generating planetary systems and simulating observations of those systems with a transit survey.  Currently, SysSim focuses on NASA's Kepler mission, but we've aimed to develop a generic framework that can be applied to other surveys (e.g., K2, TESS, PLATO, LSST, etc.).

# How to get started:
* Make sure you have julia v0.6 and git installed.
* If you want to use ssh keys instead of https authentication (to minimize typing your github password), then:
  * Setup a local ssh key using ssh-keygen
  * Tell Github about your ssh key:  Person Icon (upper right), Settings, SSH & GPG keys, New SSH Key.  Entry a name in the title box and paste the contents of `cat ~/.ssh/id_rsa.pub` into the "Key" box. Add SSH Key.  
* Run julia and install the ExoplanetsSysSim repo as a Julia package using the following command.
```
#!julia
   Pkg.clone("git@github.com:eford/ExoplanetsSysSim.jl.git")
```
* If you are using windows, you might encoutner an issue with capitalization of package names.
* If you have some issues with Blosc, you might need to follow the instructions here: https://github.com/stevengj/Blosc.jl/issues/5

It's recommended that you run the setup script below that will attempt to install two related packages (ABC, CORBITS) that are often used in combination with ExoplanetsSysSim.
```
#!julia
   include(joinpath(Pkg.dir("ExoplanetsSysSim"),"src/setup.jl"))  
```
* Optionally, run some tests, e.g. 
```
#!julia
   using ExoplanetsSysSim
   include(joinpath(Pkg.dir("ExoplanetsSysSim"),"examples/basic/test.jl"))   
```
* Change into the ExoplanetsSysSim directory (likely ~/.julia/v0.6/ExoplanetsSysSim or the similar appropriate path as the julia version increases).
* Run some simple "applications" (after exiting out of Julia and going to 
```
#!csh
cd apps
julia syssim_summary_stats.jl demo_param.in demo_ss.out
julia syssim_dist.jl demo_param.in demo_ss.out
```
* Create your own feature branch and start adding features to make SysSim even more realistic and powerful
* Write papers and cite relevant publications (e.g., Hsu et al. 2018)

NOTE: To run simulations which use the same forward model discussed in Hsu et al. (2018), please go to the "examples/hsu_etal_2018/" sub-directory.  
Alternatively, to run simulations which use the same forward model discussed in Hsu et al. (2019), please go to the "examples/dr25_gaia_fgk/" sub-directory.

# Team:
## Developers:
  * Eric Ford
  * Matthias He
  * Danley Hsu
  * Darin Ragozzine
## Other Contributors/Consultants:
  * Robert Morehead
  * Keir Ashby
  * Jessi Cisewski
  * Chad Schafer
  * Tom Loredo
  * Robert Wolpert

# Acknowledgements:
* NASA
  * Kepler Mission
  * Kepler Science Team
  * Kepler Multi-body & Transit Timing Variations Working Groups
  * Origins of Solar Systems program, award NNX14AI76G
  * Exoplanets Research Program, award NNX15AE21G
* The Pennsylvania State University
  * Dept. of Astronomy & Astrophysics
  * Center for Exoplanets & Habitable Worlds
  * Eberly College of Science
  * Institute for CyberScience
  * Center for Astrostatistics
  * Penn State Astrobiology Research Center
* Florida Institute of Technology
* University of Florida
* Statistical and Applied Mathematical Sciences Institute
