# Constrain Curvature with Strong Lenses and Complimentary Probes
This repository contains all the data and code for the project "Constraining Curvature Density Parameter by Combining Time-Delay Lenses with Other Probes: a Forecast for Next-Generation Surveys". 

Users will need to install [astropy](https://www.astropy.org/), [emcee](https://emcee.readthedocs.io/en/stable/user/install/), [corner](https://corner.readthedocs.io/en/latest/install/) and [george](https://george.readthedocs.io/en/latest/user/quickstart/) python packages in order to run the code in the [Scripts](/Scripts) directory

In the [Scripts](/Scripts) directory, users can find code for analysis on constraining curvature density parameter and other cosmological parameters using either strong gravitational lenses alone or with other complimentary probes. The analysis are done in two ways: one is non-parametric, indicated by the presence of word **"HzGP"** in the names; second is model-dependent, with a mock universe of either oLCDM or owCDM model assumed.

In the [Data](/Data) directory, users can find simulated data files for cosmological probes with the following attributes:   
**LSST-like:**  
[zlens_zsource_LSSTLike.txt](/Data/zlens_zsource_LSSTLike.txt)   
Strong gravitational lenses, source: [LSST DESC data products](https://lsstdesc.org/)   
**Pantheon-like:**  
[sys_DS17f.csv](/Data/sys_DS17f.csv)     
[syscov_WFIRST.txt](/Data/syscov_WFIRST.txt)   
Type Ia supernovae, source: [Pantheon](https://github.com/dscolnic/Pantheon)   
**Roman-like:**    
[lcparam_DS17f.csv](/Data/lcparam_DS17f.csv)   
[lcparam_WFIRST_G10.txt](/Data/lcparam_WFIRST_G10.txt)   
Type Ia supernovae, source: Hounsell et al 2018   
**DESI-like:**                  
[DESI_HZ_error.txt](/Data/DESI_HZ_error.txt)   
Baryon accoustic oscillation, source: Font-Ribera et al 2016         

As the data used are mostly simulations for next-generation cosmological surveys, this analysis provides a forcast on the ability of using next-generation strong gravitational lenses data to constrain curvature density parameter. This is a promising way to obtain a high-fidelity constraint on curvature density parameter that is indepedent from early-universe probe CMB.
