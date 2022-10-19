"""
Author: Yang Hu
1. This file gives a script to constrain cosmological parameters (H0, Omega0, Omegak, MB, etc.) in a given model 
using 3 probes: strong gravitational lenses (lenses), type Ia supernovae (SNe) and baryon accoustic oscillation
(BAO).
2. For lenses, we use simulated LSST data. For SNe, we either use simulated data of Pantheon (current generation 
survey) or Roman (next-generation survey). For BAO, we use simulated DESI data.
3. Constraints on parameters is obtained via Markov Chain Monte Carlo (MCMC) using emcee package.
4. A model of mock universe has to be chosen before doing the analysis.
"""


"""
standard imports for data analysis and astropy.cosmology to compute astrophysical quantities with ease
"""

import numpy as np
import pandas as pd
import emcee
import corner
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from astropy.cosmology import FlatLambdaCDM, FlatwCDM, LambdaCDM, wCDM



"""
Load original data. Please check that the relevant data files are named and stored in the way the code below 
specified.
"""

#load&rename lenses data
data1 = pd.read_csv("Data/zlens_zsource_LSSTLike.txt", delimiter=' ', header=None)
zlens = data1[0]
zsource = data1[1]
ddt_err = data1[2]

#load&rename SNe data
#choose either Pantheon or Roman for SNe
Pantheon = False
Roman = True
if Pantheon:
    data2 = pd.read_csv("Data/lcparam_DS17f.csv", delimiter=',', skiprows=1, header=None)
    data3 = pd.read_csv("Data/sys_DS17f.csv", delimiter=',', skiprows=1, header=None)  
    name = "Pantheon"
elif Roman:
    data2 = pd.read_csv("Data/lcparam_WFIRST_G10.txt", delimiter=' ', skiprows=1, header=None)
    data3 = pd.read_csv("Data/syscov_WFIRST.txt", delimiter=' ', skiprows=0, header=None)
    name = "Roman"
zcmb = data2[1]
#for a measure of the uncertainty of SNe measurement, we need the covariance matrix in the following way
mb_err = data2[5]
sys_err = np.array(data3)
D_stat = np.diag(mb_err**2)
C_sys = sys_err.reshape((len(data2), len(data2)))
C = D_stat + C_sys
C_inv = np.linalg.inv(C)

#load&rename BAO data
data4 = pd.read_csv("Data/DESI_HZ_error.txt", delimiter=' ', skiprows=1, header=None) 
zBAO = data4[0]
sigHz = data4[1]



"""
Sometimes we want more simulated data for lenses. Here we provide a way to genearate more LSST-like data.
The number of lensing events in the original LSST-like data is 310, an estimate of expected observed number
of events in LSST's 10-year survey baseline.
"""

#decide whether to generate data and uncertainty for lenses and whether to use them
generate_data = False
generate_uncertainty = False
save_data = False
use_data = True

#first entry for number is the total number of data points we want, second entry is the number of original data
number = 3000-310
#choose a python random seed for random processes involved
seed_no = 22

if generate_uncertainty:
    np.random.seed(seed_no)
    pu = np.random.uniform(0.06, 0.1, size=len(data1[0]))

if generate_data:
    np.random.seed(seed_no)
    rand_number = np.random.randint(0, 309, size=number)

if save_data:
    data_temp = np.array(data1)
    data_temp[:, 2]=pu
    for i in rand_number:
        data_temp = np.append(data_temp, [data_temp[i]], axis=0)
    df = pd.DataFrame(np.concatenate(([data_temp[:, 0]], [data_temp[:, 1]], [data_temp[:, 2]])).T)
    df.to_csv(r'Data/zlens_zsource_%sLSSTLike_%s.csv' % ((number+310), seed_no), index=False)

if use_data:
    data_new = pd.read_csv("Data/zlens_zsource_%sLSSTLike_%s.csv" % ((number+310), seed_no), skiprows=1, header=None)
    zlens = data_new[0]
    zsource = data_new[1]
    ddt_err = data_new[2]
else:
    zlens = data1[0]
    zsource = data1[1]
    ddt_err = data1[2]



"""
We need to choose a mock universe model to generate the "true value" of parameters which are then
used to compute the difference between true and measured (simulated) values which in turn are used to get
the constraints of parameters.
"""

#first two plots a distribution of the lenses data used
#third decides whether to add in magnitude scattering for SNe data
hist_plot = False
scatter_plot = False
magnitude_scatter = True

#set cosmology here
#models available by default are "oLCDM" and "owCDM"
cosmology = "owCDM"

#set mock cosmology
#common parameters we are interested in are Hubble's constant H0, matter density Om0, curvature density Ok0,
#equation of state parameter w and absolute magnitude of SNe Ia MB.
if cosmology == "oLCDM":
    H0_mock, Om0_mock, Ok0_mock, MB_mock = 72, 0.3, 0.0, -19.2 
    cosmo_mock = LambdaCDM(H0=H0_mock, Om0=Om0_mock, Ode0=1.-Om0_mock-Ok0_mock)
elif cosmology == 'owCDM':
    H0_mock, Om0_mock, Ok0_mock, w_mock, MB_mock = 72, 0.3, 0.00, -1, -19.2 
    cosmo_mock = wCDM(H0=H0_mock, Om0=Om0_mock, Ode0=1.-Om0_mock-Ok0_mock, w0=w_mock)

#compute mock time-delay distance for lenses
dd_mock = cosmo_mock.angular_diameter_distance(z=zlens).value
ds_mock = cosmo_mock.angular_diameter_distance(z=zsource).value
dds_mock = cosmo_mock.angular_diameter_distance_z1z2(z1=zlens, z2=zsource).value
ddt_mock = (1. + zlens) * dd_mock * ds_mock / dds_mock

#compute mock luminosity distance for SNe
dl_mock = cosmo_mock.luminosity_distance(z=zcmb).value

#compute mock Hz from BAO data
Hz_mock = cosmo_mock.H(z=zBAO).value

#add in magnitude scattering of absolute magnitude of SNe
if magnitude_scatter:
    np.random.seed(seed_no)
    if Pantheon:
        mb_mock = 5*np.log10(dl_mock)+25+MB_mock+np.random.normal(0, 0.02, size=len(zcmb))
    elif Roman: 
        mb_mock = 5*np.log10(dl_mock)+25+MB_mock+np.random.normal(0, 0.01, size=len(zcmb))
else:
    mb_mock = 5*np.log10(dl_mock)+25+MB_mock


#relevant test plots
if hist_plot:
    fig, ax = plt.subplots(1, figsize=(10, 10))
    ax.hist(zlens, bins=np.linspace(0.1, 0.8, 8))
    ax.set_title("Histogram of zlens, %s LSST-like samples" % len(zlens))
    ax.set_xlabel("zlens")
    fig.savefig("Plots/zlens_%sLSSTLike_samples_hist.png" % len(zlens))
if scatter_plot:
    fig, ax = plt.subplots(1, figsize=(10, 10))
    ax.scatter(zlens, ddt_mock)
    ax.set_title("Scatter plot of time-delay distances VS zlens, %s LSST-like samples, mock=%s,%s,%s,%s,%s" 
                 % (len(zlens), H0_mock, Om0_mock, Ok0_mock, MB_mock, cosmology)
                   )
    ax.set_xlabel("zlens")
    ax.set_ylabel("time-delay distances [Mpc]")
    fig.savefig("Plots/ddtVSzlens_%sLSSTLike_samples_%s,%s,%s,%s_%s_scatter.png"  
                 % (len(zlens), H0_mock, Om0_mock, Ok0_mock, MB_mock, cosmology)
               )


"""
Constraints are obtained by MCMC and here we define the relevant prior and likelihood functions 
"""
#decide whether to run a maximum likelihood testing before running the full lengthy MCMC analysis
ml_test = True

#use uniform priors for all parameters
def log_prior(theta):
    if cosmology == "oLCDM":
        h0, om, ok, Mb = theta
        if 0. <= h0 <= 150. and 0.05 <= om <= 0.5 and -0.5 <= ok <= 0.5 and -38.4 <= Mb <= 0.:
            return 0.0
        else:
            return -np.inf
    elif cosmology == "owCDM":
        h0, om, ok, w, Mb = theta
        if 0. <= h0 <= 150. and 0.05 <= om <= 0.5 and -0.5 <= ok <= 0.5 and -2.5 <= w <= 0.5 and -38.4 <= Mb <= 0.:
            return 0.0
        else:
            return -np.inf
    else:
        raise ValueError("I don't know the cosmology %s" % cosmology)

#use a chi-square likelihood function
def log_likelihood(theta, zlens, zsource, ddt_err, zcmb, C_inv, zBAO, sigHz): 
    if cosmology == "oLCDM":
        h0, om, ok, Mb = theta
        #check parameters are physical
        if (om < 0 or om > 1 or 1.-om-ok < 0 #or 1.-om-ok > 1.
            or np.any(ok*(1.0+zsource)**2 + om*(1.0+zsource)**3 + (1.0-om-ok) <= 0)
            or np.any(ok*(1.0+zcmb)**2 + om*(1.0+zcmb)**3 + (1.0-om-ok) <= 0)
           ):
            return -np.inf
        cosmo = LambdaCDM(H0=h0, Om0=om, Ode0=1.0-om-ok)
    elif cosmology == "owCDM":
        h0, om, ok, w, Mb = theta
        #check parameters are physical
        if (om < 0 or om > 1 or 1.-om-ok < 0 #or 1.-om-ok > 1.
            or np.any(ok*(1.0+zsource)**2 + om*(1.0+zsource)**3 + (1.0-om-ok)*(1.0+zsource)**(3*(1+w)) <= 0)
            or np.any(ok*(1.0+zcmb)**2 + om*(1.0+zcmb)**3 + (1.0-om-ok)*(1.0+zcmb)**(3*(1+w)) <= 0)
           ):
            return -np.inf
        cosmo = wCDM(H0=h0, Om0=om, Ode0=1.0-om-ok, w0=w)
    else:
        raise ValueError("I don't know the cosmology %s" % cosmology)
    dd = cosmo.angular_diameter_distance(z=zlens).value
    ds = cosmo.angular_diameter_distance(z=zsource).value
    dds = cosmo.angular_diameter_distance_z1z2(z1=zlens, z2=zsource).value
    ddt = (1. + zlens) * dd * ds / dds
    chi_sq_lenses = np.sum((ddt-ddt_mock)**2./(ddt*ddt_err)**2.)#compute chi_square for lenses
    dl = cosmo.luminosity_distance(z=zcmb).value
    mb_model = 5*np.log10(dl)+25+Mb
    del_m = mb_mock - mb_model
    chi_sq_SNe = np.dot(del_m.T, np.dot(C_inv, del_m))#compute chi_square for SNe
    Hz = cosmo.H(z=zBAO).value
    chi_sq_BAO = np.sum((Hz-Hz_mock)**2./sigHz**2.)#compute chi_square for BAO 
    return -0.5*(chi_sq_lenses+chi_sq_SNe+chi_sq_BAO)

def log_probability(theta, zlens, zsource, ddt_err, zcmb, C_inv, zBAO, sigHz):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, zlens, zsource, ddt_err, zcmb, C_inv, zBAO, sigHz)

#run a maximum likelihood test
if ml_test:
    nll = lambda *args: -log_likelihood(*args)
    if cosmology == "oLCDM":
        initial = np.array([70., 0.27, 0.02, -19.])
        soln = minimize(nll, initial, args=(zlens, zsource, ddt_err, zcmb, C_inv, zBAO, sigHz))
        H0_ml, Om_ml, Ok_ml, Mb_ml = soln.x
        print("Maximum likelihood estimates:")
        print("H0_ml = {0:.3f}".format(H0_ml))
        print("Om_ml = {0:.3f}".format(Om_ml))
        print("Ok_ml = {0:.3f}".format(Ok_ml))
        print("MB_ml = {0:.3f}".format(Mb_ml))
    elif cosmology == "owCDM":
        initial = np.array([70., 0.27, 0.02, -1.1, -19.])
        soln = minimize(nll, initial, args=(zlens, zsource, ddt_err, zcmb, C_inv, zBAO, sigHz))
        H0_ml, Om_ml, Ok_ml, w_ml, Mb_ml = soln.x
        print("Maximum likelihood estimates:")
        print("H0_ml = {0:.3f}".format(H0_ml))
        print("Om_ml = {0:.3f}".format(Om_ml))
        print("Ok_ml = {0:.3f}".format(Ok_ml))
        print("w_ml = {0:.3f}".format(w_ml))
        print("MB_ml = {0:.3f}".format(Mb_ml))


"""
presettings for MCMC
"""

nwalkers = 32
nsamples = 20000

if cosmology == "oLCDM":
    startpos = [70., 0.27, 0.02, -19.]  # H0, Om, Ok, MB
    labels = ["H0", "Om", "Ok", "MB"]
if cosmology == "owCDM":
    startpos = [70., 0.27, 0.02, -1.1, -19.]  # H0, Om, Ok, MB
    labels = ["H0", "Om", "Ok", "w", "MB"]



"""
run MCMC using emcee package with the presetting above
"""

#decide whether to run MCMC and save the obtained samples, and whether to load a certain sample file
run_mcmc = False
save_mcmc = False
load_mcmc = True

if run_mcmc:
    pos = startpos + 1e-4 * np.random.randn(nwalkers, len(startpos))
    nwalkers, ndim = pos.shape
    display("Sampling cosmological parameters in %s..." % cosmology)
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability, args=[zlens, zsource, ddt_err, zcmb, C_inv, zBAO, sigHz]
    )
    sampler.run_mcmc(pos, nsamples, progress=True);

if save_mcmc:
    samples = sampler.get_chain()
    flat_samples = sampler.get_chain(discard=5000, thin=4, flat=True) #burn-in
    sample_data = pd.DataFrame(flat_samples)
    if cosmology == "oLCDM":
        sample_data.to_csv(
            "Samples/Simulation/%sLSST+%s+DESI_mock=%s,%s,%s,%s,%s_%ix%i.csv" 
            % (len(zlens), name, H0_mock, Om0_mock, Ok0_mock, MB_mock, cosmology, nwalkers, nsamples), 
            index=False, header=labels
        )
    elif cosmology == "owCDM":
        sample_data.to_csv(
            "Samples/Simulation/%sLSST+%s+DESI_mock=%s,%s,%s,%s,%s,%s_%ix%i.csv" 
            % (len(zlens), name, H0_mock, Om0_mock, Ok0_mock, w_mock, MB_mock, cosmology, nwalkers, nsamples), 
            index=False, header=labels
        )


if load_mcmc:
    if cosmology == "oLCDM":
        flat_samples = pd.read_csv("Samples/Simulation/%sLSST+%s+DESI_mock=%s,%s,%s,%s,%s_%ix%i.csv" 
                                   % (len(zlens), name, H0_mock, Om0_mock, Ok0_mock, MB_mock, 
                                      cosmology, nwalkers, nsamples), skiprows=1, header=None
                                  )
    if cosmology == "owCDM":
        flat_samples = pd.read_csv("Samples/Simulation/%sLSST+%s+DESI_mock=%s,%s,%s,%s,%s,%s_%ix%i.csv" 
                                   % (len(zlens), name, H0_mock, Om0_mock, Ok0_mock, w_mock, MB_mock, 
                                      cosmology, nwalkers, nsamples), skiprows=1, header=None
                                  )


"""
After loading a sample file, we can plot a corner plot showing how well constraints are placed on each interested
cosmological paramters by the data.
"""

#corner plot
if cosmology == "oLCDM":
    fig = corner.corner(
        flat_samples, labels=labels, quantiles=(0.16, 0.84), show_titles=True, 
        title_fmt='.3f', use_math_text=True, truths=[H0_mock, Om0_mock, Ok0_mock, MB_mock]
    )
    fig.suptitle('Corner plot, %s LSST + %s + DESI samples, mock=%s,%s,%s,%s,%s' 
        % (len(zlens), name, H0_mock, Om0_mock, Ok0_mock, MB_mock, cosmology), y=1.02)
    fig.savefig("Plots/%sLSST+%s+DESI_mock=%s,%s,%s,%s,%s_%ix%i_corner.png" 
                % (len(zlens), name, H0_mock, Om0_mock, Ok0_mock, MB_mock, cosmology, nwalkers, nsamples)
                , bbox_inches='tight'
               )
elif cosmology == "owCDM":
    fig = corner.corner(
        flat_samples, labels=labels, quantiles=(0.16, 0.84), show_titles=True, 
        title_fmt='.3f', use_math_text=True, truths=[H0_mock, Om0_mock, Ok0_mock, w_mock, MB_mock]
    )
    fig.suptitle('Corner plot, %s LSST + %s + DESI samples, mock=%s,%s,%s,%s,%s,%s' 
        % (len(zlens), name, H0_mock, Om0_mock, Ok0_mock, w_mock, MB_mock, cosmology), y=1.02)
    fig.savefig("Plots/%sLSST+%s+DESI_mock=%s,%s,%s,%s,%s,%s_%ix%i_corner.png" 
                % (len(zlens), name, H0_mock, Om0_mock, Ok0_mock, w_mock, MB_mock, cosmology, nwalkers, nsamples)
                , bbox_inches='tight'
               )

