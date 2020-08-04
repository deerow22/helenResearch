from astropy.timeseries import BoxLeastSquares
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import scipy
import glob
import os
import pandas as pd
import lightkurve as lk

######################## create simulated light curves ########################

def one_transit(t=np.linspace(0,27,19440), 
           per=1., rp=0.1, t0=1., a=15., inc=87., ecc=0., 
           w=90., limb_dark ='nonlinear', u=[0.5,0.1,0.1,-0.1]):
        
    """
    ~Simulates a one-sector long TESS light curve with injected planet transits per input parameters.~
    Requires: batman; numpy
    Args:     t          =times at which to calculate light curve, default is one TESS sector;
              per        =orbital period;
              rp         =planet radius (in units of stellar radii);
              t0         =time of inferior conjunction);
              a          =semi-major axis (in units of stellar radii);
              inc        =orbital inclination (in degrees);
              ecc        =eccentricity;
              w          =longitude of periastron (in degrees);
              limb_dark  =limb darkening model;
              u          =limb darkening coefficients [u1, u2, u3, u4];
    
    outputs: flux array  =light curve with one injected transit at per, for use right before sim_lc to get TESS lc
    
    """
    #### maybe should make params its own fcn and split this fcn into 2....
    import batman
    params = batman.TransitParams(); params.t0 = t0; params.per = per                       
    params.rp = rp; params.a = a; params.inc = inc; params.ecc = ecc                       
    params.w = w; params.limb_dark = limb_dark; params.u = u      

    m = batman.TransitModel(params, t)    #initializes model
    flux = m.light_curve(params)        #calculates light curve
    
    return flux, m, params #flux,model,params

def sim_lc(flux,time,noise_level,plot=False, transit = None):
    '''
    ~takes flux/time arrays with noise and creates a light curve with 2min cadence for one TESS sector~
    Requires: numpy; lightkurve
    Args:  flux           -(array) flux values
           time           -(array) time values
           noise_level    -(array or int) flux uncertainity
    Returns:  LightKurve lightcurve object
    '''
    unc_level=noise_level
    flux_wnoise = np.random.randn(len(flux))*unc_level + flux
    noise = np.ones_like(flux)*unc_level
    #make it a lk object for phase-folding plots below
    lc_wnoise = lk.lightcurve.TessLightCurve(flux=flux_wnoise,time=time,flux_err=noise)
    #plot verification
    if plot ==True:
        lc_wnoise.plot()
        if transit != None:
            plt.title('Simulated Batman Light Curve with noise & orb.per: {}'.format(transit))
        else:
            plt.title('Simulated Batman Light Curve with noise')
    return lc_wnoise

def fake_data(N=1000,unc=1e-3,time_start=0,time_end=25):
    '''
    ~Generate fake data for a lightcurve~
    Requires: numpy;
    Args:       N           -(int)number of desired data points
                unc         -(int|float or array)level of uncertainity for flux
                time_start  -(int)starting time
                time_end    -(int)stopping time
    Returns:    time,flux,flux_err arrays
    '''
    sample_size = N
    time = np.linspace(time_start,time_end,sample_size)
    unc = unc
    flux = np.random.randn(len(time))*unc + 1
    flux_err = np.ones_like(flux)*unc
    
    return time,flux,flux_err

########################get real data########################

# open only sector1 lcs
def open_mylcs(tic,sector):
    '''
    ~Opens previously cleaned light curve file~
    Requires:       astropy.io.fits;
    Args:           tic         -(int) TESS TIC ID
                    sector      -(int) light curve sector number for target
    Returns:        time, flux, flux_err (numpy arrays) 
    '''
    filepath = 'data/SECONDRUN/cleaned_LightCurves/{}/sector{}_lc.fits'.format(tic,sector)
    try:
        lc = fits.open(filepath)
        lc_data = lc[1].data
        time=lc_data.TIME; flux=lc_data.FLUX; flux_err=lc_data.FLUX_ERR
    except Exception as e:
        print('TIC: {} Sector: {} couldnt be opened; encountered Error: {}'.format(tic,sector,e))
    
    return time,flux,flux_err

def kepEBopen(kic,path='mypath'):
    '''
    ~Opens Kepler EB data files I downloaded to my computer~
    Args: kic       -(int or str) id for kepler target
          path      -(str) path to csv data file, default is mypath, 
                      change it with {} where kic goes, header/delimiter 
                      specified for opening with pd.read_csv
    Returns: lightkurve Light Curve data object
    '''
    if path == 'mypath':
        try:
            kep_lcfile = pd.read_csv('data/KeplerEBs/{}_lc.csv'.format(kic),header=0,delimiter ='	')
            ktime = kep_lcfile['# bjd']
            kflux = kep_lcfile['dtr_flux']
            kfluxerr = kep_lcfile['dtr_err']
            klc = lk.lightcurve.LightCurve(flux=kflux,time=ktime,flux_err=kfluxerr)#make lk object
        except Exception as e:
            klc = "Verify KIC - couldn't open data file with mypath"
    else:
        try:
            kep_lcfile = pd.read_csv(path)
            ktime = kep_lcfile['# bjd']
            kflux = kep_lcfile['dtr_flux']
            kfluxerr = kep_lcfile['dtr_err']
            klc = lk.lightcurve.LightCurve(flux=kflux,time=ktime,flux_err=kfluxerr) #make lk object
        except Exception as e:
            klc = "Verify file path - couldn't open data file"
    return klc
    
######################## useful tools #########################
def period_grid(start=0.3,stop=27,N=25000):
    '''
    ~generates list of periods for BLS to check~
    Requires: Numpy;
    Args:          start    -(int)lowest period to check;
                   stop     -(int)highest period to check;
                   N        -(int)number of periods between interval
                            Default values optimized for single TESS sector
    Returns:       array of periods to check
    '''
    p_grid = np.linspace(start,stop,N)
    return p_grid

def duration_grid(start=.01,stop=0.29,N=5):
    '''
    ~generates list of durations for BLS to check~
    Requires: Numpy;
    Args:          start  -(int)lowest period to check;
                   stop   -(int)highest period to check;
                   N      -(int)number of periods between interval
                          Default values optimized for single TESS sector
    Returns:       array of durations to check
    '''
    dur_grid = np.linspace(start,stop,N)
    return dur_grid




######################## do things to data ########################
def tdepth(flux):
    avg = np.mean(flux) ### avg of full lc or flux with transits masked??
    dip = np.min(flux)
    depth = avg-dip
    return depth

def plot_lc(time,flux,flux_err='None',ID = 'Light Curve - no id'):
    '''
    ~Plots light curve~
    Requires:        matplotlib.pyplot; 
    Args:            time         -(np.array)
                     flux         -(np.array) 
                     flux_err     -Optional(np.array)-default is none
                     ID           -Optional(str or int) target id number for plot title;
    Returns:         Nothing, it just automatically plots the light curve
                     
    '''
    
    if flux_err != 'None':
        plt.figure(figsize=(20,10))
        plt.errorbar(time,flux,yerr=flux_err,color='lightgray')
        plt.scatter(time,flux,s=2,color='k')
        plt.xlabel('Time',fontsize=25)
        plt.ylabel('Flux',fontsize=25)
        plt.xticks(fontsize=20);plt.yticks(fontsize=20)
        plt.title('{}'.format(ID),fontsize=30)
        plt.show()
        # f = plt.gcf()
    else:
        plt.figure(figsize=(20,10))
        plt.scatter(time,flux,s=2,color='k')
        plt.xlabel('Time',fontsize=25)
        plt.ylabel('Flux',fontsize=25)
        plt.xticks(fontsize=20);plt.yticks(fontsize=20)
        plt.title('{}'.format(ID),fontsize=30)
        plt.show()
        # f = plt.gcf()
        #pass
    
    return plt.gcf()

def flatten_lc(time,flux,flux_err):
    '''
    ~Puts time,flux,fluxerr into a LightKurve Object to flatten - lk uses a wrapper
    on scipy.signal.savgol_filter.~
    
    Requires:            LightKurve;
    Args:                time      -(np.array)
                         flux      -(np.array)
                         flux_err  -(np.array)
    Returns:             Flattened time/flux/fluxerr arrays
    '''
    ##### change assert to if --return error stmt
    assert(len(time)==len(flux)==len(flux_err))#check same lenght
    lc = lk.lightcurve.LightCurve(flux=flux,time=time,flux_err=flux_err)#take generic lc
    flat_lc = lc.flatten()
    flat_time = flat_lc.time - flat_lc.time[0]
    flat_flux = flat_lc.flux
    flat_flux_err = flat_lc.flux_err
    
    return flat_time, flat_flux, flat_flux_err

############################ BLS #######################################

def bls(period_grid,duration_grid,time,flux,flux_err=0.):
    '''
    ~Runs Box Least Squares (BLS) on lightcurve data~
    Requires:        astropy.timeseries.BoxLeastSquares;
    Args:          period_grid     -(array) periods to check
                   duration_grid   -(array) durations to check
                   time            -(array) lightcurve times
                   flux            -(array) lightcurve flux
                   flux_err        -Optional(array or int) lightcurve flux error
    Returns:       BLS stats (periodogram, [power, period, depth, transit_time, duration])
    '''
    if(len(time)!=len(flux)):
        return 'ERROR: time and flux arrays NOT same length'
#     assert max(duration_grid)<min(period_grid), 'Minimum Period must be Greater than Maximum Duration'
    elif max(duration_grid)>min(period_grid):
        return 'ERROR: Minimum Period must be Greater than Maximum Duration'
    else:
        model = BoxLeastSquares(time, flux, dy=flux_err)
        periodogram = model.power(period_grid, duration_grid) 
        ppower = np.argmax(periodogram.power)
        pperiod = periodogram.period[np.argmax(periodogram.power)]
        pdepth = periodogram.depth[np.argmax(periodogram.power)]
        pduration = periodogram.duration[np.argmax(periodogram.power)]
        ptransittime = periodogram.transit_time[np.argmax(periodogram.power)]
#       #to use function attributes instead-----DONT USE --doesnt work well with loops
#         bls.power = np.argmax(periodogram.power)
#         bls.period = periodogram.period[np.argmax(periodogram.power)]
#         bls.depth = periodogram.depth[np.argmax(periodogram.power)]
#         bls.duration = periodogram.duration[np.argmax(periodogram.power)]
#         bls.transittime = periodogram.transit_time[np.argmax(periodogram.power)]
    
        return periodogram, [ppower, pperiod, pdepth, pduration, ptransittime] #bls
