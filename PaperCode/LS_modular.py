
import glob
from astropy.io import fits
import lightkurve as lk
import numpy as np
import starspot as ss

####---keeping below so doesnt break things but replacing with more modular below
def open_mylcs(tic,sector,path):
    '''
    ~Opens previously cleaned (for LS) light curve file--*might* have to be LightKurve object~
    Requires:       astropy.io.fits;
    Args:           tic         -(int) TESS TIC ID
                    sector      -(int) light curve sector number for target
                    path        -(str) path to light curve file; 
                                        Options for D.R. only:'mycvzpath' for cvz sectors,
                                        'mysecpath' for sectors 14 & 15, 'mystitchedcvzpath'
                                        with sector=any integer for stitched cvz lcs
    Returns:        time, flux, flux_err (numpy arrays) 
    '''

    if path == 'mycvzpath':
        filepath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/CleanLCs/CLEAN_CVZs/{}/sec{}_lc.fits'.format(tic,sector) #cvz targets' sectors
    elif path =='mysecpath':
        filepath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/CleanLCs/CLEAN_1415/{}/sec{}_lc.fits'.format(tic,sector) #targets in sectors 14/15
    elif path == 'mystitchedcvzpath':
        filepath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/CleanLCs/CLEAN_CVZs/{}/stitched_lc.fits'.format(tic) #cvz targets' sectors
    else:
        filepath = path
    try:
        lc = fits.open(filepath) #open file
        lc_data = lc[1].data #grab data arrays from file
        time=lc_data.TIME; flux=lc_data.FLUX; flux_err=lc_data.FLUX_ERR #isolate data arrays
    except Exception as e:
        print('TIC: {} Sector: {} couldnt be opened; encountered Error: {}'.format(tic,sector,e))
        time,flux,flux_err = 'None','None','None'
    
    return time,flux,flux_err
####---keeping above so doesnt break things but replacing with more modular below



def pathfinder(tic,sector):
    '''
    ~Determines correct path to cleaned data in PAPER_FINAL_FILES for D.R. 
    Others should amend paths used internally or bypass this function and 
    only use open_mylcs() with their path as input~
    REQUIRES: numpy as np
    ARGS:
        tic      -(int) TICID of target
        sector   -(str or int) 
                      'stitched' for stitched cvz stars (N or S)
                      or integer of single observation sector 
    RETURNS:
        path
    '''

    err_stmt = 'Sector type not understood. Specify observation sector integer or use "stitched" for stitched sectors'
    if type(sector) == type('string'):
        if sector == 'stitched':
            path = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/CleanLCs/CLEAN_CVZs/{}/stitched_lc.fits'.format(tic) #path to cvz stitched lcs
        else:
            print(err_stmt)
            path = 'None'
    elif type(sector) ==int: 
        if sector == 14 | sector ==15:
            path = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/CleanLCs/CLEAN_1415/{}/sec{}_lc.fits'.format(tic,sector) #path to sector 14 & 15 lcs
        else:
            path = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/CleanLCs/CLEAN_CVZs/{}/sec{}_lc.fits'.format(tic,sector) #path to cvz sector lcs
    else:
        print(err_stmt)
        path = 'None'
    return path

def pathfinder2(tic,stitched=False): #finds path to lcs cleaned for BLS
    '''
    ~Finds path to cleaned data in PAPER_FINAL_FILES~
    REQUIRES: numpy as np
    ARGS:
        tic      -(int) TICID of target
        stitched   -(boolean, default==False) 
                      True for stitched light curve only
                      or False for indiviual sector light curves
    RETURNS:
        list of available paths
    '''
    err_stmt = 'Sector type not understood. Specify observation sector(s) integer(s) in a list or use "stitched" for stitched sectors'
    if stitched ==True:
        fullpath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/CleanLCs/Clean_LCs_BLS/{}/stitched_lc.fits'.format(tic)
        path = glob.glob(fullpath)
    elif stitched ==False: 
        fullpath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/CleanLCs/Clean_LCs_BLS/{}/sec*_lc.fits'.format(tic)
        path = glob.glob(fullpath) #, recursive = True)
    else:
        print(err_stmt)
        path = 'None'
    return path





#more modular version of open_mylcs() & used in run_LS.py
#use these three functions instead of open_mylcs() above if need to examine file header first
def open_lcf(path):
    '''
    ~Tries to open cleaned light curve files~
    REQUIRES: lightkurve as lk
    ARGS:
        path   -(str) returned string from pathfinder function, 
                        or a string with full path and file name 
                        to a lightkurvefile object
    RETURNS:
        lcf    -(lk class object) cleaned Lightkurve LightCurveFile Object
                OR (str) 'None' if unable to locate file
    '''
    try:
        lcf = lk.open(path)
    except FileNotFoundError:
        lcf = 'None'
        pass
    return lcf

def open_lc(path):
    '''
    ~Opens cleaned (for LS) Light Curves~
    REQUIRES: lightkurve as lk
    ARGS:
        path   -(str) returned string from pathfinder function, 
                        or a string with full path and file name 
                        to a lightkurvefile object
    RETURNS:
        lc    -(lk class object) Lightkurve LightCurve Object, cleaned lc
                OR (str) 'None' if unable to locate file
    '''
    try:
        lcf = open_lcf(path)
        lc = lcf.FLUX
    except AttributeError:
        lc = 'None' 
        pass
    return lc

def header_sector(lcf):
    '''
    ~Determines sector from header of LightKurve LightCurveFile class object
    REQUIRES: lightkurve as lk
    ARgs:
        lc  -(lk class obj) lightkurve lightcurve file
    REUTRNS:
        sector of observation for given file (int)
    '''
    sector = lcf.header()['SECTOR']
    return sector
    


    
# calculates statistics

def get_rvar(flux):
    '''
    ~ Calculates Rvar from flux array~
    REQUIRES: numpy
    Args:
        flux             -(array)array of time-series flux values
    Returns:
        value of Rvar   
    '''
    R_var = np.percentile(flux, 95) - np.percentile(flux, 5)
    return R_var

def ls_measure(time,flux,flux_err=None):
    '''
    ~ Calculates Lomb-Scargle periodogram and records top 3 peaks and periods ~
    REQUIRES: starspot, numpy
    Args:
        time       -(array) light curve time array
        flux       -(array) light curve flux array
        flux_err   -(optional array or int or float) light curve flux uncertainity
    Returns:
        periods    -(list) periods corresponding to top 3 power peaks in descending order
        peaks      -(list) powers of top 3 peaks in descending order
    '''
    assert len(flux)==len(time)#find better way to do this that wont break if fails
    #creating model & gathering stats
    rotate = ss.RotationModel(time, flux, flux_err) #test if works with flux_err=None or need if/else
    #lsp
    ls_period = rotate.ls_rotation(min_period=.01 , max_period=27.) #rp1 #amend min/max periods to check for here, it will carry to freq below
    power = rotate.power #lsp y-axis array
    freq = rotate.freq #lsp x-axis array
    ps = 1./freq
    peaks = np.array([i for i in range(1, len(ps)-1) if power[i-1] < \
                power[i] and power[i+1] < power[i]]) #all lsp peaks
    peak_amps_low2high = np.sort(power[peaks])
    first_peakamp = peak_amps_low2high[-1] #amp of rp1
    second_peakamp = peak_amps_low2high[-2] #amp of rp2
    third_peakamp = peak_amps_low2high[-3] #amp of rp3
    second_rp = ps[power == second_peakamp][0] #rp2
    third_rp = ps[power == third_peakamp][0] #rp3
    periods = [ls_period,second_rp,third_rp]
    powers = [first_peakamp,second_peakamp,third_peakamp]
    return periods, powers
    




