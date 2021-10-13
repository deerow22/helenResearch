
import glob
from astropy.io import fits
import lightkurve as lk
import numpy as np
import starspot as ss

# open only sector1 lcs---keeping so doesnt break things but replacing with more modular
def open_mylcs(tic,sector,path='mypath'): #### path to cvzs only for now
    '''
    ~Opens previously cleaned light curve file--*might* have to be LightKurve object~
    Requires:       astropy.io.fits;
    Args:           tic         -(int) TESS TIC ID
                    sector      -(int) light curve sector number for target
    Returns:        time, flux, flux_err (numpy arrays) 
    '''

    if path == 'mycvzpath':
        filepath = '/Volumes/Seagate-stars/SECONDRUN/cleaned_LightCurves/{}/sector{}_lc.fits'.format(tic,sector) #cvz targets' sectors
    elif path =='mypath':
        filepath = '/Volumes/Seagate-stars/SECTORS/Sector_{}/Sec{}_cleaned/{}/lc.fits'.format(sector,sector,tic) #targets in sectors 14/15
    else:
        filepath = path
    try:
        lc = fits.open(filepath)
        lc_data = lc[1].data
        time=lc_data.TIME; flux=lc_data.FLUX; flux_err=lc_data.FLUX_ERR
    except Exception as e:
        print('TIC: {} Sector: {} couldnt be opened; encountered Error: {}'.format(tic,sector,e))
        time,flux,flux_err = 'None','None','None'
    
    return time,flux,flux_err

def pathfinder(tic,sector):
    '''
    ~Determines correct path to cleaned data in FINAL_RUN~
    REQUIRES: numpy as np
    ARGS:
        tic      -(int) TICID of target
        sector   -(str or int) 
                      'stitched' for stitched cvz stars (N or S)
                      or integer of single observation sector 
    RETURNS:
        path
    '''
    north = np.arange(13,27,1)
    south = np.arange(1,14,1)
    err_stmt = 'Sector type not understood. Specify observation sector(s) integer(s) in a list or use "stitched" for stitched sectors'
    if type(sector) == type('string'):
        if sector == 'stitched':
            path = '/Volumes/Seagate-stars/Final_Run/Clean_LCs_forBLS/{}/stitched_lc.fits'.format(tic)
        else:
            print(err_stmt)
            path = 'None'
    elif type(sector) ==int: 
        path = '/Volumes/Seagate-stars/Final_Run/Clean_LCs_forBLS/{}/sec{}_lc.fits'.format(tic,sector)
    else:
        print(err_stmt)
        path = 'None'
    return path

def pathfinder2(tic,stitched=False):
    '''
    ~Finds path to cleaned data in FINAL_RUN~
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
        fullpath = '/Volumes/Seagate-stars/Final_Run/Clean_LCs_forBLS/{}/stitched_lc.fits'.format(tic)
        path = glob.glob(fullpath)
    elif stitched ==False: 
        fullpath = '/Volumes/Seagate-stars/Final_Run/Clean_LCs_forBLS/{}/sec*_lc.fits'.format(tic)
        path = glob.glob(fullpath) #, recursive = True)
    else:
        print(err_stmt)
        path = 'None'
    return path






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
    ~Trys to open cleaned Light Curves~
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
    #ls
    ls_period = rotate.ls_rotation(min_period=.01 , max_period=27.) #rp1 #amend min/max periods to check for here, it will carry to freq below
    power = rotate.power
    freq = rotate.freq
    ps = 1./freq
    peaks = np.array([i for i in range(1, len(ps)-1) if power[i-1] < \
                power[i] and power[i+1] < power[i]])
    peak_amps_low2high = np.sort(power[peaks])
    first_peakamp = peak_amps_low2high[-1] #amp of rp1
    second_peakamp = peak_amps_low2high[-2] #amp of rp2
    third_peakamp = peak_amps_low2high[-3] #amp of rp3
    second_rp = ps[power == second_peakamp][0] #rp2
    third_rp = ps[power == third_peakamp][0] #rp3
    periods = [ls_period,second_rp,third_rp]
    powers = [first_peakamp,second_peakamp,third_peakamp]
    return periods, powers
    



