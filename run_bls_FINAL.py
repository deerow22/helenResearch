from astropy.table import Table
from astropy.timeseries import BoxLeastSquares
import lightkurve as lk
from astropy.io import ascii
from astropy.io import fits
import numpy as np

########functions##########

def open_lc(tic,path):
    '''
    Opens lightkurvefile object as a lightkurve file...
    ...a work around b/c cleaned_lcs saved as lkfiles lose the flux attribute
    
    Requirements: astropy.io.fits; lightkurve;
    inputs:   tic (int or str)
              path (str); '"mypath" is for me only, all others must specify filepath'
    outputs:  lightkurve lightcurve object
              additional output - writes text file of tics unable to open & exceptions (this path is hardcoded so only works for me)
    '''
    try:
        if path == 'mypath':
            filepath = 'data/SECONDRUN/cleaned_LightCurves/{}/lc.fits'.format(tic)
        else:
            filepath = path
        openfile = fits.open(filepath)
        filedata = openfile[1].data
        lc = lk.lightcurve.TessLightCurve(flux=filedata['FLUX'],time=filedata['TIME'],flux_err=filedata['FLUX_ERR'])
    except Exception as e:
        lc = 'None for TIC: {}'.format(tic)
        print('TIC: {} failed to open as lk lc'.format(tic))
        stmt = 'TIC: {} had error {}\n'.format(tic,e)
        with open("data/FOURTHRUN/BLS_couldnt_openfile.txt","a") as text_file:
            text_file.write("\n")
            text_file.write(stmt)
    return lc


def bls_stats(tic,lc):
    '''
    Runs BoxLeastSquares(BLS) on a light curve and returns BLS statistics.
    
    REQUIRES: numpy, astropy.timeseries.BoxLeastSquares,
    Args:   tic id      -(int or str)
            lc          -(lightkurve light curve object)
            plot        -(boolean) optional, if True plots 
                          periodogram and phase folded light curve 
                          and binned phase folded light curve
                    
    outputs: bls_stats  -(dictionary) Box Least Squares statistics 
                          of 'power', 'period', 'depth', 'duration', 'transit_time' 
    '''
    assert(str(type(lc)) == "<class 'lightkurve.lightcurve.TessLightCurve'>");
    ###gets cleaned data
    clean = lc.flatten()
    time = clean.time - clean.time[0]
    flux = clean.flux
    fluxerr = clean.flux_err
    ###gets periodogram based on bls
    from astropy.timeseries import BoxLeastSquares
    model = BoxLeastSquares(time, flux, dy=fluxerr)
    period_grid = np.linspace(.5,20,5000)
    duration_grid = np.linspace(.01,.49,3)
    periodogram = model.power(period_grid, duration_grid) 
    ppower = np.argmax(periodogram.power)
    pperiod = periodogram.period[np.argmax(periodogram.power)]
    depth = periodogram.depth[np.argmax(periodogram.power)]
    duration = periodogram.duration[np.argmax(periodogram.power)]
    tt = periodogram.transit_time[np.argmax(periodogram.power)]
    ###writes data to files
    bls_y = periodogram.power
    bls_x = periodogram.period
    try:
        np.save('data/FOURTHRUN/data_arrs/{}/blsx_periods'.format(tic),bls_x)
        np.save('data/FOURTHRUN/data_arrs/{}/blsy_powers'.format(tic),bls_y)
    except Exception as err:
        print('failure saving arrays')
        stmt2 = 'TIC: {} had error {}\n'.format(tic,err)
        with open("data/FOURTHRUN/BLS_couldnt_savearrays.txt","a") as text_file2:
            text_file2.write("\n")
            text_file2.write(stmt2)

    return tic, ppower, pperiod, depth, tt, duration #bls_stats #power, period, depth, transit_time, duration,periodogram




##########start on targets##########

tics = np.load('data/all_dled_tics.npy')
# tics = [149603524, 9999999, 25063396, 149603524, 9999999] #for testing



t = Table()
t['tic'] = [0]; t['power'] = [0.]; t['period']=[0.];t['depth']=[0.];t['transit_time']=[0.];t['duration']=[0.]
for count,tic in enumerate(tics):
#     count += 8821
    print('Starting on TIC: {} ....................which is {} of {}'.format(tic, count,len(tics)))
    lc = open_lc(tic,'mypath')
    if str(type(lc)) == "<class 'lightkurve.lightcurve.TessLightCurve'>":
        print('Starting BLS')
        stats = bls_stats(tic,lc)
        t.add_row((stats[0],stats[1],stats[2],stats[3],stats[4],stats[5]))
        ascii.write(t, 'data/FOURTHRUN/BLS_stats_table3.csv',format='csv', overwrite=True)
        print('Finished BLS stats')
    else:
        pass