import lightkurve as lk
import glob
import os
import numpy as np
import starspot as ss
from astropy.table import Table
import pandas as pd
from starspot import sigma_clipping


###defing fcn to clean and stitch full lightcurve #############################
def prep_lcfs(tic,path=None):
    """ 
    Locates all TESS lightcurve files with filenames formatted from a mast bulk dl.
    Does minor cleaning & stitching, then plots final lightcurve.
    REQUIRES: lightkurve, glob, os


    
    Args:
        tic: TIC identification number (integer or string).
        path: (optional) computer path to file location (string) (excluding filename).
        
    Returns:
        cleaned & stitched lightkurve object.
    """
#locates files
    if path == None: #if only need filename
        fullpath = glob.glob('*{}-*-s_lc.fits'.format(tic)) #to use wildcard*
    elif path == 'mypath': #path on my computer
        #fullpath = glob.glob('data/LightCurves/*{}-*-s_lc.fits'.format(tic)) #to use wildcard*
        fullpath = glob.glob('data/Added_LightCurves/*{}-*-s_lc.fits'.format(tic)) #to catch those new added tics
    else: #to change path for other computer
        pathstart = path #user defined path to datafile on their computer  
        pathstart = str(pathstart) #make a string in case user forgets to but think that gives an err anyway
        pathend = pathstart +'*{}-*-s_lc.fits'.format(tic) #stitches path & filename
        fullpath= glob.glob(pathend) #to use wildcard* 
#collects files into a lightkurve class object
    lcfs =[] 
    sectorspresent = []
    probs=[]
    mylist = [1,2,3,4,5,6,7,8,9,10,11,12,13]
    for i in mylist:
        for file in fullpath:
            try:
            	lcfile = lk.open(file) #open only works one file at a time
            	mystring = str(type(lcfile))
            	if mystring[34:-2] == 'TessLightCurveFile': #guards against'TessTargetPixelFile'& more
                	hdr = lcfile.header()
                	sector = hdr['SECTOR']
                	if i == sector:
                		lcfs.append(lcfile) #collect
                		sectorspresent.append(sector)
                	else:
                		pass
            except ValueError:
            	pass
    if len(lcfs)==0: 
        print('Unable to locate any files for TIC {}; verify path/files exist'.format(tic))
        cleaned = -99
        number_ofsectors=0 #use to match with tic list and remove tics with no info
        usedsectors = 0
    else:
        number_ofsectors = int(len(lcfs))
        usedsectors = sectorspresent #0
#cleans & stitches #& plots
        lcfiles = lk.collections.LightCurveFileCollection(lcfs) #making list into class collection
        stitched = lcfiles.PDCSAP_FLUX.stitch() #this detrends/normalizes each sector before stitching together all lightcurves
        #could mask bad bits here then stitch, instead of stitch above
        nonans = stitched.remove_nans() #preps for periodogram etc
        cleaned = nonans.remove_outliers(sigma=5) #removes cosmic rays & flares mostly
##############################################################################################
#need to generalize for saving
        filename = 'data/SECONDRUN/cleaned_LightCurves/{}/lc.fits'.format(tic)
        os.makedirs(os.path.dirname(filename), exist_ok=True) #verify/make folder with tic_id as the name
        cleaned.to_fits(filename,flux_column_name = 'flux',overwrite=True); #save cleaned to file
#############################################################################################
    if number_ofsectors == 0:
        print('nothing found for: TIC',tic)
        return cleaned, number_ofsectors, usedsectors
    else:
        return cleaned, number_ofsectors, usedsectors #target final data table containing flux, time, flux_err 



#open target list
#### requires numpy 
#cvz_tics = np.load('data/unique_cvz_tics.npy')
#cvz_tics = np.load('data/cool_cvz_tics.npy')
cvz_tics = np.load('data/additional_cvz_tics.npy')

#measure rotation periods  
### requires starrotate, Table     
tics=[]
rvar_smoothed=[]
rvar_orig=[]
ls=[]
firstamps=[]
secondamps=[]
thirdamps=[]
ls2=[]
ls3=[]
acf=[]
acf2=[]
acf3=[]
acfamp1=[]
acfamp2=[]
acfamp3=[]
pdm=[]
stitched_sectors=[]

for count, tic in enumerate(cvz_tics):
	print('currently cleaning #', count)
	print('currently cleaning #', count)
	print('currently cleaning #', count)
	print('currently cleaning #', count)
	print('                 aka this star has TIC:',tic, '  and count:  ',count)
	print('                 aka this star has TIC:',tic, '  and count:  ',count)
	#test
	lc, totalsectors, sectorslist = prep_lcfs(tic,'mypath')
	if str(type(lc)) == "<class 'lightkurve.lightcurve.TessLightCurve'>":
#tic order
		tics.append(tic)
		np.save('data/SECONDRUN/tic_order_ADDITIONS',tics)
#
#rvar-original
		R_var = np.percentile(lc.flux, 95) - np.percentile(lc.flux, 5)
		rvar_orig.append(R_var)
		np.save('data/SECONDRUN/rvar_original_ADDITIONS',rvar_orig)
    	
# start Ruth tutorial code starspot
# 		Calculate the median so that we can median-normalize.
# 		med = np.median(lc.flux)
# 
# 		Do an initial sigma clip to remove big outliers.
# 		m = sigma_clipping.sigma_clip(lc.flux/med - 1, nsigma=6)
# 		x, y, yerr = lc.time[m], lc.flux[m]/med - 1, lc.flux_err[m]/med
# 
# 		Then a sigma clip using a Sav-Gol filter for smoothing
# 		smooth, mask = sigma_clipping.filter_sigma_clip(x, y, window_length=199)
# 		time, flux, flux_err = x[mask], y[mask], yerr[mask]
# end Ruth tutorial code starspot
# 
# rvar-smoothed
# 		R_var2 = np.percentile(flux, 95) - np.percentile(flux, 5)
# 		rvar_smoothed.append(R_var2)
# 		np.save('data/SECONDRUN/rvar_smoothed_ADDITIONS',rvar_smoothed)
# 
# 		rotate = ss.RotationModel(time, flux, flux_err) 
# ls rp
# 		ls_period = rotate.ls_rotation()
# 		print('Lomb-Scargle Measurement is {} day rotation'.format(ls_period))
# 		ls.append(ls_period)
# 		np.save('data/SECONDRUN/ls_rps_ADDITIONS', ls)
# 		#getting peak-heights-ls
# 		power = rotate.power
# 		freq = rotate.freq
# 		ps = 1./freq
# 		peaks = np.array([i for i in range(1, len(ps)-1) if power[i-1] < \
# 						power[i] and power[i+1] < power[i]])
# 		peak_amps_low2high = np.sort(power[peaks])
# #amplitudes of ls top 3 peaks/rps
# 		first_peakamp = peak_amps_low2high[-1]
# 		firstamps.append(first_peakamp)
# 		np.save('data/SECONDRUN/ls_amps_ADDITIONS',firstamps)
# 		second_peakamp = peak_amps_low2high[-2]
# 		secondamps.append(second_peakamp)
# 		np.save('data/SECONDRUN/second_ls_amp_ADDITIONS',secondamps)
# 		third_peakamp = peak_amps_low2high[-3]
# 		thirdamps.append(third_peakamp)
# 		np.save('data/SECONDRUN/third_ls_amp_ADDITIONS',thirdamps)
# #second/third highest ls rp
# 		second_rp = ps[power == second_peakamp][0]
# 		ls2.append(second_rp)
# 		np.save('data/SECONDRUN/second_ls_rp_ADDITIONS',ls2) #ls-rp2
# 		third_rp = ps[power == third_peakamp][0]
# 		ls3.append(third_rp)
# 		np.save('data/SECONDRUN/third_ls_rp_ADDITIONS',ls3) #ls-rp3
# 
		
#total num sectors
		stitched_sectors.append(totalsectors)
		np.save('data/SECONDRUN/stitched_sectors_ADDITIONS',stitched_sectors)
		print('sectors:',totalsectors,'which sectors:',sectorslist)
	else:
		print('Not correct type OR nothing found for TIC:', tic, 'at count:',count,'type:', type(lc))
		continue
# 		
# 		
# 		
# ##########################2nd half#######################################################		
# 		
# 		
# 	if str(type(lc)) == "<class 'lightkurve.lightcurve.TessLightCurve'>":
# 		print('TIC-acf:', tic)
# 		##rotate = ss.RotationModel(time,flux,flux_err)
# 		##tess_cadence = 1./24./30.
# #start Ruth tutorial code starspot
# 		# Calculate the median so that we can median-normalize.
# 		med = np.median(lc.flux)
# 
# 		# Do an initial sigma clip to remove big outliers.
# 		m = sigma_clipping.sigma_clip(lc.flux/med - 1, nsigma=6)
# 		x, y, yerr = lc.time[m], lc.flux[m]/med - 1, lc.flux_err[m]/med
# 
# 		# Then a sigma clip using a Sav-Gol filter for smoothing
# 		smooth, mask = sigma_clipping.filter_sigma_clip(x, y, window_length=199)
# 		time, flux, flux_err = x[mask], y[mask], yerr[mask]
# #end Ruth tutorial code starspot
# 
# #acf rps & amps
# 		##acf_period = rotate.acf_rotation(tess_cadence)
# 		interval = 0.00138889
# 		lags1, acf1 = simple_acf(time,flux,interval,smooth=9,window_length=99,polyorder=3)
# 		m = lags1 > 0
# 		x2 = lags1[m]
# 		y2 = acf1[m]
# 		
# 		peaks = np.array([i for i in range(1, len(y2)-1) if y2[i-1] < y2[i] and \
# 						y2[i+1] < y2[i]])
# 		if len(peaks) >= 1:
# 			x_peaks = x2[peaks]
# 			y_peaks = y2[peaks]
# 			inds = np.argsort(y_peaks)
# 			xpeaks, ypeaks = x_peaks[inds][::-1], y_peaks[inds][::-1]
# 			acf_rp = xpeaks[0]
# 			acf_amp1 = ypeaks[0]
# 			lags = lags1[m]
# 			acf = acf1[m]
# 			if len(peaks) >= 2:
# 				acfrp2 = xpeaks[1]
# 				acf_amp2 = ypeaks[1]
# 			else:
# 				acfrp2 = 0
# 				acf_amp2 = -9999.
# 			if len(xpeaks) >= 3:
# 				acfrp3 = xpeaks[2]
# 				acf_amp3 = ypeaks[2]
# 			else:
# 				acfrp3 = 0
# 				acf_amp3 = -9999.
# 		else:
# 			acf_rp = 0
# 			acfrp2 = 0
# 			acfrp3 = 0
# 			acf_amp1 = -9999.
# 			acf_amp2 = -9999.
# 			acf_amp3 = -9999.
# 		print('ACF Measurement is {} day rotation'.format(acf_rp))
# 		acf.append(acf_rp)
# 		np.save('data/SECONDRUN/acf_rps_ADDITIONS',acf)
# 		acfamp1.append(acf_amp1)
# 		np.save('data/SECONDRUN/acf_amps_ADDITIONS',acfamp1)
# 		acf2.append(acfrp2)
# 		np.save('data/SECONDRUN/second_acf_rps_ADDITIONS',acf2)
# 		acfamp2.append(acf_amp2)
# 		np.save('data/SECONDRUN/second_acf_amp_ADDITIONS',acfamp2)
# 		acf3.append(acfrp3)
# 		np.save('data/SECONDRUN/third_acf_rps_ADDITIONS',acf3)
# 		acfamp3.append(acf_amp3)
# 		np.save('data/SECONDRUN/third_acf_amp_ADDITIONS',acfamp3)
# 	else:
# 		continue
# 		
# 
#     
#     
# stats = Table([tics,stitched_sectors,rvar_orig,rvar_smoothed,ls,ls2,ls3,firstamps,secondamps, \
# 			thirdamps,acf,acf2,acf3, acfamp1,acfamp2,acfamp3],names=('ID','TOTALnum_sectors', \
# 			'Rvar','Rvar-smoothedlc','ls_rp1','ls_rp2','ls_rp3','ls_amp1','ls_amp2','ls_amp3', \
# 			'acf_rp1','acf_rp2','acf_rp3','acf_amp1','acf_amp2','acf_amp3'))
# stats.write('data/SECONDRUN/stats_table_ADDITIONS.fits', format='fits',overwrite=True)