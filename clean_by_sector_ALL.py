import lightkurve as lk
import glob
import numpy as np
import pandas as pd
import os
import starspot as ss

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
        fullpath = glob.glob('/Volumes/Seagate-stars/LightCurves/*{}-*-s_lc.fits'.format(tic))
    else: #to change path for other computer
        pathstart = path #user defined path to datafile on their computer  
        pathstart = str(pathstart) #make a string in case user forgets to but think that gives an err anyway
        pathend = pathstart +'*{}-*-s_lc.fits'.format(tic) #stitches path & filename
        fullpath= glob.glob(pathend) #to use wildcard* 
#collects files into a lightkurve class object
    sectors =[] 
    periods1 = []
    periods2=[]
    periods3=[]
    powers1=[]
    powers2=[]
    powers3=[]
    rvars =[]
    acf_firstrp = []
    acf2 = []
    acf3 = []
    acfamp1 = []
    acfamp2 = []
    acfamp3 = []
    if len(fullpath) >0:
        for file in fullpath:
        	try:
        		lcfile = lk.open(file) #open only works one file at a time
        	except FileNotFoundError:
        		pass
        	else:
        		mystring = str(type(lcfile))
        		if mystring[34:-2] == 'TessLightCurveFile': #guards against'TessTargetPixelFile'& more
        			#sectors
        			hdr = lcfile.header();
        			sector = hdr['SECTOR'];
        			sectors.append(sector)
        			#cleaning
        			lc = lcfile.PDCSAP_FLUX.normalize()
        			clean_lc = lc.remove_outliers(sigma=5)
        			#acf
        			interval = 0.00138889
        			lags1, acf1 = ss.simple_acf(clean_lc.time,clean_lc.flux,interval,smooth=9,window_length=99,polyorder=3)
        			m = lags1 > 0
        			x2 = lags1[m]
        			y2 = acf1[m]
        			peaks = np.array([i for i in range(1, len(y2)-1) if y2[i-1] < y2[i] and \
        			y2[i+1] < y2[i]])
        			print('mytest1: what is lenght of peaks?', len(peaks))
        			print('mytest2: now whats lenght?',len(peaks))
        			if len(peaks) >= 1:
        				x_peaks = x2[peaks]
        				y_peaks = y2[peaks]
        				inds = np.argsort(y_peaks)
        				xpeaks, ypeaks = x_peaks[inds][::-1], y_peaks[inds][::-1]
        				acf_rp = xpeaks[0]
        				acf_amp1 = ypeaks[0]
        				lags = lags1[m]
        				acf = acf1[m]
        				if len(peaks) >= 2:
        					acfrp2 = xpeaks[1]
        					acf_amp2 = ypeaks[1]
        				else:
        					acfrp2 = 0
        					acf_amp2 = -9999.
        				if len(xpeaks) >= 3:
        					acfrp3 = xpeaks[2]
        					acf_amp3 = ypeaks[2]
        				else:
        					acfrp3 = 0
        					acf_amp3 = -9999.
        			#ls with ss
        			time = clean_lc.time
        			time = clean_lc.time
        			flux = clean_lc.flux
        			flux_err = clean_lc.flux_err
        			rotate = ss.RotationModel(time, flux, flux_err)
        			first_rp = rotate.ls_rotation() #rp1
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
        			#rvar
        			#flux = clean_lc.flux
        			R_var = np.percentile(flux, 95) - np.percentile(flux, 5)
        			rvars.append(R_var)
        			
        		else:
        			sector = -9999.
        			rvar = -99999.
        			sectors.append(sector)
        			rvars.append(rvar)
        			acf_rp = 0
        			acfrp2 = 0
        			acfrp3 = 0
        			acf_amp1 = -9999.
        			acf_amp2 = -9999.
        			acf_amp3 = -9999.
        			first_rp = 0
        			second_rp=0
        			third_rp =0
        			first_peakamp=-9999.
        			second_peakamp=-9999.
        			third_peakamp=-9999.
        	print('ACF Measurement is {} day rotation'.format(acf_rp))
        	acf_firstrp.append(acf_rp)
        	acfamp1.append(acf_amp1)
        	acf2.append(acfrp2)
        	acfamp2.append(acf_amp2)
        	acf3.append(acfrp3)
        	acfamp3.append(acf_amp3)
        	periods1.append(first_rp)
        	periods2.append(second_rp)
        	periods3.append(third_rp)
        	powers1.append(first_peakamp)
        	powers2.append(second_peakamp)
        	powers3.append(third_peakamp)
    else:
    	sectors = []
    if len(sectors) == 0:
    	stmt = 'found nothing for TIC: {}'.format(tic)
    else:
    	tracktics = np.repeat(tic,len(sectors))
    	stats = {'sector':sectors,'rvar_persec':rvars,'ls_rp1':periods1,'ls_rp2':periods2, \
    			'ls_rp3':periods3,'amplitude1':powers1,'amplitude2':powers2,'amplitude3':powers3, \
    			'acf_rp1':acf_firstrp, 'acf_rp2':acf2,'acf_rp3':acf3,'acf_amp1':acfamp1,'acf_amp2':acfamp2, \
    			'acf_amp3':acfamp3}
    	print('Heres my stats dictionary:', stats)
    	df = pd.DataFrame(data=stats)
    	df.to_csv('data/SECONDRUN/cleaned_LightCurves/{}/stats.csv'.format(tic),index=False)
    	stmt = 'Cleaned {} sectors for TIC: {}'.format(len(sectors),tic)
    return stmt, tic, len(sectors)







secs = []
ids=[]
#tics = np.load('data/cool_cvz_tics.npy')
#tics = np.load('data/bad_ticlistsAND_mast/additional_cvz_tics.npy')
tics = np.load('data/all_dled_tics.npy')

for count, tic in enumerate(tics):
	print('Starting Tic: {}'.format(tic))
	print('                             Working on count:', count)
	print('                             Working on count:', count)
	print('                             Working on count:', count)
	print('                             Working on count:', count)
	statement, tid, totsecs = prep_lcfs(tic,'mypath')
	ids.append(tid)
	secs.append(totsecs)
	print(statement)
datas = {'IDS':ids ,'NUMsectors':secs}
ddf = pd.DataFrame(data=datas)
ddf.to_csv('data/THIRDRUN/numofsectors_ALL.csv',index=False)