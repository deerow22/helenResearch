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
        fullpath = glob.glob('data/LightCurves/*{}-*-s_lc.fits'.format(tic)) #to use wildcard*
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
    if len(fullpath) >0:
        for file in fullpath:
        	try:
        		lcfile = lk.open(file) #open only works one file at a time
        	except FileNotFoundError:
        		pass
        	mystring = str(type(lcfile))
        	if mystring[34:-2] == 'TessLightCurveFile': #guards against'TessTargetPixelFile'& more
        		#sectors
        		hdr = lcfile.header();
        		sector = hdr['SECTOR'];
        		sectors.append(sector)
        		#cleaning
        		lc = lcfile.PDCSAP_FLUX.normalize()
        		clean_lc = lc.remove_outliers(sigma=5)
        		#ls with ss
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
        		periods1.append(first_rp)
        		periods2.append(second_rp)
        		periods3.append(third_rp)
        		powers1.append(first_peakamp)
        		powers2.append(second_peakamp)
        		powers3.append(third_peakamp)
        		#rvar
        		flux = clean_lc.flux
        		R_var = np.percentile(flux, 95) - np.percentile(flux, 5)
        		rvars.append(R_var)
        		#saving clean lc
        		savepath = 'data/SECONDRUN/cleaned_LightCurves/{}/sector{}_lc.fits'.format(tic,sector);
        		os.makedirs(os.path.dirname(savepath), exist_ok=True); #verify/make folder with tic_id as the name
        		clean_lc.to_fits(path=savepath,overwrite=True, flux_column_name='FLUX');
        	else:
        		pass
    else:
    	sectors = []
        
    if len(sectors) == 0:
        stmt = 'found nothing for TIC: {}'.format(tic)
    else:
        stats = {'sector':sectors,'rvar_persec':rvars,'ls_rp1':periods1,'ls_rp2':periods2, \
                 'ls_rp3':periods3,'amplitude1':powers1,'amplitude2':powers2,'amplitude3':powers3}
        
        df = pd.DataFrame(data=stats)
        df.to_csv('data/SECONDRUN/cleaned_LightCurves/{}/stats_ADDITIONS.csv'.format(tic),index=False)
        stmt = 'Cleaned {} sectors for TIC: {}'.format(len(sectors),tic)

    return stmt, tic, len(sectors)
    
secs = []
ids=[]
#tics = np.load('data/cool_cvz_tics.npy')
tics = np.load('data/additional_cvz_tics.npy')

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
datas = {'IDS':ids,'NUMsectors':secs}
ddf = pd.DataFrame(data=datas)
ddf.to_csv('data/THIRDRUN/numofsectors_ADDITIONS.csv',index=False)