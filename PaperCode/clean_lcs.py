import lightkurve as lk
import glob
import os


def locate_files(tic,path=None):
    '''
    ~ Locates TESS lightcurve files with filenames formatted from a mast bulk download.~
    REQUIRES: glob
    Args: 
        tic            -(int or str)TESS TIC ID
        path           -(str) path on computer to file(s) location
    Returns:
        list of path strings for all files found with specified tic
    '''
    if path == None: #if only need filename
        fullpath = glob.glob('*{}-*-s_lc.fits'.format(tic)) #to use wildcard*
    else: #user defined path to datafile on their computer
        pathstart = path   
        pathstart = str(pathstart) #make a string in case user forgets to but think that gives an err anyway
        pathend = pathstart +'*{}-*-s_lc.fits'.format(tic) #stitches path & filename
        fullpath= glob.glob(pathend) #to use wildcard* 
    return fullpath

def sector_ordered_files(filepaths,sectors_list,tic):
    '''
    ~ Opens multiple lightcurve files as Light Kurve 'TessLightCurveFile' Objects and sorts them by sector number.
                    Warning: this function only good for TESS ~
    REQUIRES: lightkurve
    Args:
        filepaths           -(str) list of paths to files' location on computers
        sectors_list        -(array of integers) list of expected sectors to check for sorting
                                Warning: file will be discarded if sector in file header is NOT present in sectors_list
    Returns:
        sorted list of lightkurve lightkurvefile objects, list of sectors present in returned lightkurve file objects
    '''
    lcfs =[] 
    sectorspresent = []
    for i in sectors_list:
        for file in filepaths:
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
    	lcfs = -99
    	sectorspresent = -99
    else:
    	pass
    
    return lcfs, sectorspresent


def clean_files(lightkurvefiles):
    '''
    ~ Stitches LightKurveFiles objects if more than one, then selects PDCSAP_FLUX, removes nans, and removes outliers with 5sigma tolerance~
    REQUIRES: lightkurve
    Args: 
        lightkurvefiles      -(lk object)TESS light curve files for ONE 2min cadence target at a time
    Returns:
        cleaned lightkurve lightcurve object
    '''
    if type(lightkurvefiles) is not list: lightkurvefiles = [ lightkurvefiles ] #enables list len check, then use [0] to access files after if/else
    lcfs = lightkurvefiles 
    
    if len(lcfs) ==1: #dont stitch or put in collection if single sector
        detrend = lcfs[0].PDCSAP_FLUX.normalize() #this detrends/normalizes each sector
        nonans = detrend.remove_nans() #preps for periodogram etc
    else: 
        #cleans & stitches 
        lcfiles = lk.collections.LightCurveFileCollection(lcfs) #making list into class collection
        stitched = lcfiles.PDCSAP_FLUX.stitch() #this detrends/normalizes each sector before stitching together all lightcurves
        #could mask bad bits here then stitch, instead of stitch above
        nonans = stitched.remove_nans() #preps for periodogram etc
        
    cleaned = nonans.remove_outliers(sigma=5) #removes cosmic rays & flares mostly
    return cleaned
    
def clean_files_forbls(lightkurvefiles):
    '''
    ~ Optimized for pre-BLS analysis. Stitches LightKurveFiles objects if more than one, then selects PDCSAP_FLUX, 
        removes nans, and removes upwards outliers with 5sigma tolerance while preserving transit depths~
    REQUIRES: lightkurve
    Args: 
        lightkurvefiles      -(lk object)TESS light curve files for ONE 2min cadence target at a time
    Returns:
        cleaned lightkurve lightcurve object
    '''
    if type(lightkurvefiles) is not list: lightkurvefiles = [ lightkurvefiles ] #enables list len check, then use [0] to access files after if/else
    lcfs = lightkurvefiles 
    
    if len(lcfs) ==1: #dont stitch or put in collection if single sector
        detrend = lcfs[0].PDCSAP_FLUX.normalize() #this detrends/normalizes each sector
        nonans = detrend.remove_nans() #preps for periodogram etc
    else: 
        #cleans & stitches 
        lcfiles = lk.collections.LightCurveFileCollection(lcfs) #making list into class collection
        stitched = lcfiles.PDCSAP_FLUX.stitch() #this detrends/normalizes each sector before stitching together all lightcurves
        #could mask bad bits here then stitch, instead of stitch above
        nonans = stitched.remove_nans() #preps for periodogram etc
        
    cleaned = nonans.remove_outliers(sigma_upper=5, sigma_lower=float('inf')) #removes cosmic rays & flares mostly but perserves transits
    
    return cleaned

def save_lc(data,savepath=None):
    '''
    ~ Writes lightkurve object data to file with specified path, creating directories if needed by path~
    REQUIRES: os, lightkurve
    Args:
        data        -(lk lightcurve object) 
        savepath        -(str) path with new file name for saving on computer
    Returns:
        nothing - just saves the file to specified location
    '''
    #need to generalize for saving
    filename = savepath
    os.makedirs(os.path.dirname(filename), exist_ok=True) #verify/make folder with tic_id as the name
    data.to_fits(filename,flux_column_name = 'flux',overwrite=True); #save cleaned to file


