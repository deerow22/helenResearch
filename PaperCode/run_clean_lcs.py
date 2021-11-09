#working on cleaning script to have all in one place & ensure normalized

import time as imp_t
from tqdm import tqdm
import pandas as pd
import numpy as np
import clean_lcs

'''
PURPOSE: cleans raw lcs & saves them DIFFERENT from run_clean_lcs_forBLS.py b/c code below
		does have lower outlier removal for lomb-scargle analysis but this is bad for BLS 
		anaylsis b/c removes transit bottoms
INSTRUCTIONS: replace all paths throughout code to direct code to raw data and specify
		where to save cleaned data and file name desired. Current raw data file names
		based on MAST bulk download naming conventions (see clean_lcs.py locate_files() 
		function for details).

current limitations: 
- can only stitch north or south (but no stars should be in both)
- MUST do sectors 14 and 15 individually unless target is in the CVZ (b/c of my paths)
current possibilities:
- stitch all N/S sectors for one star
- clean single sectors as long as star is observed in N/S CVZ or Sectors 14 or 15

'''
####################################################################################
#########################make changes here##########################################

###always needed
externalpath = '/Volumes/Seagate-stars/' #when generalized this is sec_path to raw datafiles (w/o fn)



# #first run
# df = pd.read_csv('{}Final_Run/SCVZ_targetlist.csv'.format(externalpath)) #south
# sectors_list = [1,2,3,4,5,6,7,8,9,10,11,12,13] #south
# save_path = externalpath + 'PAPER_FINAL_FILES/CleanLCs/CLEAN_CVZs/'
# stitch = True
# 
# #second run
# df = pd.read_csv('{}Final_Run/NCVZ_targetlist.csv'.format(externalpath)) #north
# sectors_list = [14,15,16,17,18,19,20,21,22,23,24,25, 26] #north cvz
# save_path = externalpath + 'PAPER_FINAL_FILES/CleanLCs/CLEAN_CVZs/'
# stitch = True
# 
# #third run
# df = pd.read_csv('{}Final_Run/secs_14_15_targetlist.csv'.format(externalpath)) #14 & 15
# sectors_list = [14]
# save_path = externalpath + 'PAPER_FINAL_FILES/CleanLCs/CLEAN_1415/'
# stitch = False
# 
# #fourth run
# df = pd.read_csv('{}Final_Run/secs_14_15_targetlist.csv'.format(externalpath)) #14 & 15
# sectors_list = [15]
# save_path = externalpath + 'PAPER_FINAL_FILES/CleanLCs/CLEAN_1415/'
# stitch = False
# 
#fifth run
df = pd.read_csv('{}Final_Run/SCVZ_targetlist.csv'.format(externalpath)) #south
sectors_list = [1,2,3,4,5,6,7,8,9,10,11,12,13] #south
save_path = externalpath + 'PAPER_FINAL_FILES/CleanLCs/CLEAN_CVZs/'
stitch = False

# #sixth run
# df = pd.read_csv('{}Final_Run/NCVZ_targetlist.csv'.format(externalpath)) #north
# sectors_list = [16,17,18,19,20,21,22,23,24,25, 26] #north w/o 14/15 so will only do cvz/s since 14/15 done before
# save_path = externalpath + 'PAPER_FINAL_FILES/CleanLCs/CLEAN_CVZs/'
# stitch = False


###always needed
tic_list = df['TIC'].to_numpy()
tic_list = tic_list #[12598::] #b/c ncvz-stitched stopped early due to bug

####################################################################################
####################################################################################



##get paths to raw data, ideally this is generalized to above to one sec_path to all raw files
south = np.arange(1,14,1)
north = np.arange(14,27,1)
if len(sectors_list) ==1:
    #if single sector
    sector = sectors_list[0] #this line MANDATES sector 15 be done individually unless target is in cvz
    if (sector==14) |(sector==15):
        sec_path = externalpath + 'SECTORS/Sector_{}/Sec{}_rawLCs/'.format(sector,sector)
    elif np.isin(sector,north)==True:
        sec_path = externalpath + 'North_CVZ/raw_lcs/'
    elif np.isin(sector,south)==True:
        sec_path = externalpath + 'SCVZ_raw_LightCurves/'
    else:
        err_stmt ='Sector of observation outside range of integers 1-26; or no files found'
        print(err_stmt)
else:
    #if multiple sectors
    sector_arr = np.array(sectors_list)
    if np.isin(sector_arr,north).all()==True: #note must do sec 14 & sec 15 individually unless target in cvz
        sec_path = externalpath + 'North_CVZ/raw_lcs/'
    elif np.isin(sector_arr,south).all()==True:
        sec_path = externalpath + 'SCVZ_raw_LightCurves/'
    else:
        err_stmt ='WARNING: May not have cleaned all desired sectors! \nCannot mix sectors from Northern & Southern Hemispheres. Please use single sectors or integers 1-13 for south; 14-16 for north'
        print(err_stmt)
        #note this doens't stop it from running, it just uses last target's sec_path so if some secs are there will clean them
        
    

    
##cleans file(s) per tic & saves
for tic in tqdm(tic_list):
    imp_t.sleep(1)
    ##use sec_path to find raw datafiles  
    allfilepaths = clean_lcs.locate_files(tic,path=sec_path) #search for all paths
    ##open raw lcs for available sectors 
    lcfiles, secs_here = clean_lcs.sector_ordered_files(allfilepaths,sectors_list,tic) #only open files within sector_list
    ##clean & save lcs
    if stitch == True:
        if lcfiles != -99:
            savepath = save_path + '{}/stitched_lc.fits'.format(tic)
            cleanlc = clean_lcs.clean_files(lcfiles)
            clean_lcs.save_lc(cleanlc,savepath=savepath)
            print('save successful for stitched TIC {}'.format(tic))
        else:
            pass
    elif lcfiles == -99 : #catches if nofilefound before len(of an integer) throws an error
        pass
    elif len(lcfiles) > 1: #multiple sectors
        for file in lcfiles: 
            if file != -99:
                realsec = file.header()['SECTOR'] #changes b/c of sorting so this gets accurate one
                savepath = save_path + '{}/sec{}_lc.fits'.format(tic,realsec)
                cleanlc = clean_lcs.clean_files(file)
                clean_lcs.save_lc(cleanlc,savepath=savepath)
                print('save successful for TIC {} sector {}'.format(tic,realsec))
            else:
                pass
    else: #single sector
        if lcfiles != -99:
                realsec = lcfiles[0].header()['SECTOR'] #prob not needed b/c single sector but provides extra security
                savepath = save_path + '{}/sec{}_lc.fits'.format(tic,realsec)
                cleanlc = clean_lcs.clean_files(lcfiles[0])
                clean_lcs.save_lc(cleanlc,savepath=savepath)
                print('save successful for TIC {} sector {}'.format(tic,realsec))
        else:
            pass
        
print(' F I N I S H E D ')   