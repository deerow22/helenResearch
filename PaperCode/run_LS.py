import pandas as pd
import time as imp_t
from tqdm import tqdm
import numpy as np
import LS_modular 

'''
INSTRUCTIONS:
this script will run a lomb-scargle periodogram (LSP) analysis on all desired targets.
currently set to receive targets from 1) NCVZ 2)SCVZ 3)Sectors 14 & 15. Then combines
all targets into one list and saves a dataframe of results from the LSP. Results include 
highest 3 peaks' rotation periods and power amplitudes.
Has ability to analyze single sectors of data or stitched sectors of data. For stitched 
light curves set 'sectors' variable to 'stitched' or for single sectors set variable equal
to a list of desired sector(s) integer(s) value(s). You will need to change 'externalpath' 
variable to match computer system - this is location of data files for target lists, target cleaned light
curves, and where statistics dataframes will save to. 
'''

#load target lists
externalpath = '/Volumes/Seagate-stars/Final_Run' #change as needed #also used for saving stats dfs
ndf = pd.read_csv('{}/NCVZ_targetlist.csv'.format(externalpath)) #north cvz
sdf = pd.read_csv('{}/SCVZ_targetlist.csv'.format(externalpath)) #south cvz
df = pd.read_csv('{}/secs_14_15_targetlist.csv'.format(externalpath)) #all sector 14 &15
north = ndf['TIC'].to_numpy() #needed for sector naming
south = sdf['TIC'].to_numpy() #needed to check tic hemisphere
#merge target lists for efficiency
stitched_targets = pd.concat([ndf,sdf])
single_sector_targets = pd.concat([ndf,sdf,df])

############################################################################
########################make changes here ##################################

##choose targets here

#run 1
# tic_list = stitched_targets['TIC'].to_numpy()
# sectors = 'stitched'
##OR##
#run 2
tic_list = single_sector_targets['TIC'].to_numpy()
sectors=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]


##for testing
# tic_list = [199682037,  33733169, 123,4]
# sectors = 'stitched'
# sectors=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
############################################################################
############################################################################
#empty lists to save data to
tics=[]; rvar=[]; secs=[]
lsamp1=[]; lsamp2=[]; lsamp3=[]
ls1=[]; ls2=[]; ls3=[]
count=0
try:
    if sectors == 'stitched': #NCVZ & SCVZ
        #get stats for targets
        for tic in tqdm(tic_list):
            count+=1
            imp_t.sleep(1)
            #find data
            mypath = LS_modular.pathfinder(tic,sectors)
            lc = LS_modular.open_lc(mypath)
            if lc != 'None':
                #get rvar & sector
                Rvar = LS_modular.get_rvar(lc.flux)
                if np.isin(tic,north)==True:
                    sector = 'NCVZ-stitched'
                elif np.isin(tic,south)==True:
                    sector = "SCVZ-stitched"
                else:
                    sector = 'None' #safety check
                #get ls
                rps, amps = LS_modular.ls_measure(lc.time,lc.flux,lc.flux_err)
                #append results
                rvar.append(Rvar);ls1.append(rps[0]); ls2.append(rps[1]); ls3.append(rps[2])
                lsamp1.append(amps[0]); lsamp2.append(amps[1]); lsamp3.append(amps[2])
                tics.append(tic);secs.append(sector)
            else:
                pass
        #save data to a df
        mydata = {'TIC':tics, 'Sector':secs, 'rvar':rvar,'ls-1':ls1, 'ls-2':ls2, 'ls-3':ls3,
                 'lsamp-1':lsamp1, 'lsamp-2':lsamp2, 'lsamp-3':lsamp3}
        df = pd.DataFrame(mydata)
        df.to_csv('{}/stitched_ls_statsdf.csv'.format(externalpath),index=False)


    elif type(sectors) == type(['list']): #single sectors
        #do ls for targets
        for tic in tqdm(tic_list):
            count+=1
            imp_t.sleep(1)
            paths = []
            for sec in sectors:
                mypath = LS_modular.pathfinder(tic,sec)
                paths.append(mypath)
            for p in paths:
                lcf = LS_modular.open_lcf(p)
                if lcf != 'None':
                    #open data
                    lc = lcf.FLUX
                    #get rvar & sector
                    Rvar = LS_modular.get_rvar(lc.flux)
                    sector = lcf.header()['SECTOR']
                    #get lsp
                    rps, amps = LS_modular.ls_measure(lc.time,lc.flux,lc.flux_err)
                    #append results
                    rvar.append(Rvar);ls1.append(rps[0]); ls2.append(rps[1]); ls3.append(rps[2])
                    lsamp1.append(amps[0]); lsamp2.append(amps[1]); lsamp3.append(amps[2])
                    tics.append(tic);secs.append(sector)
                else:
                    pass
        #save data to a df
        mydata = {'TIC':tics, 'Sector':secs, 'rvar':rvar,'ls-1':ls1, 'ls-2':ls2, 'ls-3':ls3,
                 'lsamp-1':lsamp1, 'lsamp-2':lsamp2, 'lsamp-3':lsamp3}
        df = pd.DataFrame(mydata)
        df.to_csv('{}/ls_statsdf.csv'.format(externalpath),index=False)


    else:
        print('Sectors not understood, check formatting')
        
except: #emergency save if unexpected error
    mydata = {'TIC':tics, 'Sector':secs, 'rvar':rvar,'ls-1':ls1, 'ls-2':ls2, 'ls-3':ls3,
                 'lsamp-1':lsamp1, 'lsamp-2':lsamp2, 'lsamp-3':lsamp3}
    df = pd.DataFrame(mydata)
    df.to_csv('{}/emergencysave_ls_statsdf.csv'.format(externalpath),index=False)
    
    print('Unknown Error happened but saved progress before tic {} at index {}'.format(tic,count-1))
print('F I N I S H E D ')