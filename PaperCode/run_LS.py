import pandas as pd
import time as imp_t
from tqdm import tqdm
import numpy as np
import LS_modular 


###to use, change 6 total paths for opening target lists and saving df of stats

#load target lists
externalpath = '/Volumes/Seagate-stars/' #always needed
ndf = pd.read_csv('{}PAPER_FINAL_FILES/target_lists/NCVZ_targetlist.csv'.format(externalpath)) #north cvz
sdf = pd.read_csv('{}PAPER_FINAL_FILES/target_lists/SCVZ_targetlist.csv'.format(externalpath)) #south cvz
df = pd.read_csv('{}PAPER_FINAL_FILES/target_lists/secs_14_15_targetlist.csv'.format(externalpath)) #all sector 14 &15
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
# savepath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/stats/stitched_ls_statsdf.csv' #for final df with stats
##OR##
#run 2
tic_list = single_sector_targets['TIC'].to_numpy()
sectors=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
savepath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/stats/ls_statsdf.csv' #for final df with stats
#use with both runs
emergencysavepath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/stats/emergencysave_ls_statsdf.csv' #to save stats already finished if any errors occur


##for testing purposes
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
        df.to_csv('{}'.format(savepath),index=False)


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
                    sector = str(sector) #so final df has consistent data types for sec column
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
        df.to_csv('{}'.format(savepath),index=False)


    else:
        print('Sectors not understood, check formatting')
        
except: #emergency save if unexpected error
    mydata = {'TIC':tics, 'Sector':secs, 'rvar':rvar,'ls-1':ls1, 'ls-2':ls2, 'ls-3':ls3,
                 'lsamp-1':lsamp1, 'lsamp-2':lsamp2, 'lsamp-3':lsamp3}
    df = pd.DataFrame(mydata)
    df.to_csv('{}'.format(emergencysavepath),index=False)
    
    print('Unknown Error happened but saved progress before tic {} at index {}'.format(tic,count-1))
print('F I N I S H E D ')