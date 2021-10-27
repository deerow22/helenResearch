import bls_modular
import LS_modular 
import pandas as pd
from astropy.io import ascii
from astropy.table import Table
import time as imp_t
from tqdm import tqdm
import numpy as np



##opening kepler eb files


#load target lists
externalpath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/' #always needed
ndf = pd.read_csv('{}target_lists/NCVZ_targetlist.csv'.format(externalpath)) #north cvz
sdf = pd.read_csv('{}target_lists/SCVZ_targetlist.csv'.format(externalpath)) #south cvz
df = pd.read_csv('{}target_lists/secs_14_15_targetlist.csv'.format(externalpath)) #all sector 14 &15
north = ndf['TIC'].to_numpy() #needed for sector naming
south = sdf['TIC'].to_numpy() #needed to check tic hemisphere
#merge target lists for efficiency
stitched_targets = pd.concat([ndf,sdf])
single_sector_targets = pd.concat([ndf,sdf,df])
##################################################make choices here
#run 1
tic_list = stitched_targets['TIC'].to_numpy()#[4039::]
sectors = 'stitched'
##OR##
#run 2
# tic_list = single_sector_targets['TIC'].unique()
# sectors='singles'

#for testing
# tic_list = [199682037,  33733169, 123,4] #NCVZ; SCVZ; none;none
# sectors = 'stitched'
# sectors = 'singles'
##################################################


##bls needs
pgrid =  bls_modular.period_grid(log='on') #np.concatenate((next1,abc))
dgrid = bls_modular.duration_grid(start=.01,stop=.19,N=5)
Secsecs=[]
Secids=[];Secpowers=[];Secperiods=[];Secdepths=[];Secdurs=[];Sectts=[]

if sectors == "stitched":
    for tic in tqdm(tic_list):
        imp_t.sleep(1)
        mypath = LS_modular.pathfinder2(tic,stitched=True)
        if len(mypath) >0: #avoid empty path list
            lc = LS_modular.open_lc(mypath[0])
            t = lc.time; f=lc.flux; fe=lc.flux_err
            flat_time, flat_flux, flat_fluxerr = bls_modular.flatten_lc(t,f,fe) #flatten [time,flux,fluxerr]
            periodogram, stats = bls_modular.bls(pgrid,dgrid,flat_time, flat_flux, flat_fluxerr)
            np.save('{}CleanLCs/Clean_LCs_BLS/{}/stitched_bls_pg_period'.format(externalpath,tic),periodogram.period)
            np.save('{}CleanLCs/Clean_LCs_BLS/{}/stitched_bls_pg_power'.format(externalpath,tic),periodogram.power)
            if np.isin(tic,north)==True:
                sector = 'NCVZ-stitched'
            elif np.isin(tic,south)==True:
                sector = "SCVZ-stitched"
            else:
                sector = 'None' #safety check
            Secids.append(tic);Secsecs.append(sector)
            Secpowers.append(stats[0]);Secperiods.append(stats[1]);Secdepths.append(stats[2]);
            Secdurs.append(stats[3]);Sectts.append(stats[4])
            #save arrays as they come in case of crash - but all in final table
            fn = '{}stats/blsarrs_stitched/'.format(externalpath)
            np.save('{}BLS_ticids.npy'.format(fn),Secids)
            np.save('{}BLS_powers.npy'.format(fn),Secpowers)
            np.save('{}BLS_periods.npy'.format(fn),Secperiods)
            np.save('{}BLS_depths.npy'.format(fn),Secdepths)
            np.save('{}BLS_durations.npy'.format(fn),Secdurs)
            np.save('{}BLS_transit_times.npy'.format(fn),Sectts)

        else:
            print(f'TIC {tic} had no stitched file found')
    #make a table
    print('making final table: bls_stats_{}.csv'.format(sectors))
    sec_table = Table([Secids,Secsecs,Secperiods,Secpowers,Secdepths,Secdurs,Sectts],
        names=('TIC','Sector','period_bls','power_bls','depth_bls','dur_bls','tt_bls'))
    ascii.write(sec_table, '{}stats/bls_stats_stitched.csv'.format(externalpath),format='csv', overwrite=True) 
    
    print(' F  I  N  I  S  H  E  D')
    
    
elif sectors == 'singles': #single sectors
    for tic in tqdm(tic_list):
        imp_t.sleep(1)
        mypaths = LS_modular.pathfinder2(tic,stitched=False)
        if len(mypaths) >0:
            for p in mypaths:
                lcf = LS_modular.open_lcf(p) #open file
                sector = lcf.header()['SECTOR']
                lc = lcf.FLUX
                t=lc.time; f=lc.flux; fe=lc.flux_err
                flat_time, flat_flux, flat_fluxerr = bls_modular.flatten_lc(t,f,fe) #flatten [time,flux,fluxerr]
                periodogram, stats = bls_modular.bls(pgrid,dgrid,flat_time, flat_flux, flat_fluxerr)
                np.save('{}CleanLCs/Clean_LCs_BLS/{}/sec{}_bls_pg_period'.format(externalpath,tic,sector),periodogram.period)
                np.save('{}CleanLCs/Clean_LCs_BLS/{}/sec{}_bls_pg_power'.format(externalpath,tic,sector),periodogram.power)
                Secids.append(tic);Secsecs.append(sector)
                Secpowers.append(stats[0]);Secperiods.append(stats[1]);Secdepths.append(stats[2]);
                Secdurs.append(stats[3]);Sectts.append(stats[4])
                    #save arrays as they come in case of crash - but all in final table
                fn = '{}stats/blsarrs_sectors/'.format(externalpath)
                np.save('{}BLS_ticids'.format(fn),Secids)
                np.save('{}BLS_sectors'.format(fn),Secsecs)
                np.save('{}BLS_powers'.format(fn),Secpowers)
                np.save('{}BLS_periods'.format(fn),Secperiods)
                np.save('{}BLS_depths'.format(fn),Secdepths)
                np.save('{}BLS_durations'.format(fn),Secdurs)
                np.save('{}BLS_transit_times'.format(fn),Sectts)
        else:
            print(f'No sectors found for TIC {tic}')
    sec_table = Table([Secids,Secsecs,Secperiods,Secpowers,Secdepths,Secdurs,Sectts],
        names=('TIC','Sector','period_bls','power_bls','depth_bls','dur_bls','tt_bls'))

    ascii.write(sec_table, '{}stats/bls_stats_sectors.csv'.format(externalpath,sector,sector),format='csv', overwrite=True) 
    print(' F  I  N  I  S  H  E  D')
    
    
else:
    print('stitch setting error, use sectors="singles" or sectors="stitched" ')
    