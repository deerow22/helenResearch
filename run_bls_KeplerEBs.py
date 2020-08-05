import bls_modular 
import pandas as pd
from astropy.io import ascii
from astropy.table import Table
import time as imp_t
from tqdm import tqdm
import numpy as np



##opening kepler eb files
kep_table = pd.read_csv('data/KeplerEBs/KeplerEB_bls_stats.csv') #a subselection with <15 day periods
kics = kep_table['KIC'].to_numpy()
trueper = kep_table['period_truth'].to_numpy()


##bls needs
# pgrid = bls_modular.period_grid(start=.3,stop=27,N=25000)
abc = np.logspace(0,1.45,15000)
next1 = np.linspace(0.3,1,10000)
pgrid = np.concatenate((next1,abc))

dgrid = bls_modular.duration_grid(start=.01,stop=.29,N=5)
KEBids=[];KEBtrueper=[];KEBpowers=[];KEBperiods=[];KEBdepths=[];KEBdurs=[];KEBtts=[];KEBexpectedrp=[]

count = 0
for kic in tqdm(kics):
	imp_t.sleep(1)
	# print('starting {} out of {}'.format(count,len(kics)))
	klc = bls_modular.kepEBopen(kic) #open data
	flat_ktime, flat_kflux, flat_kfluxerr = bls_modular.flatten_lc(klc.time,klc.flux,klc.flux_err) #flatten
	# print('starting bls')
	periodogram, kstats = bls_modular.bls(pgrid,dgrid,flat_ktime, flat_kflux, flat_kfluxerr)
	np.save('data/KeplerEBs/bls_pg_period_{}'.format(kic),periodogram.period)
	np.save('data/KeplerEBs/bls_pg_power_{}'.format(kic),periodogram.power)
	#kstats order ===[ppower, pperiod, pdepth, pduration, ptransittime] 
	KEBids.append(kic);KEBtrueper.append(trueper[count])
	count +=1
	KEBpowers.append(kstats[0]);KEBperiods.append(kstats[1]);KEBdepths.append(kstats[2]);
	KEBdurs.append(kstats[3]);KEBtts.append(kstats[4])
#make a table
print('making final table: KeplerEB_bls_stats2')
kep_table2 = Table([KEBids,KEBtrueper,KEBperiods,KEBpowers,KEBdepths,KEBdurs,KEBtts],
	names=('KIC','period_truth','period_bls','power_bls','depth_bls','dur_bls','tt_bls'))

ascii.write(kep_table2, 'data/KeplerEBs/KeplerEB_bls_stats2.csv',format='csv', overwrite=True)  