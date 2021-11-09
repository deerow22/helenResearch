#this should work to create a curl for just the additions


import lightkurve as lk
import Useful_functions as uf
import pandas as pd
import numpy as np

#open data to get additional targets
curlsavepath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/shell_scripts/dl_2021_additions.sh'
old_master_df = pd.read_csv('/Volumes/Seagate-stars/Final_Run/test_sets/master_df.csv')
updated_targets = pd.read_csv('/Volumes/Seagate-stars/PAPER_FINAL_FILES/target_lists/all_final_targets.csv')
print('old master df length:',len(old_master_df))
a = np.unique(old_master_df['TIC'].to_numpy(int))
b = np.unique(updated_targets['ID'].to_numpy(int))
# print(type(a),type(b),type(a[0]),type(b[0]))
old_targets_missing, new_targets_toadd = uf.returnNotMatches(a, b)
print(len(new_targets_toadd), 'targets added in updated lists')


#search for available files and store info into curl string
curlscript = []
for num,tic in enumerate(new_targets_toadd):
    print('~~~starting {} out of {}'.format(num,len(new_targets_toadd)))
    ticid = 'TIC {}'.format(tic)
    lc_search = lk.search_lightcurvefile(ticid)
    mast_table = lc_search.table
    sectors = list(mast_table['sequence_number']) #int
    idx = mast_table['#'] #object
    for count,sec in enumerate(sectors):
        if sec <=26:
            i = int(idx[count])
            firststr = 'curl -C - -L -o '
            middlestr = str(mast_table['obs_id'][i]) +'_lc.fits '
            webaddy = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=" 
            laststr = webaddy +str(mast_table['dataURL'][i])
            script = firststr + middlestr + laststr 
            curlscript.append(script)
            print('finished sector {}'.format(sec))
        else:
            pass
curlscript=np.array(curlscript)
print('final curlscript has shape:',curlscript.shape)


#write curl string to script

with open ('{}'.format(curlsavepath), 'w') as rsh:
    rsh.write('''\
#! /bin/bash''')
    for count,script in enumerate(curlscript):
        rsh.write('''
{}'''.format(script))
print('F I N I S H E D')

