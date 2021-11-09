   
import numpy as np
import starspot as ss
from starspot import sigma_clipping
import pandas as pd 
import lightkurve as lk
import os
import matplotlib.pyplot as plt
import LS_modular
from matplotlib.gridspec import GridSpec


#change 4 lines with ############# at end to customize

targets_1415 = pd.read_csv('/Volumes/Seagate-stars/PAPER_FINAL_FILES/target_lists/secs_14_15_targetlist.csv')['TIC'].to_numpy() ######################
targets_NCVZ = pd.read_csv('/Volumes/Seagate-stars/PAPER_FINAL_FILES/target_lists/NCVZ_targetlist.csv')['TIC'].to_numpy() ######################
targets_SCVZ = pd.read_csv('/Volumes/Seagate-stars/PAPER_FINAL_FILES/target_lists/SCVZ_targetlist.csv')['TIC'].to_numpy() ######################
all_CVZ = np.concatenate([targets_NCVZ,targets_SCVZ])

train_tics = np.load('/Volumes/Seagate-stars/PAPER_FINAL_FILES/all_training_tics.npy')
##############################################################################
########################## make choices here ##########################
##############################################################################

# #option-1
# stitched = True
# pathstart = 'stitched'
# target_list = all_CVZ #or 'targets_NCVZ' or 'targets_SCVZ'

# #option-2
# stitched = False 
# pathstart = 'mycvz'
# target_list = all_CVZ #or 'targets_NCVZ' or 'targets_SCVZ'

# #option-3
# stitched = False
# pathstart = 'secs1415'
# target_list = targets_1415

#or customize target list
stitched = False #not available for stitched light curves
pathstart = 'custom'
mymask = [i not in train_tics for i in targets_SCVZ]
SCVZ_notin_train = targets_SCVZ[mymask]
target_list =  SCVZ_notin_train #train_tics[282::]
##############################################################################

#to store calculated stats
tic_list = [] #to verify picking order stays same
pdms_err =[] 
pdms =[]


def runthis(target_list,stitched):
    
    for count,tic in enumerate(target_list):
        print('starting work on tic:',tic)
        print('starting work on tic:',tic)
        print('starting work on tic:',tic)
        print('starting work on tic:',tic)
        print('starting work on tic:',tic)
        print('         starting work on tic:',tic, ' and with count:  ',count,'of ',len(target_list))
        print('         starting work on tic:',tic, ' and with count:  ',count)
        
#determine paths to lcfiles
        if stitched == True:
            paths = LS_modular.pathfinder(tic,sector=pathstart)
        elif stitched == False:
            if pathstart == 'mycvz':
                paths = LS_modular.find_secfiles(tic,pathstart)
            elif pathstart == 'secs1415':
                paths = LS_modular.find_secfiles(tic,pathstart)
            elif pathstart == 'custom': #takes longer b/c needs determine path & files avaliable
                if np.isin(tic,all_CVZ)==True:
                    newpathstart = 'mycvz'
                elif np.isin(tic,targets_1415)==True:
                    newpathstart = 'secs1415'
                else:
                    print('Error could not determine path for TIC {} based on target lists available'.format(tic))
                    continue
                paths = LS_modular.find_secfiles(tic,newpathstart)
            else: #for use by other computers
                paths = LS_modular.find_secfiles(tic,pathstart)
        else:
            print('Use boolean to specifiy if want stitched light curves or individual sector light curves')
            break

#iterate through available paths per tic to open data
        for p in paths:
            lcf = LS_modular.open_lcf(p)
            lc = lcf.FLUX #make into a lc object to get data   
            sector = LS_modular.header_sector(lcf)
            flux = lc.flux
            flux_err = lc.flux_err
            time = lc.time
#creating model & gathering stats---same as code for stats used in ls_measure() from LS_modular.py
            rotate = ss.RotationModel(time, flux, flux_err)

#ls
            ls_period = rotate.ls_rotation(min_period=.01 , max_period=27.)
            power = rotate.power
            freq = rotate.freq
            ps = 1./freq
            peaks = np.array([i for i in range(1, len(ps)-1) if power[i-1] < \
                        power[i] and power[i+1] < power[i]])
            peak_amps_low2high = np.sort(power[peaks])
            second_rp = ps[power == peak_amps_low2high[-2]][0]
            third_rp = ps[power == peak_amps_low2high[-3]][0]

#folding lc
            folded_lc = lc.fold(period=ls_period,t0=time[0])
            folded_lc2 = lc.fold(period=second_rp,t0=time[0])

#plots
            #set up figure for plots
            fig = plt.figure(figsize=(15,20),tight_layout=True)
            gs1 = GridSpec(5, 2)
            ax1 = fig.add_subplot(gs1[0:2, :])
            ax2 = fig.add_subplot(gs1[2:4, :])
            ax3 = fig.add_subplot(gs1[4, :-1])
            ax4 = fig.add_subplot(gs1[4, -1])            
            plt.subplots_adjust(hspace=0.5,wspace=0.25)

            #light curve
            ax1.scatter(time,flux,color='k',s=.3,label='cleaned light curve')
            ax1.set_xlabel('Time [days]')
            ax1.set_ylabel('Relative Flux')
            ax1.set_title('Sector {} Light Curve for TIC:{}'.format(sector,tic),fontsize=30);
            ax1.legend(prop={'size': 12})

            #ls periodogram
            ax2.plot(-np.log10(freq), power, "k", zorder=0)
            ax2.axvline(np.log10(ls_period), color="C1", lw=4, alpha=0.5,
                            zorder=1,label=('{} days'.format("%.4f" % ls_period)))
            ax2.axvline(np.log10(second_rp),lw=4,alpha=0.5,zorder=2,linestyle='--',color='cyan',label=('{} days'.format("%.4f" % second_rp)))
            ax2.axvline(np.log10(third_rp),lw=4,alpha=0.5,zorder=3,linestyle=(0, (1, 10)),color='g',label=('{} days'.format("%.4f" % third_rp)))
            ax2.set_xlabel("log10(Period [days])")
            ax2.set_ylabel("Power");
            ax2.set_title('Lomb-Scargle Periodogram for TIC:{}'.format(tic),fontsize=30)
            ax2.legend(prop={'size': 15})

            #folded lc
            ax3.scatter(folded_lc.phase,folded_lc.flux,color='k',s=.2)
            ax3.set_xlabel("Phase")
            ax3.set_ylabel("Flux")
            ax3.set_title('Folded at LS-1 ({}) Period TIC:{}'.format("%.2f" % ls_period,tic),fontsize=30)

            #folded lc at 2nd peak
            ax4.scatter(folded_lc2.phase,folded_lc2.flux,color='k',s=.2)
            ax4.set_xlabel("Phase")
            ax4.set_ylabel("Flux")
            ax4.set_title('Folded at LS-2 ({}) Period TIC:{}'.format("%.2f" % second_rp,tic),fontsize=30)

#save figure
            savepath = '/Volumes/Seagate-stars/PAPER_FINAL_FILES/Plots/{}/sec{}_lc_ls_folded.png'.format(tic,sector)#######################
            os.makedirs(os.path.dirname(savepath), exist_ok=True) #creates folders if necessary
            plt.savefig(savepath)
            plt.close()
            print('ENDING  TIC:',tic)

            

        
#### actually run this file
if __name__=="__main__":
    runthis(target_list,stitched)


