import bls_modular
import numpy as np
import matplotlib.pyplot as plt

def test_fake_data():
    test1 = bls_modular.fake_data(10,.1,0,20)
    test2 = bls_modular.fake_data()
    assert(len(test1[0])==len(test1[1])==len(test1[2]))
    assert(len(test2[0])==len(test2[1])==len(test2[2]))


def test_open_mylcs():
    tic1= 140801898; sector1=1
    t1,f1,fe1 = open_mylcs(tic1,sector1) #passes
    assert len(t1)==len(f1)==len(fe1) #passes
    tic2 = 1234; sector2 = 1
    print('Expecting print statement about tic 1234:')
    t2,f2,fe2 = open_mylcs(tic2,sector2)
    assert t2 == 'None'
    
def test_kepEBopen():
    kic1 = 123; kic2 = 5376552
    path2 = 'no/correct/path'
    test1 = kepEBopen(kic1,path=path2)
    test2 = kepEBopen(kic2)
    test3 = kepEBopen(kic1)
    assert test1 == "Verify file path - couldn't open data file"
    assert str(type(test2)) == "<class 'lightkurve.lightcurve.LightCurve'>"
    assert test3 == "Verify KIC - couldn't open data file with mypath"
    
    
def test_plot_lc():
    #generate fake data
    sample_size = 1000
    time = np.linspace(0,25,sample_size)
    unc = 1e-3
    flux = np.random.randn(len(time))*unc + 1
    flux_err = np.ones_like(flux)*unc
    ID2 = '2ND TEST: with ErrorBars';ID3 = 3,'RD TEST: with funny title';
    ID4 = '4TH TEST: with one constant error';fluxerr2 = 2
    
    bls_modular.plot_lc(time,flux)#1st plot-defaults
    plt.close() #to avoid 2nd empty figure afterward
    bls_modular.plot_lc(time,flux,flux_err,ID2)#2nd plot
    plt.close()
    bls_modular.plot_lc(time,flux,ID=ID3)#3rd plot
    plt.close()
    bls_modular.plot_lc(time,flux,ID=ID4,flux_err=fluxerr2)#4th plot
    plt.close()
    print('Expecting 4 plots with various differences from test_plot_lc() fcn')
    

def test_flatten_lc():
    #generate fake data
    sample_size = 1000
    time = np.linspace(0,25,sample_size)
    unc = 1e-3
    flux = np.random.randn(len(time))*unc + 1
    flux_err = np.ones_like(flux)*unc
    #run fcn test
    one,two,three = bls_modular.flatten_lc(time,flux,flux_err)
    assert(len(one)==len(two)==len(three)) #checks all arrays same length
    assert(one[0]==0.) #checks time

def test_period_grid():
    start1=20;stop1=.5;N1=10;start2=1;stop2=10;N2=1
    test1 = bls_modular.period_grid(start1,stop1,N1)
    test2 = bls_modular.period_grid(start2,stop2,N2)
    test3 = bls_modular.period_grid()
    assert(test1[0]>test1[-1])
    assert(len(test2)==1)
    assert(len(test3)==5000)
    assert(test3[0]<test3[-1])
    
def test_duration_grid():
    start1=20;stop1=.5;N1=100;start2=1;stop2=10;N2=1
    test1 = bls_modular.duration_grid(start1,stop1,N1)
    test2 = bls_modular.duration_grid(start2,stop2,N2)
    test3 = bls_modular.duration_grid()
    assert(test1[0]>test1[-1])
    assert(len(test1)==100)
    assert(len(test2)==1)
    assert(len(test3)==3)
    assert(test3[0]<test3[-1])


def test_bls():
    #generate fake data
    time,flux,flux_err = bls_modular.fake_data()
    #generate grids
    periodgrid1=bls_modular.period_grid();periodgrid2=bls_modular.period_grid(5,10,100)
    dur_grid1=bls_modular.duration_grid();dur_grid2=bls_modular.duration_grid(1,7,4)
    periodgrid3=bls_modular.period_grid(1.,10,100);dur_grid3=bls_modular.duration_grid(.001,.99,3)
    #run bls test
    test1 = bls_modular.bls(periodgrid1,dur_grid1,time,flux,flux_err) #defaults
    test2 = bls_modular.bls(periodgrid2,dur_grid2,time,flux,flux_err)# failure b/c pgrid<dgrid
    test3 = bls_modular.bls(periodgrid3,dur_grid3,time,flux,0.1) # success
    test4 = bls_modular.bls(periodgrid1,dur_grid1,time[5::],flux,flux_err)
#     print('one:',test1,'two:',test2,'three:',test3)
    assert len(test1)==2
    assert test2=='ERROR: Minimum Period must be Greater than Maximum Duration'
    assert len(test3)==2
    assert test4=='ERROR: time and flux arrays NOT same length'

test_fake_data()
test_open_mylcs()
test_kepEBopen()
test_plot_lc()
test_flatten_lc()
test_period_grid()
test_duration_grid()
test_bls()