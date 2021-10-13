import clean_lcs

def test_locate_files():
    sector = 14
    testpath = '/Volumes/Seagate-stars/Sectors/Sector_{}/Sec{}_rawLCs/'.format(sector,sector)
    tic = 7582634
    tic2 = 7582633 
    fullpath = clean_lcs.locate_files(tic,path=testpath)
    fullpath2 = clean_lcs.locate_files(tic2,path=testpath)
    assert isinstance(fullpath,list)
    assert fullpath[0][-10:-1] == '-s_lc.fit'
    assert len(fullpath2) == 3

    
def test_sector_ordered_files():
    filepaths = ['/Volumes/Seagate-stars/Sectors/Sector_14/Sec14_rawLCs/tess2019198215352-s0014-0000000007582634-0150-s_lc.fits']
    sectors_list = [14]
    tic = 7582634
    lcfs, sectorspresent = clean_lcs.sector_ordered_files(filepaths,sectors_list,tic)
    assert sectorspresent[0] == 14
    assert isinstance(lcfs, list)
    
def test_clean_files():
    filepaths = ['/Volumes/Seagate-stars/Sectors/Sector_14/Sec14_rawLCs/tess2019198215352-s0014-0000000007582634-0150-s_lc.fits']
    sectors_list = [14]
    tic = 7582634
    lcfs, sectorspresent = clean_lcs.sector_ordered_files(filepaths,sectors_list,tic)
    cleaned = clean_lcs.clean_files(lcfs)
    assert str(cleaned) == 'TessLightCurve(TICID: 7582634)'
    
test_locate_files()
test_sector_ordered_files()
test_clean_files()