#!/usr/bin/env python
"""
2017/04/03: Eddited to handdle missing files.
"""
import netCDF4
import numpy as np
from cpnLES_MSCARTlib import POLCARTdset 

missing_cnt=0

def read_mscart(file_name):
    global missing_cnt
    print(file_name)
    try:
        MC = netCDF4.Dataset(file_name, 'r')
    except FileNotFoundError:
        MC = netCDF4.Dataset(rpl_file,'r')
        missing_cnt=missing_cnt+1
        print('-------MISSING!!!------  '+str(missing_cnt))
       
           

    nbat = len(MC.dimensions[u'nbat'])
    Time = MC.variables['MeanTiming'][:]
    Mean = MC.variables['MeanPRad'][:]
    RMSE = MC.variables['RMSEPRad'][:]

    MC.close()

    return (nbat, Time, Mean, RMSE)

def sts_mscart(file_prefix, nexp):
   global glb_nbat
   (nbat, Time, Mean, RMSE) = read_mscart(file_prefix+'_1.nc')
   (nsrc, npr, npy, npx, nstokes) = Mean.shape

   CMean = np.zeros((nsrc, npr, npy, npx, nstokes))
   CRMSE = np.zeros((nsrc, npr, npy, npx, nstokes))
   CTime = 0.0

   for i in range(nexp):
      (nbat, Time, Mean, RMSE) = read_mscart(file_prefix+'_'+str(i+1)+'.nc')
      CMean = CMean + Mean
      CRMSE = CRMSE + RMSE**2 * (nbat - 1) + Mean**2
      CTime = CTime + Time * nbat

   CMean = CMean / nexp
   CRMSE = np.sqrt((CRMSE / nexp - CMean**2) / (nexp * nbat - 1))
   CTime= CTime/nexp
   glb_nbat=nbat

   return (CTime, CMean, CRMSE)

def write_mscart(file_name, Time, Mean, RMSE):
   (nsrc, npr, npy, npx, nstokes) = Mean.shape

   MC = netCDF4.Dataset(file_name, 'w')
   MC.createDimension('nsrc', nsrc)
   MC.createDimension('npr' , npr )
   MC.createDimension('npy' , npy )
   MC.createDimension('npx' , npx )
   MC.createDimension('nstokes', nstokes)
   
   MC.createVariable('MeanTiming', np.float32) 
   MC.createVariable('MeanPRad', np.float32, ('nsrc','npr','npy','npx','nstokes',))
   MC.createVariable('RMSEPRad', np.float32, ('nsrc','npr','npy','npx','nstokes',))
   MC.createVariable('Nbat',np.int16)
   
   MC.variables['MeanTiming'][:] = Time
   MC.variables['MeanPRad'][:] = Mean
   MC.variables['RMSEPRad'][:] = RMSE
   MC.variables['Nbat'][:]=runs*glb_nbat
   MC.close()

   return None

def combine_mscart(file_prefix):
   (Time, Mean, RMSE) = sts_mscart(dpath+file_prefix, runs)
   write_mscart(dpath+file_prefix+'.nc', Time, Mean, RMSE)
   print(dpath+file_prefix+'.nc Saved!! ')
   return None

dpath='results/b2p13/'
#fnames=['OP_dharma_008036_full_3_26_MSCART_SSA020_SAA030_VAA000_NPH2e6',\
#        'OP_dharma_008036_full_3_26_MSCART_SSA020_SAA030_VAA180_NPH2e6']
fnames=['fractal_cld_b2p13re12ve05_x40km_MSCART_SZA100_SAA000_VAA000plus_NPH1e5']
runs=100
glb_nbat=[]
for i in np.arange(0,np.size(fnames)):
    rpl_file=dpath+fnames[i]+'_3.nc'
    missing_cnt=0
    combine_mscart(fnames[i])
    print(str(missing_cnt)+' of missing files were replaced by '+rpl_file)

#saving as hdf5
f_name=fnames[0]
fracSZAx=POLCARTdset('fracSZAx','/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Toy_clouds/Fractal_Cloud/')
fracSZAx.readMSCARTplus(f_name+'.nc',fdpath=dpath,step=True)
fracSZAx.savePOLCARThdf5(fracSZAx.fdpath+fracSZAx.fname.split('.',1)[0]+'.hdf5',pc_format='fractal_cloud',action='spot3D_revision')


#combine_mscart('dharma_008036_full_MSCART_SZA20_SAA30_V180_3_26')

#combine_mscart('C3_1_FMC_0')
#combine_mscart('C3_1_FMC_1')
#combine_mscart('C3_1_FMC_2')
#combine_mscart('C3_1_FMC_3')
#combine_mscart('C3_1_FMC_4')
#combine_mscart('C3_1_FMC_5')
#combine_mscart('C3_1_FMC_6')
#combine_mscart('C3_1_FMC_7')
#combine_mscart('C3_1_FMC_8')
#
#combine_mscart('C3_2_FMC_0')
#combine_mscart('C3_2_FMC_1')
#combine_mscart('C3_2_FMC_2')
#combine_mscart('C3_2_FMC_3')
#combine_mscart('C3_2_FMC_4')
#combine_mscart('C3_2_FMC_5')
#combine_mscart('C3_2_FMC_6')
#combine_mscart('C3_2_FMC_7')
#combine_mscart('C3_2_FMC_8')
