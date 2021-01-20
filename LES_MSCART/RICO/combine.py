"""
2017/04/03: Eddited to handdle missing files.
"""
import netCDF4, sys
import numpy as np
from cpnLES_MSCARTlib import POLCARTdset 

dpath    =str(sys.argv[1]) #'results/b0p860/'
fnames   =[str(sys.argv[2])]#'DYCOMS2_dharma_008036_b0p860_yx_MSCART_SZA140_SAA030_VAA000plus_NPH1e4'
runs     =int(sys.argv[3])
pc_format=str(sys.argv[4])#LES/step/fractal
action   =str(sys.argv[5])
missing_cnt=0
invalid_cnt=0

if sys.version_info[0]>2:
    print('For version earlier than 3: check exception to RuntimeError')            
def read_mscart(file_name):
    global missing_cnt,missing_files,invalid_cnt,invalid_files
    print(file_name)
    try:
        MC = netCDF4.Dataset(file_name, 'r')
    except RuntimeError:
        MC = netCDF4.Dataset(rpl_file,'r')
        print('-------MISSING!!!------  '+str(missing_cnt))
        if missing_cnt==0:
            missing_files=[file_name]
        else:
            missing_files=np.append(missing_files,file_name)           
        missing_cnt=missing_cnt+1       

    nbat = len(MC.dimensions[u'nbat'])
    Time = MC.variables['MeanTiming'][:]
    Mean = MC.variables['MeanPRad'][:]
    RMSE = MC.variables['RMSEPRad'][:]

    MC.close()
    if np.sum(Mean[:,:,:,0]>10)>0:
        print('Invalid high radiance value found in '+file_name+' !!')
        if invalid_cnt==0:
            invalid_files=[file_name]
        else:
            invalid_files=np.append(invalid_files,file_name)
        invalid_cnt+=1

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



glb_nbat=[]
for i in np.arange(0,np.size(fnames)):
    rpl_file=dpath+fnames[i]+'_2.nc'
    missing_cnt=0
    combine_mscart(fnames[i])
    print(str(missing_cnt)+' of missing files were replaced by '+rpl_file)

#saving as hdf5

if invalid_cnt==0:
    print('No files with invalid radiances!!')
    f_name=fnames[0]
    case=POLCARTdset('LES','/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/')
    if len(case.VAA)==2:
        #Principal plane. Combine VAA 0 and 180 together
        case.readMSCARTplus(f_name+'.nc',fdpath=dpath,step=True)
        case.savePOLCARThdf5(case.fdpath+case.fname.split('.',1)[0]+'.hdf5',pc_format=pc_format,action=action)
    else:
        #Has more 2 VAAs.
        case.cc3D='mulVAA'
        case.readMSCARTmulVAA(f_name+'.nc',fdpath=dpath)
        case.savePOLCARThdf5(case.fdpath+case.fname.split('.',1)[0]+'.hdf5',pc_format=pc_format,action=action)
else:
    print('Following files have invalid radiances!!')
    print(invalid_files)

print('Invalid job numbers:')
if invalid_cnt>0 and missing_cnt>0:
    missing_job=np.zeros(missing_files.size,dtype=int)
    invalid_job=np.zeros(invalid_files.size,dtype=int)
    missing_files=missing_files.reshape(missing_files.size)
    invalid_files=invalid_files.reshape(invalid_files.size)
    for i in np.arange(0,missing_files.size):
        missing_job[i]=int(missing_files[i].split('.',1)[0].rsplit('_',1)[1])
    for i in np.arange(0,invalid_files.size):
        invalid_job[i]=int(invalid_files[i].split('.',1)[0].rsplit('_',1)[1])
    redos=np.append(missing_job,invalid_job)
    uredos=np.unique(redos)
    uredos.sort()
    print("%d missing runs found!!"%(np.size(uredos)))
    print(uredos)
    
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
