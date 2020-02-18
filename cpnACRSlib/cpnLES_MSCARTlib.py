#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri Apr  7 16:01:02 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
python library for LES_MSCART

05/24/2017:
    collectColsMSCART : To collect C-by-C MSCART runs and create single POLCARTdset 
    object that contains all bins. 
05/25/2017:
    LES_filed object to handle MSCART "Field files". Will be usefull when change
    the physical properties of the scene (like adding aerosols etc.)
    0 and 180 VAA. Negative VZA => 180.0 Positive VZA=>0
06/29/2017:
    Change iqu and iquRMS to handly 0 VZA.
    According to the step cloud case, I think we should not swap the xy axises.
    Thus, removed all "transpose" parts both for column and 3d readings.
08/17/2017:
    savePOLCARThdf5 updated to handle 1D data (for the step cloud case)
08/21/2017:
    Going to make SZA: Source zenith angle and SAA=Source Azimuth Angle
    Source=>"The direction of the solar radiaton" not the "Sun"
    Added scat_ang function.
08/22/2017:
    readMSCARTplus adapted to handle 1D step cloud (column) runs.
    savePOLCARThdf5 also adapted to handle cc3D='step_3D' and cc3D='step_1D'(not tested yet)
10/17/2017:
    avg_resolution() added to average over different resolutions.   
11/20/2017:
    An exception handled when "ScatA" not exist in the hdf5 data.
11/23/2017:
    xy_extent added to iqu() and ip(). 
    pad=0.2 added to iqu() and ip() to avoid overlaping colorbars and x labels.
    set_ScatA() function added.
06/04/2018:
    DHARMA_onmp added. (To convert DHARMA LES outputs to MSCART input field files)    
03/05/2019:
    find_reVW() added. (bidirectional weighting function Platnick 2000 )
    
"""
import numpy as np
import netCDF4
import h5py, os, sys
#from pyhdf.SD import SD
import time, copy
import matplotlib.pyplot as plt
from textwrap import wrap
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from cpnCommonlib import movingaverage2D
from cpnMielib import MieSet

def add_cloud_mask(dict_atr,cldM,dictionary=False):
    '''
    A genralized function to add a cloud mask.
    Be careful about the content. All the attributes that satisfy certain condtions 
    will be masked. So make sure there are no any array that SHOULD NOT BE MASKED!!
    -----------------------------------------------------------------------------------
    dict_atr + dictionary=False: dict_atr must be an object. Will search over all attributes.
    dict_atr + dictionary=True : dict_atr must be a dictionary
    Find all the ndarrays and mask,
        (1) n=2 arrays [cldM]=np.nan
        (2) last two dimensions of n=3 arrays [:,cldM]=np.nan 
    '''
    print('Applying cloud mask...')
    arrs={'two':[],'lastPair':[],'firstPair':[]}
    if dictionary:
        variables=dict_atr
    else:
        variables=vars(dict_atr)
    for key in variables.keys():
        if type(variables[key]) is np.ndarray:
            if variables[key].ndim==2:
                arrs['two']      =np.append(arrs['two'],key)
            if variables[key].ndim==3 and variables[key].shape[0]==variables[key].shape[1]:
                arrs['firstPair'] =np.append(arrs['firstPair'],key)
            if variables[key].ndim==3 and variables[key].shape[1]==variables[key].shape[2]:
                arrs['lastPair']=np.append(arrs['lastPair'],key)
    print('Masking two dimensions:')
    for key in arrs['two']:
        print('\t'+key)
        variables[key][np.invert(cldM)]=np.nan
    print('Masking last two dimensions:')
    for key in arrs['lastPair']:
        print('\t'+key)
        variables[key][:,np.invert(cldM)]=np.nan
    print('Masking first two dimensions:')
    for key in arrs['firstPair']:
        print('\t'+key)
        variables[key][np.invert(cldM),:]=np.nan
class POLCARTdset(object):
    def __init__(self,dset,nmldpath):
        self.nmldpath=nmldpath
        self.dset=dset
        self.cc3D=[]#1D 3D or step
        self.SZA=[]#source zenith angle
        self.SAA=[]#source azimuth angle
        self.VZA=[]
        self.VAA=[]
        self.NPH=[]
        self.Nbat=[]
        self.MeanPRad=[]
        self.MeanTiming=[]
        self.RMSEPRad=[]
        self.pnt=[]
        

    def readMSCARTplus(self,fname,fdpath=None,bm=False,clm=False,prnt=True,\
                       step=False):
        '''
        fname= 'filename.nc'
        bm=True for benchmakrk data from Dr. Wang
        clm = True for single column data
        step = True if the cloud field span along x axis (like the step cloud)
        plus means observation along complete veiwing plane (2 azimuth angles)
        '''
        if fdpath!=None:
            data=netCDF4.Dataset(fdpath+fname,'r')
            self.fdpath=fdpath
        else:
            data=netCDF4.Dataset(fname,'r')
        self.fname=fname
        self.MeanTiming=np.squeeze(data.variables['MeanTiming'][:])
        MeanPRad=np.squeeze(data.variables['MeanPRad'][:])
        RMSEPRad=np.squeeze(data.variables['RMSEPRad'][:])
        if bm:
            self.Nbat=500*10 # number of batches for benchmark case
        else:
            try:
                self.Nbat=data.variables['Nbat'][:]
            except KeyError:
                if prnt:
                    print(self.fname+' is a single batch file')
                self.Nbat=1
#        self.NPH=fname[-6:-3]
        self.NPH=(fname.split('_NPH',1)[1]).split('.',1)[0]

#        OP_dharma_008036_full_3_26_MSCART_SSA20_SAA30_VAA0.nml
#        fr=open(self.nmldpath+fname[:-10]+'.nml','r')
        fr=open(self.nmldpath+fname.split('_NPH',1)[0]+'.nml','r')
        d2=fr.readlines()[:]
        fr.close()
        self.SZA=float(d2[42].split('= ',1)[1].split(',',1)[0])
        self.SAA=float(d2[43].split('= ',1)[1].split(',',1)[0])
        n_view_the=int(d2[48].split('= ',1)[1].split(',',1)[0])
        n_view_phi=int(d2[49].split('= ',1)[1].split(',',1)[0])
        mn_view_the=float(d2[50].split('= ',1)[1].split(',',1)[0])
        mx_view_the=float(d2[51].split('= ',1)[1].split(',',1)[0])
        vamin=float(d2[52].split('= ',1)[1].split(',',1)[0])
        vamax=float(d2[53].split('= ',1)[1].split(',',1)[0])
        self.VAA=np.linspace(vamin,vamax,n_view_phi)
        self.VZA=np.hstack((np.flipud(np.linspace(mn_view_the,mx_view_the,n_view_the)*np.cos(np.deg2rad(self.VAA[1]))),np.linspace(mn_view_the,mx_view_the,n_view_the)*np.cos(np.deg2rad(self.VAA[0]))))
        self.ScatA=np.hstack((np.flipud(scat_ang(self.SZA,np.linspace(mn_view_the,mx_view_the,n_view_the),self.SAA-self.VAA[1])),\
                              scat_ang(self.SZA,np.linspace(mn_view_the,mx_view_the,n_view_the),self.SAA-self.VAA[0])))
        if clm :
            MeanPRad=np.vstack((np.flipud(MeanPRad[1::2,:]),MeanPRad[0::2,:]))
            RMSEPRad2=np.vstack((np.flipud(RMSEPRad[1::2,:]),RMSEPRad[0::2,:]))
#            self.MeanPRad=np.transpose(MeanPRad,(1,0))
#            self.RMSEPRad=np.transpose(RMSEPRad2,(1,0))
            self.MeanPRad=np.transpose(MeanPRad,(0,1))
            self.RMSEPRad=np.transpose(RMSEPRad2,(0,1))
            if step:
                self.cc3D='step_1D'
                self.MeanPRad=np.transpose(MeanPRad,(0,1))
                self.RMSEPRad=np.transpose(RMSEPRad2,(0,1))
            else:
                self.cc3D="1D"
                self.MeanPRad=np.transpose(MeanPRad,(0,1))
                self.RMSEPRad=np.transpose(RMSEPRad2,(0,1))
                self.pnt=np.zeros(2,dtype='int')
                self.pnt[1]=int(self.fname.split('binsy',1)[1].split('x',1)[0])
                self.pnt[0]=int(self.fname.split('binsy',1)[1].split('x',1)[1].split('_',1)[0])
                #MeanPRAd dimenstions
        else:
            if step:
                self.MeanPRad=np.vstack((np.flipud(MeanPRad[1::2,:,:]),MeanPRad[0::2,:,:]))
                self.RMSEPRad=np.vstack((np.flipud(RMSEPRad[1::2,:,:]),RMSEPRad[0::2,:,:]))
                self.cc3D="step_3D"
            else:
                MeanPRad2=np.vstack((np.flipud(MeanPRad[1::2,:,:,:]),MeanPRad[0::2,:,:,:]))
                RMSEPRad2=np.vstack((np.flipud(RMSEPRad[1::2,:,:,:]),RMSEPRad[0::2,:,:,:]))
                #MeanPRAd dimenstions
#                n_v_the, npy, npx, nstokes=MeanPRad.shape
#                self.MeanPRad=np.transpose(MeanPRad2,(0,2,1,3))
#                self.RMSEPRad=np.transpose(RMSEPRad2,(0,2,1,3))
                self.MeanPRad=np.transpose(MeanPRad2,(0,1,2,3))
                self.RMSEPRad=np.transpose(RMSEPRad2,(0,1,2,3))
                self.cc3D="3D"
        data.close()
            #number of sources, num of viewing polar angles, 
    def readMSCART(self,fname,fdpath=None,bm=False,clm=False):
        #read MSCART benchmark results (from Dr. Wang)
        #clm : True if a single column.
        if fdpath!=None:
            data=netCDF4.Dataset(fdpath+fname,'r')
            self.fdpath=fdpath
        else:
            data=netCDF4.Dataset(fname,'r')
        self.fname=fname
        self.MeanTiming=np.squeeze(data.variables['MeanTiming'][:])
        MeanPRad=np.squeeze(data.variables['MeanPRad'][:])
        RMSEPRad=np.squeeze(data.variables['RMSEPRad'][:])
        if bm:
            self.Nbat=500*10 # number of batches for benchmark case
        else:
            self.Nbat=data.variables['Nbat'][:]
        self.NPH=fname[-6:-3]

#        OP_dharma_008036_full_3_26_MSCART_SSA20_SAA30_VAA0.nml
        try:
            print('Trying to read '+self.nmldpath+fname[:-10]+'.nml')
            fr=open(self.nmldpath+fname[:-10]+'.nml','r')
        except IOError:
            print(self.nmldpath+fname[:-10]+'.nml NOT FOUND!!')
            print('Trying '+self.nmldpath+fname)
            fr=open(self.nmldpath+fname+'.nml','r')
        d2=fr.readlines()[:]
        fr.close()
        self.SZA=float(d2[42].split('= ',1)[1].split(',',1)[0])
        self.SAA=float(d2[43].split('= ',1)[1].split(',',1)[0])
        n_view_the=int(d2[48].split('= ',1)[1].split(',',1)[0])
        mn_view_the=float(d2[50].split('= ',1)[1].split(',',1)[0])
        mx_view_the=float(d2[51].split('= ',1)[1].split(',',1)[0])
        self.VAA=float(d2[52].split('= ',1)[1].split(',',1)[0])
        if self.VAA==0 or self.VAA==180:
            self.VZA=np.linspace(mn_view_the,mx_view_the,n_view_the)*np.cos(np.deg2rad(self.VAA))
            self.ScatA=scat_ang(self.SZA,np.linspace(mn_view_the,mx_view_the,n_view_the),self.SAA-self.VAA)
        else:
            self.VZA=np.linspace(mn_view_the,mx_view_the,n_view_the)
            
        #MeanPRAd dimenstions
        if not(clm):
            n_v_the, npy, npx, nstokes=MeanPRad.shape
            self.MeanPRad=np.transpose(MeanPRad,(0,2,1,3))
            self.RMSEPRad=np.transpose(RMSEPRad,(0,2,1,3))
        else:
            self.MeanPRad=MeanPRad
            self.RMSEPRad=RMSEPRad
        data.close()
        #number of sources, num of viewing polar angles, 
    def readMCPOL(self,fname,fdpath=None):
        #read MCPOL results from Dans
        self.MeanTiming=0
        self.Nbat=0
        self.fname=fname
        if fdpath==None:
            f=h5py.File(fname,'r')
        else:
            self.fdpath=fdpath
            f=h5py.File(fdpath+fname,'r')
        I=f['I'][:]/100
#        Ip=f['Ip'][:]
        Q=f['Q'][:]/100
        U=f['U'][:]/100
#        V=f['V'][:]
        self.VAA=np.squeeze(f['VAA'][:])
        self.VZA=np.squeeze(f['VZA'][:])
        dI=f['dI'][:]/100
#        dIp=f['dIp'][:]
        dQ=f['dQ'][:]/100
        dU=f['dU'][:]/100
#        dV=f['dV'][:]
        self.scat=f['scat'][:]
        Mrad=np.squeeze([[I],[Q],[U]])
        nstokes,npx,npy,n_v_the=Mrad.shape
        self.MeanPRad=np.transpose(Mrad,(3,2,1,0))
        self.RMSEPRad=np.transpose(np.squeeze([[dI],[dQ],[dU]]),(3,2,1,0))
        
    def savePOLCARThdf5(self,out_name,dpath=None,pc_format=None,action=None):
        '''
        out_name=output file name fname_band.hdf5       
        '''
        if dpath==None:
            if os.path.isfile(out_name):
                print(out_name+' File already exist! Check!')
            else:
                f=h5py.File(out_name,'w')
        else:
            if os.path.isfile(dpath+out_name):
                print(dpath+out_name+' file already exist! Check!')
            else:
                f=h5py.File(dpath+out_name,'w')
        if pc_format==None:
            f.attrs['POLCART_format']='unknown'#MSCART_1D,MCPOL,MSCART_3D
        else:
            f.attrs['POLCART_format']=pc_format
        if action==None:
            f.attrs['Action']='unknown'#'To_tackle_parallax'
        else:
            f.attrs['Action']=action
        
        if self.cc3D=='step_3D' or self.cc3D=='step_1D':
            PCentry=f.create_dataset('MeanPRad',data=self.MeanPRad)
            PCentry.dims[0].label='VZA'
            PCentry.dims[1].label='xgrid'
            PCentry.dims[2].label='IQUV'
            PCentry.attrs['units']='Radiance'
            PCentry.attrs["long_name"]='Mean_Polarized_Radiance'
            
            PC=f.create_dataset('RMSEPRad',data=self.RMSEPRad)
            PC.dims[0].label='VZA'
            PC.dims[1].label='xgrid'
            PC.dims[2].label='IQUV'
            PC.attrs['units']='Radiance'
            PC.attrs["long_name"]='RMSE_Polarized_Radiance'
        else:
            PCentry=f.create_dataset('MeanPRad',data=self.MeanPRad)
            PCentry.dims[0].label='VZA'
            PCentry.dims[1].label='xgrid'
            PCentry.dims[2].label='ygrid'
            PCentry.dims[3].label='IQUV'
            PCentry.attrs['units']='Radiance'
            PCentry.attrs["long_name"]='Mean_Polarized_Radiance'

            PC=f.create_dataset('RMSEPRad',data=self.RMSEPRad)
            PC.dims[0].label='VZA'
            PC.dims[1].label='xgrid'
            PC.dims[2].label='ygrid'
            PC.dims[3].label='IQUV'
            PC.attrs['units']='Radiance'
            PC.attrs["long_name"]='RMSE_Polarized_Radiance'

        PC=f.create_dataset('SAA',data=np.array([self.SAA]))
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='Source(ray)_Azimuth_Angle'
        
        PC=f.create_dataset('SZA',data=np.array([self.SZA]))
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='Source(ray)_Zenith_Angle'

        PC=f.create_dataset('VZA', data=self.VZA)
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='Viewing_Zenith_Angle_negative_VAA180'

        PC=f.create_dataset('VAA',data=self.VAA)
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='Viewing_Azimuth_Angle'
        
        PC=f.create_dataset('ScatA',data=self.ScatA)
        PC.attrs['units']='degrees'
        PC.attrs['long_name']='(Apparent)Scattering_Angle'

        PC=f.create_dataset('MeanTiming',data=np.array([self.MeanTiming]))
        PC.attrs['units']='mins'
        PC.attrs['long_name']='Mean_time_per_batch'
        if sys.version_info[0]>2:
            PC=f.create_dataset('NPH',data=np.string_([self.NPH]))
        else:
            PC=f.create_dataset('NPH',data=np.array[self.NPH])
        PC.attrs['units']='counts'
        PC.attrs['long_name']='Num_of_photons_per_batch'

        PC=f.create_dataset('Nbat',data=np.array([self.Nbat]))
        PC.attrs['units']='counts'
        PC.attrs['long_name']='Number_of_batches'

        PC=f.create_dataset('File_stamp/cc3D',data=self.cc3D)
        PC.attrs['units']='None'
        PC.attrs['long_name']='Column_by_column_or_3D'

        PC=f.create_dataset('File_stamp/dset',data=self.dset)
        PC.attrs['units']='None'
        PC.attrs['long_name']='POLCARTdset_obj_name'

        PC=f.create_dataset('File_stamp/fname',data=self.fname)
        PC.attrs['units']='None'
        PC.attrs['long_name']='Origina_data_file_name'

        PC=f.create_dataset('File_stamp/nmldpath',data=self.nmldpath)
        PC.attrs['units']='None'
        PC.attrs['long_name']='Original_data_namelist_path'
    
        f.close()
        print(out_name+' Saved!!')
        
    def readPOLCARThdf5(self,filename,dpath=None):
        if dpath==None:
            f=h5py.File(filename,'r')
        else:
            f=h5py.File(dpath+filename,'r')
            self.hdfpath=dpath
        self.cc3D=f['File_stamp/cc3D'].value
        self.dset=f['File_stamp/dset'].value
        self.fname=f['File_stamp/fname'].value
        self.nmldpath=f['File_stamp/nmldpath'].value
        self.MeanPRad=f['MeanPRad'][:]
        self.MeanTiming=f['MeanTiming'][0]
        self.NPH=f['NPH'][0]
        self.Nbat=f['Nbat'][0]
        self.RMSEPRad=f['RMSEPRad'][:]
        self.SAA=f['SAA'][0]
        self.SZA=f['SZA'][0]
        self.VAA=f['VAA'][:]
        self.VZA=f['VZA'][:]
        try:
            self.ScatA=f['ScatA'][:]
        except KeyError:
            print('WARNING!!!!!! NO "ScatA" in '+filename)
            print('Check SZA,VZA,SAA,VAA and then run set_ScatA() to get appropriate scattering angles.')
        if type(self.fname) is np.ndarray:
            self.fname=self.fname[0].astype(str)
        f.close()
    def find_obVZA_ix(self,VZA):
        #Return the index that corresponds to the give VZA
        if VZA==-20.0:
            return 40
        elif VZA==0.0:
            return 61
        else:
            ix=np.where(np.isclose(self.VZA,VZA))
            print('%%%VZA indices for observations%%%%')
            print('VZA indices:'+str(ix))
            print('VZA: '+str(self.VZA[ix]))
            obVZA_ix=input('Select VZA index:')
            return obVZA_ix
    def avg_resolution(self,ftr,RT_dim='1D'):
        '''
        Averaging radiance over given factor. (Not a moving average)
        1D field
        ---------
            ex. If we have radiance field with 1000 x elements, if ftr=5 we'll have
            a new radiance field with 1000/5 x. Every 5 elements will be averaged.
        2D field
        ---------
            ex. If we have a 2D field, 2D averaging.(Not implemeted yet)
        obj: POLCARTdset object
        ftr: factor to do averaging.
        '''
        shape=self.MeanPRad.shape
        newMPRad=np.zeros((shape[0],shape[1]/ftr,shape[2]))
        for i in range(0,self.MeanPRad.shape[0]):
            MPRad=np.squeeze(self.MeanPRad[i,:,:])
            if np.remainder(shape[1],ftr)==0:
                B=MPRad.reshape(-1,ftr,shape[2])
            else:
                rmv_ix=MPRad.shape[0]-MPRad.shape[0]/ftr*ftr
                B=MPRad[0:-rmv_ix].reshape(-1,ftr,shape[2])
            newMPRad[i,:,:]=np.mean(B,1)
        return newMPRad
    def movingavg_res(self,ftr,RT_dim='2D'):
        '''
        Find the moving average of the 2D field based on the given factor ftr
        ftr: factor to average (the length of the one side of the averaging domain)
        RT_dim: '2D' or '1D'
        '''
        newMPRad=np.zeros_like(self.MeanPRad,dtype=float)
        for i in np.arange(0,self.VZA.size):
            for l in np.arange(0,4):
                newMPRad[i,:,:,l]=movingaverage2D(self.MeanPRad[i,:,:,l],window=ftr)
        return newMPRad
    def set_ScatA(self,):
        mn_view_the=np.min(np.abs(self.VZA))
        mx_view_the=np.max(np.abs(self.VZA))
        n_view_the=np.size(self.VZA)/2
        self.ScatA=np.hstack((np.flipud(scat_ang(self.SZA,np.linspace(mn_view_the,mx_view_the,n_view_the),self.SAA-self.VAA[1])),\
                    scat_ang(self.SZA,np.linspace(mn_view_the,mx_view_the,n_view_the),self.SAA-self.VAA[0])))
  
    def put_scat_ang(self,ax2):
        '''
        Put the corresponding scattering angle ticks and label on the top of the plot
        '''
        ax2t=ax2.twiny()
        ax2ticks=ax2.get_xticks()
        ax2tticks=ax2ticks
        ax2t.set_xticks(ax2tticks)
        ax2t.set_xbound(ax2.get_xbound())
        ax2t.set_xticklabels(tick_fun(ax2tticks,self.SZA,self.SAA))
        ax2.set_xlabel('VZA',size=10)
        ax2t.set_xlabel('Apparent Scat. Angle',size=10)
    def find_DolP(self,iqu=None):
        '''
        Calculating degree of linear polarization
        return:
            DolP,lP
        '''
        self.lP=np.sqrt(self.MeanPRad[:,1]**2+self.MeanPRad[:,2]**2)
        self.DolP=self.lP/self.MeanPRad[:,0]
        if not(iqu==None):
            lP=np.sqrt(iqu[1]**2+iqu[2]**2)
            DolP=lP/iqu[0]
            return DolP,lP
    def remove_redundant_nadir(self,removNadir_i=60):
        '''
        When complete principle-plane gemetries are considered, it would be 
        two sets of VZA (0 to 60) for each VAA=0 and 180. So the nadir geometry
        is considered twice. 
        This is to remove columns correspond to the additional angle zero.
        '''
        self.VZA     =np.delete(self.VZA,removNadir_i)
        self.MeanPRad=np.delete(self.MeanPRad,removNadir_i,0)
        self.RMSEPRad=np.delete(self.RMSEPRad,removNadir_i,0)
        if 'ScatA' in vars(self):
            self.ScatA=np.delete(self.ScatA,removNadir_i)            
    def read_from_nml(self,ix):
        '''
        Reading values from *.nml files. Confirm indexes when uses!!!
        54: vec_PRad_zlv(1)
        '''
        fl = open(self.nmldpath+self.fname.split('_NPH',1)[0]+'.nml')
        data = fl.readlines()[:]
        fl.close()
        line = [s.split(',\n',1)[0] for s in data][ix]
        print(line)
        return float(line.split('=',1)[1].replace(' ',''))        
class LES_field(object):
    def __init__(self,fname,dpath=None):
        self.fname=fname
        self.dpath=dpath
    def readLES_field(self,):
        if self.dpath==None:
            data=netCDF4.Dataset(self.fname,'r')
        else:
            data=netCDF4.Dataset(self.dpath+self.fname,'r')
        self.file_format=data.file_format
        self.nxp=len(data.dimensions[u'nxp'])
        self.nyp=len(data.dimensions[u'nyp'])
        self.nzp=len(data.dimensions[ u'nzp'])
        self.nx=len(data.dimensions[u'nx'])
        self.ny=len(data.dimensions[u'ny'])
        self.nz=len(data.dimensions[u'nz'])
        self.nz3=len(data.dimensions[u'nz3'])
        self.iz3l=len(data.dimensions[u'iz3l'])
        self.npar1=len(data.dimensions[u'npar1'])
        self.npar3=len(data.dimensions[u'npar3'])
        self.nkd=len(data.dimensions[u'nkd'])
        self.nstokes=len(data.dimensions['nstokes'])
        self.nang=len(data.dimensions[u'nang'])
        self.npf=len(data.dimensions[u'npf'])
        self.npsfc=len(data.dimensions[u'npsfc'])

        self.xgrd=data.variables['xgrd'][:]
        self.ygrd=data.variables['ygrd'][:]
        self.zgrd=data.variables['zgrd'][:]
        self.extp1d=data.variables['extp1d'][:]
        self.omgp1d=data.variables['omgp1d'][:]
        self.jpfp1d=data.variables['jpfp1d'][:]
        self.rpfp1d=data.variables['rpfp1d'][:]
        self.extp3d=data.variables['extp3d'][:]
        self.omgp3d=data.variables['omgp3d'][:]
        self.jpfp3d=data.variables['jpfp3d'][:]
        self.rpfp3d=data.variables['rpfp3d'][:]
        self.absg1d=data.variables['absg1d'][:]
        self.wkd=data.variables['wkd'][:]
        self.scatAng=data.variables['scatAng'][:]
        self.PMatVal=data.variables['PMatVal'][:]
        self.tsfc2d=data.variables['tsfc2d'][:]
        self.jsfc2d=data.variables['jsfc2d'][:]
        self.psfc2d=data.variables['psfc2d'][:]
        data.close()
    def plot_PM(self,fig1_ttl='PM',normP=True,show=True):
        normP11=np.ones_like(self.scatAng)
        figs=list()
        ttls=list()
        for i in np.arange(0,self.PMatVal.shape[0]):
            fig1,ax1=plt.subplots(3,2,sharex=True)
            if normP:
                normP11=np.squeeze(self.PMatVal[i,:,0])
                ax1[0,1].set_title('P22/P11')
                ax1[1,0].set_title('P33/P11')
                ax1[1,1].set_title('P44/P11')
                ax1[2,0].set_title('P12/P11')
                ax1[2,1].set_title('P34/P11')
            else:
                ax1[0,1].set_title('P22')
                ax1[1,0].set_title('P33')
                ax1[1,1].set_title('P44')
                ax1[2,0].set_title('P12')
                ax1[2,1].set_title('P34')
            ax1[0,0].plot(self.scatAng,np.squeeze(self.PMatVal[i,:,0]));ax1[0,0].set_title('P11')
            ax1[0,0].set_yscale('log')
            ax1[0,1].plot(self.scatAng,np.squeeze(self.PMatVal[i,:,1])/normP11);
            ax1[1,0].plot(self.scatAng,np.squeeze(self.PMatVal[i,:,2])/normP11);
            ax1[1,1].plot(self.scatAng,np.squeeze(self.PMatVal[i,:,3])/normP11);
            ax1[2,0].plot(self.scatAng,np.squeeze(self.PMatVal[i,:,4])/normP11);
            ax1[2,1].plot(self.scatAng,np.squeeze(self.PMatVal[i,:,5])/normP11);
 
          
            ax1[2,0].set_xlabel('Scattering angle')
            ax1[2,1].set_xlabel('Scattering angle')
            fig1.suptitle(fig1_ttl+'_%d'%i)
            if show:
                fig1.show()
            else:
                figs.append(fig1)
                ttls.append(fig1_ttl+'_%d'%i)
        if not(show):
            return figs,ttls
    def plot_Tau(self,v=None):
        '''
        Plot optical depth
        '''
        xcens=(self.xgrd[1:]+self.xgrd[0:-1])/2
        ycens=(self.ygrd[1:]+self.ygrd[0:-1])/2
        fig1,ax1=plt.subplots(subplot_kw={'aspect':'equal'})
        fig1_ttl=self.fname.split('.',1)[0]+'_COT'
        fig1.suptitle(fig1_ttl)
        tau=self.get_Tau()
        if v is None:
            ctf=ax1.contourf(xcens,ycens,tau,cmap=plt.cm.jet)
        else:
            ctf=ax1.contourf(xcens,ycens,tau,v,cmap=plt.cm.jet)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.2)
        if v is None:
            fig1.colorbar(ctf,cax=cax,label=r'$\tau$')
        else:
            fig1.colorbar(ctf,cax=cax,ticks=np.arange(v.min(),v.max()+1/10,10),label=r'$\tau$')
        ax1.set_xlabel('km')
        ax1.set_ylabel('km')
        fig1.show()
        return fig1,fig1_ttl
    def get_Tau(self,):
        Bex=self.extp3d#km^-1
        DZ=self.zgrd[1:]-self.zgrd[0:-1]#km
        return np.einsum('ijkl,j->kl',Bex,DZ)
    def writeLES_field(self,f_name,dpath=None):
        '''
        Keep f_name enter manually to avoid overite existing files.
        f_name='filename.nc'
        dpath='filepath/'
        '''
        if dpath==None:
            data = netCDF4.Dataset(f_name, 'w',format='NETCDF3_CLASSIC')
        else:
            data = netCDF4.Dataset(dpath+f_name, 'w',format='NETCDF3_CLASSIC')
        self.fname=f_name
        data.createDimension('nxp',self.nxp)
        data.createDimension('nyp',self.nyp)
        data.createDimension('nzp',self.nzp)
        data.createDimension('nx',self.nx)
        data.createDimension('ny',self.ny)
        data.createDimension('nz',self.nz)
        data.createDimension('nz3',self.nz3)
        data.createDimension('iz3l',self.iz3l)
        data.createDimension('npar1',self.npar1)
        data.createDimension('npar3',self.npar3)
        data.createDimension('nkd',self.nkd)
        data.createDimension('nstokes',self.nstokes)
        data.createDimension('nang',self.nang)
        data.createDimension('npf',self.npf)
        data.createDimension('npsfc',self.npsfc)
        
        data.createVariable('xgrd',np.float32,('nxp')) 
        data.createVariable('ygrd',np.float32,('nyp'))
        data.createVariable('zgrd',np.float32,('nzp'))
        data.createVariable('extp1d',np.float32,('npar1', 'nz'))
        data.createVariable('omgp1d',np.float32,('npar1', 'nz'))
        data.createVariable('jpfp1d',np.float32,('npar1', 'nz'))
        data.createVariable('rpfp1d',np.float32,('npar1', 'nz'))
        data.createVariable('extp3d',np.float32,('npar3', 'nz3', 'ny', 'nx'))
        data.createVariable('omgp3d',np.float32,('npar3', 'nz3', 'ny', 'nx'))
        data.createVariable('jpfp3d',np.float32,('npar3', 'nz3', 'ny', 'nx'))
        data.createVariable('rpfp3d',np.float32,('npar3', 'nz3', 'ny', 'nx'))
        data.createVariable('absg1d',np.float32,('nkd', 'nz'))
        data.createVariable('wkd',np.float32,('nkd'))
        data.createVariable('scatAng',np.float32,('nang')) 
        data.createVariable('PMatVal',np.float32,('npf', 'nang', 'nstokes'))
        data.createVariable('tsfc2d',np.float32,('ny', 'nx'))
        data.createVariable('jsfc2d',np.float32,('ny', 'nx'))
        data.createVariable('psfc2d',np.float32,('ny', 'nx', 'npsfc')) 
      
        data.variables['xgrd'][:]=self.xgrd
        data.variables['ygrd'][:]=self.ygrd
        data.variables['zgrd'][:]=self.zgrd
        data.variables['extp1d'][:]=self.extp1d
        data.variables['omgp1d'][:]=self.omgp1d
        data.variables['jpfp1d'][:]=self.jpfp1d
        data.variables['rpfp1d'][:]=self.rpfp1d
        data.variables['extp3d'][:]=self.extp3d
        data.variables['omgp3d'][:]=self.omgp3d
        data.variables['jpfp3d'][:]=self.jpfp3d
        data.variables['rpfp3d'][:]=self.rpfp3d
        data.variables['absg1d'][:]=self.absg1d
        data.variables['wkd'][:]=self.wkd
        data.variables['scatAng'][:]=self.scatAng 
        data.variables['PMatVal'][:]=self.PMatVal
        data.variables['tsfc2d'][:]=self.tsfc2d
        data.variables['jsfc2d'][:]=self.jsfc2d
        data.variables['psfc2d'][:]=self.psfc2d 
   
        data.close()
        print(f_name+' SAVED!')
        
def collectColsMSCART(fname00,fdpath,dsNew,out_name):
    #dsNew: new POLCARTdset object
    #fname00 'OP_dharma_008036_full_3_26_binsy0x0_MSCART_SSA020_SAA030_VAA000plus_NPH1e5.nc' (0,0 bin)
    #fdpath= 'results/OP_dharma_008036_full_3_26_bins/' (datapath)
    dsNew.MeanPRad=np.zeros((122,128,128,4))
    dsNew.RMSEPRad=np.zeros((122,128,128,4))
    dsNew.fname=fname00[0:32]+'x'+fname00[35:]
    for i in np.arange(0,128,1):
        for j in np.arange(0,128,1):
            fname=fname00[0:32]+str(i)+'x'+str(j)+fname00[35:]
            ds1D=POLCARTdset('ds1D','OP_dharma_008036_full_3_26_bins/')
            print(fname)
            ds1D.readMSCARTplus(fname,fdpath=fdpath,clm=True,prnt=False)
            if i==0 and j==0:
                dsNew.SZA=ds1D.SZA
                dsNew.SAA=ds1D.SAA
                dsNew.NPH=ds1D.NPH
                dsNew.Nbat=ds1D.Nbat
                dsNew.MeanTiming=ds1D.MeanTiming
                dsNew.VZA=ds1D.VZA
                dsNew.VAA=ds1D.VAA
                dsNew.cc3D=ds1D.cc3D
            dsNew.MeanPRad[:,ds1D.pnt[0],ds1D.pnt[1],:]=ds1D.MeanPRad.T
            dsNew.RMSEPRad[:,ds1D.pnt[0],ds1D.pnt[1],:]=ds1D.RMSEPRad.T
    dsNew.savePOLCARThdf5(out_name)

def scat_ang(SZA,VZA,RAA):
    '''
    Return the scattering angle in degrees
    SZA,VZA : Source(ray)_Zenith_Angle, Viewing_Zenith_Angle (in degrees)
    RAA : Source Azimuth Angle - (Apparent)Scattering Azimuth
    '''
#    RAA=0-180-RAA
# scat_a=cos(VZA)*cos(SZA)+sin(VZA)*sin(SZA)*cos(SAA-VAA)
    SZA=np.deg2rad(SZA);
    VZA=np.deg2rad(VZA)
    RAA=np.deg2rad(RAA)
    
    scat_a=np.arccos(np.cos(SZA)*np.cos(VZA)+np.sin(SZA)*np.sin(VZA)*np.cos(RAA))
    scat_a=np.rad2deg(scat_a)
    return scat_a
    
def getKokhanovsky_PM(ac='aero'):
    '''
    ac:'aoro' or 'cld'
    Read Kokhanovsky Phase matrix from a text file and return the Pij
    '''
    if ac=='aero':
        fl=open('/umbc/xfs1/zzbatmos/users/charaj1/PDA_python_wrapper/Kokhanovsky_aero/Kokhanovsky_benchmark_aerosol.txt','r')
    elif ac=='cld':
        fl=open('/umbc/xfs1/zzbatmos/users/charaj1/PDA_python_wrapper/Kokhanovsky_cld/Kokhanovsky_benchmark_cloud.txt','r')
    else:
        print('Give \'aero\' or \'cld\' as ac')
    data=fl.readlines()
    header=data[0:6]
    values=np.loadtxt(data[6:])
    fl.close()
    ang=values[:,0]
    P11=values[:,1]
    P22=values[:,2]*P11
    P33=values[:,3]*P11
    P44=values[:,4]*P11
    P12=values[:,5]*P11
    P34=values[:,6]*P11

    return ang,P11,P22,P33,P44,P12,P34

class LES_case(object):
    '''
    3D LES RT simulations with the input MSCART field
    --------------------------------------------------
    - Old __init__ function was changed. If you are getting errors when using with many input parameters, that 
    could be the reason. Old function is self.read_rt_and_field_old() now.
    - New __init__ only takes one input parameter which defines the complete case.
    name: filename including the path -or- case name
        case name ex. DYCOMS2_120_b0p860, filenameExp=False
        filename +path ex. name=['path_and_filename_3d','path_and_filename_1d'], filenameExp=True
    '''
    def __init__(self,name,filenameExp=False):
        self.rotate1D='not_done'
        self.cloud_mask='not_done'
        
        if filenameExp:
            fname_3d=name[0]
            fname_1D=name[1]
        else:
            fname_3d=self.get_file_name(name+'_3D')
            fname_1D=self.get_file_name(name+'_1D')
        self.band=fname_3d.split('results/',1)[1].split('/',1)[0]
        #3D RT
        result_dir3d='results/'+self.band+'/'
        base_dir3d=fname_3d.split('results',1)[0]
        self.RT=POLCARTdset('LES',base_dir3d)
        MSCARThdf=fname_3d.split(self.band+'/',1)[1]
        self.RT.readPOLCARThdf5(MSCARThdf,dpath=base_dir3d+result_dir3d)
        self.RT_field=LES_field(MSCARThdf.split('_MSCART',1)[0]+'.nc',dpath=base_dir3d)
        self.RT_field.readLES_field()
        #1D RT
        
        base_dir1d=fname_1D.split('1Druns',1)[0]
        self.RT1D=POLCARTdset('LES1D',base_dir1d+'1Druns/LES'+self.band+'_bins/')
        MSCARThdf1D=fname_1D.split('LES'+self.band+'_bins/',1)[1]
        self.RT1D.readPOLCARThdf5(MSCARThdf1D,base_dir1d+'1Druns/results/LES'+self.band+'_bins/')
        self.rotate_1D_domain()
        self.xcens=(self.RT_field.xgrd[1:]+self.RT_field.xgrd[0:-1])/2
        self.ycens=(self.RT_field.ygrd[1:]+self.RT_field.ygrd[0:-1])/2  
        if type(self.RT.fname) is np.ndarray:
            self.RT.fname=self.RT.fname[0].astype(str)
        if type(self.RT1D.fname) is np.ndarray:
            self.RT1D.fname=self.RT1D.fname[0].astype(str)
        self.RT.remove_redundant_nadir()
        self.RT1D.remove_redundant_nadir()
    def read_rt_and_field_old(self,MSCARThdf,base_dir,result_dir=None,RT1Dname=None,band=None):
        '''
        Old __init__() function of this class. Atleast hdf file name and based_dir should be given.
        -----------------------------------------------------
        MSCARThdf: *.hdf5 name (with extension)
        base_dir: Directory of the 3D runs
        result_dir(None): Give 3D result dir from the base_dir o/w tries to create one based on conventions.
        w1D(False): True if want to read corresponding 1D simulations
        band(None): '0p860' string o/w tries to create one based on conventions.

        '''
        self.rotate1D='not_done'
        if band==None:
            self.band=MSCARThdf.split('_b',1)[1].split('_',1)[0]
        else:
            self.band=band
        if result_dir==None:
            result_dir='results/b'+self.band+'/'
        else:
            result_dir=result_dir
        self.RT=POLCARTdset('LES',base_dir)
        self.RT.readPOLCARThdf5(MSCARThdf,dpath=base_dir+result_dir)
        self.RT_field=LES_field(MSCARThdf.split('_MSCART',1)[0]+'.nc',dpath=base_dir)
        self.RT_field.readLES_field()
        if not(RT1Dname==None):
            self.RT1D=POLCARTdset('LES1D',base_dir+'1Druns/LESb'+self.band+'_bins/')
            MSCARThdf1D=RT1Dname#MSCARThdf.split('_SZA')[0]+'_1D_bins_'+MSCARThdf.split('MSCART_',1)[1]
            self.RT1D.readPOLCARThdf5(MSCARThdf1D,base_dir+'1Druns/results/LESb'+self.band+'_bins/')
            self.rotate_1D_domain()
        self.xcens=(self.RT_field.xgrd[1:]+self.RT_field.xgrd[0:-1])/2
        self.ycens=(self.RT_field.ygrd[1:]+self.RT_field.ygrd[0:-1])/2  
        if type(self.RT.fname) is np.ndarray:
            self.RT.fname=self.RT.fname[0].astype(str)
        if type(self.RT1D.fname) is np.ndarray:
            self.RT1D.fname=self.RT1D.fname[0].astype(str)
        self.RT.remove_redundant_nadir()
        self.RT1D.remove_redundant_nadir()
    def get_file_name(self,case):
        '''
        Return the appropriate file name
        When a new run was done, change the filename here. Hopefully, the others will manage it.
        '''
        filename={}
        filename['DYCOMS2_120_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/results/b0p860/DYCOMS2_dharma_008036_b0p860_MSCART_SZA120_SAA000_VAA000plus_NPH2e6.hdf5'
        filename['DYCOMS2_140_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/results/b0p860/DYCOMS2_dharma_008036_b0p860_MSCART_SZA140_SAA000_VAA000plus_NPH2e6.hdf5'
        filename['DYCOMS2_160_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/results/b0p860/DYCOMS2_dharma_008036_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH2e6.hdf5'
        filename['DYCOMS2_120_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/results/b2p13/DYCOMS2_dharma_008036_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH2e6.hdf5'
        filename['DYCOMS2_140_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/results/b2p13/DYCOMS2_dharma_008036_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH2e6.hdf5'
        filename['DYCOMS2_160_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/results/b2p13/DYCOMS2_dharma_008036_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH2e6.hdf5'
        filename['DYCOMS2_120_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/1Druns/results/LESb0p860_bins/DYCOMS2_dharma_008036_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['DYCOMS2_140_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/1Druns/results/LESb0p860_bins/DYCOMS2_dharma_008036_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['DYCOMS2_160_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/1Druns/results/LESb0p860_bins/DYCOMS2_dharma_008036_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['DYCOMS2_120_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/1Druns/results/LESb2p13_bins/DYCOMS2_dharma_008036_b2p13_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['DYCOMS2_140_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/1Druns/results/LESb2p13_bins/DYCOMS2_dharma_008036_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['DYCOMS2_160_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/1Druns/results/LESb2p13_bins/DYCOMS2_dharma_008036_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
    
        filename['RICO_120_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b0p860/RICO_dharma_005044_b0p860_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['RICO_140_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b0p860/RICO_dharma_005044_b0p860_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['RICO_160_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b0p860/RICO_dharma_005044_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5'
        #filename['RICO_120_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b2p13/RICO_dharma_005044_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['RICO_120_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b2p13/RICO_dharma_005044_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH2e6.hdf5'
        #filename['RICO_140_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b2p13/RICO_dharma_005044_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['RICO_140_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b2p13/RICO_dharma_005044_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH2e6.hdf5'
        #filename['RICO_160_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b2p13/RICO_dharma_005044_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['RICO_160_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/results/b2p13/RICO_dharma_005044_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH2e6.hdf5'
        filename['RICO_120_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/1Druns/results/LESb0p860_bins/RICO_dharma_005044_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['RICO_140_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/1Druns/results/LESb0p860_bins/RICO_dharma_005044_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['RICO_160_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/1Druns/results/LESb0p860_bins/RICO_dharma_005044_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['RICO_120_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/1Druns/results/LESb2p13_bins/RICO_dharma_005044_b2p13_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['RICO_140_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/1Druns/results/LESb2p13_bins/RICO_dharma_005044_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['RICO_160_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/1Druns/results/LESb2p13_bins/RICO_dharma_005044_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
    
        filename['ATEXc_120_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/results/b0p860/ATEXc_dharma_007877_b0p860_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXc_140_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/results/b0p860/ATEXc_dharma_007877_b0p860_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXc_160_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/results/b0p860/ATEXc_dharma_007877_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXc_120_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/results/b2p13/ATEXc_dharma_007877_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXc_140_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/results/b2p13/ATEXc_dharma_007877_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXc_160_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/results/b2p13/ATEXc_dharma_007877_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5'
        
        filename['ATEXc_120_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/1Druns/results/LESb0p860_bins/ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXc_140_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXc/1Druns/results/LESb0p860_bins/ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        #filename['ATEXc_140_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/1Druns/results/LESb0p860_bins/ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXc_160_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXc/1Druns/results/LESb0p860_bins/ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
        #filename['ATEXc_160_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/1Druns/results/LESb0p860_bins/ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXc_120_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/1Druns/results/LESb2p13_bins/ATEXc_dharma_007877_b2p13_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXc_140_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXc/1Druns/results/LESb2p13_bins/ATEXc_dharma_007877_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        #filename['ATEXc_140_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/1Druns/results/LESb2p13_bins/ATEXc_dharma_007877_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXc_160_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXc/1Druns/results/LESb2p13_bins/ATEXc_dharma_007877_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e6.hdf5'
		#filename['ATEXc_160_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXc/1Druns/results/LESb2p13_bins/ATEXc_dharma_007877_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
    
        filename['ATEXp_120_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/results/b0p860/ATEXp_dharma_013067_b0p860_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXp_140_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/results/b0p860/ATEXp_dharma_013067_b0p860_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXp_160_b0p860_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/results/b0p860/ATEXp_dharma_013067_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXp_120_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/results/b2p13/ATEXp_dharma_013067_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXp_140_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/results/b2p13/ATEXp_dharma_013067_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXp_160_b2p13_3D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/results/b2p13/ATEXp_dharma_013067_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5'
        filename['ATEXp_120_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXp/1Druns/results/LESb0p860_bins/ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        #filename['ATEXp_120_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/1Druns/results/LESb0p860_bins/ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXp_140_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXp/1Druns/results/LESb0p860_bins/ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        #filename['ATEXp_140_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/1Druns/results/LESb0p860_bins/ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXp_160_b0p860_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/1Druns/results/LESb0p860_bins/ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXp_120_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/1Druns/results/LESb2p13_bins/ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXp_140_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/1Druns/results/LESb2p13_bins/ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5'
        filename['ATEXp_160_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXp/1Druns/results/LESb2p13_bins/ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e6.hdf5'
        #filename['ATEXp_160_b2p13_1D']='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/1Druns/results/LESb2p13_bins/ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5'
        
        return filename[case]
    def rotate_1D_domain(self,):
        '''
        Transpose MeanPRad and RMSEPRad arrays of the 1D results to be matched with 3D
        '''
        if self.rotate1D=='not_done':
            self.rotate1D='done'
            self.RT1D.MeanPRad=self.RT1D.MeanPRad.swapaxes(1,2)
            self.RT1D.RMSEPRad=self.RT1D.RMSEPRad.swapaxes(1,2)
            #print('1D domain rotated')
        elif self.rotate1D=='done':
            print('Already rotated')
    def find_3DimpFac(self,VZAi):
        '''
        Find 3D impact factor
            f=(3D-1D)/1D
        VZAi: viewing zenith angle index (usually 61 for zero degrees)
        '''
        self.bias3D1D =self.RT.MeanPRad[VZAi,:]-self.RT1D.MeanPRad[VZAi,:]
        self.impFac3D =self.bias3D1D/self.RT1D.MeanPRad[VZAi,:]
    def get_bin_defs(vci):
        '''
        To get pre-defined bin arrays for MSCART LES reflectance simulations
        vci: Give index (ix) as the following table
            Contourf color bars and v (VC)
            +++++++++++++++++++++++++++++++++++++++
            DYCOMS2
            +++++++++++++++++++++++++++++++++++++++
            ---------------------------------------
            band   |  SZA  |  VZA  |  ix  |  SAA  |
            ---------------------------------------
            0p860  |  140  |  000  |   1  |  030  |
            2p13   |  140  |  000  |   2  |  030  |
            0p860  |  120  |  000  |   3  |  030  |
            2p13   |  120  |  000  |   4  |  030  |
            ---------------------------------------
            2p13   |  140  |  000  |   5  |  000  |
            0p860  |  140  |  000  |   6  |  000  |
            +++++++++++++++++++++++++++++++++++++++
            ATEXc
            +++++++++++++++++++++++++++++++++++++++
            ---------------------------------------
            0p860  |  120  |  000  |   7  |  000  |
            2p13   |  120  |  000  |   8  |  000  |
        RTdim: '1D' or '3D' RT transfer string
        '''
        VC={1:{'vIR':np.linspace(0   ,1,50) ,'vQR':np.linspace(-.05,0,50)    ,'vUR':np.linspace(0.0,0.1,50) ,\
               'vIe':np.linspace(0.0,2.0,20),'vQe':np.linspace(0.0,1.00,20)  ,'vUe':np.linspace(0.0,1.0,20),\
               'cIR':np.arange(0,1.1,0.25)  ,'cQR':np.arange(-.05,0.011,0.02),'cUR':np.arange(0.0,0.11,0.05),\
               'cIe':np.arange(0.0,2.1,0.5) ,'cQe':np.arange(0.0,1.1,0.5)    ,'cUe':np.arange(0.0,1.1,0.5)},
            2:{'vIR':np.linspace(0.1 ,0.3,50) ,'vQR':np.linspace(-.02,0,50)    ,'vUR':np.linspace(0.02,0.03,50) ,\
               'vIe':np.linspace(0.75,2.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,2,20),\
               'cIR':np.arange(0.1,0.31,0.1)  ,'cQR':np.arange(-.02,0.001,0.01),'cUR':np.arange(0.02,0.031,0.005),\
               'cIe':np.arange(0.75,2.1,0.25) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,2.1,0.5)},\
            3:{'vIR':np.linspace(0.2 ,0.8,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(-0.05,0.05,50) ,\
               'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
               'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(-0.05,0.051,0.05),\
               'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
            4:{'vIR':np.linspace(0.1 ,0.3,50) ,'vQR':np.linspace(-.005,0,50)   ,'vUR':np.linspace(-0.05,0.05,50) ,\
               'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,2,20),\
               'cIR':np.arange(0.1,0.31,0.1)  ,'cQR':np.arange(-.005,.001,.002),'cUR':np.arange(-0.05,0.051,0.05),\
               'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,2.1,0.5)},\
            5:{'vIR':np.linspace(0.2 ,0.4,50) ,'vQR':np.linspace(-.05,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
               'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)  ,'vUe':np.linspace(0.0,2,20),\
               'cIR':np.arange(0.2,0.41,0.1)  ,'cQR':np.arange(-.05,.01,.02),'cUR':np.arange(-0.05,0.051,0.05),\
               'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)   ,'cUe':np.arange(0,2.1,0.5)},\
            6:{'vIR':np.linspace(0.2 ,0.82,50) ,'vQR':np.linspace(-.1,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
               'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)  ,'vUe':np.linspace(0.0,2,20),\
               'cIR':np.arange(0.2,0.82,0.2)  ,'cQR':np.arange(-.1,.01,.05),'cUR':np.arange(-0.05,0.051,0.05),\
               'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)   ,'cUe':np.arange(0,2.1,0.5)},\
            7:{'vIR':np.linspace(0.2 ,0.8,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
               'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
               'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.010),\
               'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
            8:{'vIR':np.linspace(0.0 ,0.4,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
               'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
               'cIR':np.arange(0.0,0.41,0.1)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.010),\
               'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)}    
            }        
        return VC[vci]
    def add_cloud_mask(self,cldM):
        '''
        Mask 'cloud free' pixels
        cldM: 2D array. 'True' are cloudy pixels
        '''
        if self.cloud_mask=='not_done':
            self.RT.MeanPRad[:,np.invert(cldM),:]=np.nan
            self.RT.RMSEPRad[:,np.invert(cldM),:]=np.nan
            self.RT1D.MeanPRad[:,np.invert(cldM),:]=np.nan
            self.RT1D.RMSEPRad[:,np.invert(cldM),:]=np.nan
            self.cloud_mask='done'
        else:
            print('Already applied a cloud mask!')

class DHARMA(object):
    '''
    To read DHARMA outputs
    fname: filename (with extension)
    fpath(None): path to the file
    '''
    def __init__(self,fname,fpath=None):
        self.fname=fname
        if not(fpath==None):
            self.fpath=fpath
            self.readDHARMA(self.fname,self.fpath)
        else:
            self.readDHARMA(self.fname)
    def readDHARMA(self,fname,fdpath=None):
        if not(fdpath==None):
            filename=fdpath+fname
        else:
            filename=fname
        f = netCDF4.Dataset(filename,'r')
        self.x=f.variables['x'][:]            #(nx)E location of grid box center(m) 
        self.y=f.variables['y'][:]            #(ny)N location of grid box center(m) 
        self.z=f.variables['z'][:]            #(nz)altitude of grid box center(m)
        self.zw=f.variables['zw'][:]          #(nzp1)altitude of grid box boundary(m)
        self.r_drops=f.variables['r_drops'][:]#(nbin)radius of mass-weighted droplet bin center(um)
        self.rbound_drops=f.variables['rbound_drops'][:]#(nbinp1)radius at boundary between droplet bins(um)
        self.time=f.variables['time'][:]      #elapsed simulation time(s)
        try:
            self.rho_air=f.variables['rho_air'][:]#(nz)air density(kg/m^3)
        except KeyError:
            print('rho_air does not exist in '+filename)
        self.P=f.variables['P'][:]            #(nz)air pressure(mb)
        self.T=(f.variables['T'][:]).swapaxes(1,2)#(nz,nx,ny)air temperature(K) 
        self.qv=(f.variables['qv'][:]).swapaxes(1,2)#(nz,nx,ny)water vapor mixing ratio(g/kg) 
        self.dN_drops=(f.variables['dN_drops'][:]).swapaxes(2,3)#(nbin,nz,ny,nx)number concentration of activated droplets in size bin(cm^-3) 
        try:
            self.w=(f.variables['w'][:]).swapaxes(1,2)#(nzp1,ny,nx)vertical wind(m/s)
        except KeyError:
            print('w does not exist in '+filename)
        f.close()
    def selectDomain(self,x_bound,y_bound,z_bound):
        '''
        to select a particular domain and save as a new file.
        For the selected Towering Cumulus case for Dan and collaborators.
        '''
        x=self.x[x_bound[0],x_bound[1]]
        y=self.y[y_bound[0],y_bound[1]]


class DHARMA_onmp(object):
    '''
    Optical and microphysical properties of the DHARMA output
    '''
    def __init__(self,fname,fpath=None):
        self.cloud_mask='not_done'
        self.fname=fname
        self.DHARMA=DHARMA(self.fname,fpath)
        
    def readDHARMA(self,):
        self.x=self.DHARMA.x
        self.y=self.DHARMA.y
        self.dz=self.DHARMA.zw[1:]-self.DHARMA.zw[0:-1]#m
        self.lwc,self.lwp=self.findLWP(self.DHARMA)
        self.BextQe2p09=self.find_Bext(self.DHARMA)#m^-1
        self.tauQe2p09=np.einsum('i,ijk->jk',self.dz,self.BextQe2p09)
    def readMie(self,mie_name,mie_path=None):
        if mie_path==None:
            c_mie=MieSet(mie_name)#cloud_mie object
        else:
            c_mie=MieSet(mie_name,mie_path)
        c_mie.readMie()
        self.c_mie=c_mie       
        print('self.c_mie created!')
    def setupMieOut_crude(self,mie_name):
        '''
        To read Mie output and setup required variables to generate MSCART field_file.
        Have to run Mie code first and know the name and the location of the Mie output *.nc
        
        Consider only the bin center radii for Mie calculations. Since the bin withs are considerably large,
        the best approach would be consider finer radii separations for Mie calculations and average over the
        each bin width in order to get the required optical proerties (Qe,alb,Pij) of each bin. This approach is
        implemeted in setupMieOut().
        -------------------------------------------------------------
        mie_name: string name of the Mie output file without *.nc
        return: extp3d,omgp3d,P11,P33,P12,P34,cloud_mie
        '''
        case=self
        self.readMie(mie_name)
        extp3d=np.einsum('ij,i,iklm->jiklm',self.c_mie.qe,case.DHARMA.r_drops**2,case.DHARMA.dN_drops)*np.pi*1e-6#m^-1
        s1=case.DHARMA.dN_drops.shape    
        s2=self.c_mie.wvl.size
        omgp3d=np.ones((s2,s1[0],s1[1],s1[2],s1[3]),dtype=float)
        omgp3d=np.einsum('ij,jiklm->jiklm',self.c_mie.alb,omgp3d)
        P11=self.c_mie.P11.swapaxes(0,1)
        P33=self.c_mie.P33.swapaxes(0,1)
        P12=self.c_mie.P12.swapaxes(0,1)
        P34=self.c_mie.P34.swapaxes(0,1)    
        return extp3d,omgp3d,P11,P33,P12,P34,self.c_mie
        
    def get_pDefr_l(self,lesCname):
        '''
        Predefined radii bin edges for each LES case
        ----------------------------------------------
        An example of finding bin edges
            MieName='DYCOM2_dharma_008036_mie_470_860_2p13'
            cloud_mie=Mie.MieSet(MieName)
            cloud_mie.readMie()
            r=cloud_mie.d/2
            lesBinEd=np.zeros_like(ATEXc.DHARMA.rbound_drops,dtype=int)
            for i in np.arange(0,ATEXc.DHARMA.rbound_drops.size):
                lesBinEd[i]=(abs(r-ATEXc.DHARMA.rbound_drops[i])).argmin()
        lesCname: DYCOMS2,TowerSc,TowerSc2,ATEXc,ATEXp,RICO
        '''
        pDefr_l={'DYCOMS2':np.array([0,8,18,31,47,68,93,126,166,218,282,364,467,596,759,964,\
              1223,1549,1959,2477,3128,3949,4984,6288,7930,9999],dtype=int),\
                 'TowerSc':np.array([0,    0,    1,    2,    2,    3,    4,    5,    7,    8,   10,
         12,   15,   18,   22,   27,   32,   38,   45,   54,   64,   77,
         91,  108,  128,  152,  180,  213,  252,  298,  353,  417,  493,
        583,  690,  815,  964, 1139, 1347, 1592, 1882, 2224, 2628, 3106,
       3671, 4339, 5127, 6059, 7160, 8461, 9999],dtype=int),\
                'TowerSc2':np.array([   0,    1,    3,    5,    7,   10,   13,   16,   21,   26,   32,
         39,   47,   57,   69,   83,   99,  119,  142,  169,  201,  239,
        284,  337,  399,  422,  449,  481,  519,  564,  617,  680,  754,
        842,  945, 1068, 1212, 1383, 1585, 1823, 2105, 2438, 2831, 3296,
       3845, 4494, 5261, 6167, 7238, 8504, 9999]),\
                'ATEXc':np.array([   0,    8,   18,   31,   47,   68,   93,  126,  166,  218,  282,
        364,  466,  596,  759,  964, 1223, 1549, 1959, 2476, 3128, 3949,
       4984, 6287, 7930, 9999]),\
                'ATEXp':np.array([   0,    8,   18,   31,   47,   68,   93,  126,  166,  218,  282,
        364,  466,  596,  759,  964, 1223, 1549, 1959, 2476, 3128, 3949, 4984, 6287, 7930, 9999]),
                'RICO':np.array([   0,    1,    2,    3,    5,    7,    9,   12,   16,
                     21,   28,   36,   46,   59,   75,   95,  121,  153,
                    194,  245,  309,  391,  493,  622,  784,  989, 1247,
                   1572, 1982, 2497, 3147, 3966, 4998, 6298, 7936, 9999])
                }#r bin edge indices
        return pDefr_l[lesCname]
        
    def setupMieOut(self,mie_name,lesBinEd=None,lesCname=None,mie_path=None):
        '''
        To read Mie output and setup required variables to generate MSCART field_file.
        Have to run Mie code first and know the name and the location of the Mie output *.nc
        Have to provide either lesBinEd or lesCname.
        
        Optical properties will be averaged over each bin width separation (given as lesBinEd array or specified as lesCname) 
        and output the appropriate "bin optical proerties" to setup MSCART field file
        ------------------------------------------------------------------------------------
        mie_name: String name of the Mie output file without *.nc
        lesBinEd: Array of the LES radii bin edges
            Indices of r that are nearest to the bind edge values in DHARMA.rbound_drops
        lesCname: String name to get predefined radii bin edges for specific case.
            lesCname='DYCOMS2' for DYCOMS2 dharma_008036 bin edges
            lesCname='TowerSc' for Towering Sc case (Dan and Ann case)
        return:
            extp3d,omgp3d,aP11,aP33,aP12,aP34,c_mie
        '''
        self.readMie(mie_name,mie_path)
        qe=self.c_mie.qe
        al=self.c_mie.alb
        r=self.c_mie.d/2
        P11,P33,P12,P34=self.c_mie.P11,self.c_mie.P33,self.c_mie.P12,self.c_mie.P34
        if lesCname==None and lesBinEd==None:
            print('Assign either lesBinEd or lesCname')
        elif not(lesBinEd==None):
            r_l=lesBinEd
        elif not(lesCname==None):
            print('Make sure to give correct lesCname.')
            try:
                r_l=self.get_pDefr_l(lesCname)
            except KeyError:
                print('KeyError: Do not have predefined radii bin edges for this LES case yet!!')
                print('An example to get binedges is shown in self.get_pDefr_l function help.')
        #Bext= sum(Qe(r_i)*pi*r_i^2*dN_i)
        qe_avg=np.zeros((qe.shape[1],r_l.size-1),dtype=float)
        alb_av=np.zeros_like(qe_avg,dtype=float)
        aP11=np.zeros((qe.shape[1],r_l.size-1,self.c_mie.ang.size),dtype=float)
        aP33=np.zeros_like(aP11,dtype=float)
        aP12=np.zeros_like(aP11,dtype=float)
        aP34=np.zeros_like(aP11,dtype=float)
        for i in np.arange(0,r_l.size-1,1):
            lt=r_l[i]
            rt=r_l[i+1]
            qe_avg[:,i]=np.einsum('ij,i->j'        ,qe[lt:rt,:],r[lt:rt]**2)/np.sum(r[lt:rt]**2)
            alb_av[:,i]=np.einsum('ij,i,ij->j'     ,qe[lt:rt,:],r[lt:rt]**2,al[lt:rt,:],)/np.einsum('i,ij->j',r[lt:rt]**2,qe[lt:rt])
            pn11       =np.einsum('ij,i,ij,ijk->jk',qe[lt:rt,:],r[lt:rt]**2,al[lt:rt,:],P11[lt:rt,:])
            pn33       =np.einsum('ij,i,ij,ijk->jk',qe[lt:rt,:],r[lt:rt]**2,al[lt:rt,:],P33[lt:rt,:])
            pn12       =np.einsum('ij,i,ij,ijk->jk',qe[lt:rt,:],r[lt:rt]**2,al[lt:rt,:],P12[lt:rt,:])
            pn34       =np.einsum('ij,i,ij,ijk->jk',qe[lt:rt,:],r[lt:rt]**2,al[lt:rt,:],P34[lt:rt,:])
            pd         =np.einsum('ij,i,ij->j'     ,qe[lt:rt,:],r[lt:rt]**2,al[lt:rt,:])
            aP11[:,i,:]=np.einsum('ij,i->ij'      ,pn11,1/pd)
            aP33[:,i,:]=np.einsum('ij,i->ij'      ,pn33,1/pd)
            aP12[:,i,:]=np.einsum('ij,i->ij'      ,pn12,1/pd)
            aP34[:,i,:]=np.einsum('ij,i->ij'      ,pn34,1/pd)
            
        extp3d=np.einsum('ji,i,iklm->jiklm',qe_avg,self.DHARMA.r_drops**2,self.DHARMA.dN_drops)*np.pi*1e-3#km^-1
        s1=self.DHARMA.dN_drops.shape    
        s2=self.c_mie.wvl.size
        omgp3d=np.ones((s2,s1[0],s1[1],s1[2],s1[3]),dtype=float)
        omgp3d=np.einsum('ji,jiklm->jiklm',alb_av,omgp3d)
        return extp3d,omgp3d,aP11,aP33,aP12,aP34,self.c_mie
    def readNsetFromRSPMieOut(self,):
        '''
        To read Dan's RSP 50 bin Mie calculations for ToweringSc case.
        '''
        filename='/umbc/xfs1/zzbatmos/users/charaj1/microphysics/TSc/sngl_bin_blk_tCU_DHARMA_25bin_RSP_B'
        aP11,aP33,aP12,aP34=[np.zeros((50,1801),dtype=float) for _ in range(4)]
        qe_avg,alb_av=[np.zeros(50,dtype=float) for _ in range(2)]
        for i in np.arange(0,50):
            f=h5py.File(filename+'%02d.h5'%(i+1),'r')
            aP11[i,:]=f['P11_blk_avg'][:]
            aP33[i,:]=f['P33_blk_avg'][:]
            aP12[i,:]=f['P12_blk_avg'][:]
            aP34[i,:]=f['P34_blk_avg'][:]
            qe_avg[i]=f['Qe_blk_avg'][0]
            alb_av[i]=f['al_blk_avg'][0]
        ang=f['ang'][:]
        wvl=np.asarray(0.865)
        extp3d=np.einsum('i,i,iklm->iklm',qe_avg,self.DHARMA.r_drops**2,self.DHARMA.dN_drops)*np.pi*1e-3#km^-1
        s1=self.DHARMA.dN_drops.shape    
        omgp3d=np.ones((s1[0],s1[1],s1[2],s1[3]),dtype=float)
        omgp3d=np.einsum('i,iklm->iklm',alb_av,omgp3d)
        
        return extp3d,omgp3d,aP11,aP33,aP12,aP34,wvl,ang
    def setup_dTau(self,mie_name,lesCname,mie_path=None,lesBinEd=None):
        '''
        Use Mie results to get COT as the vertical axis.
        Will setup self.zTau[band,z,x,y]
        mie_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/'
        mie_name='DYCOM2_dharma_008036_mie_470_860_2p13'
        see self.get_pDefr_l for lesCname/lesBinEd
        '''
        qe_avg,_=self.get_qe_avg(mie_name,lesCname,mie_path,lesBinEd)
        extp3d=np.einsum('ji,i,iklm->jiklm',qe_avg,self.DHARMA.r_drops**2,self.DHARMA.dN_drops)*np.pi*1e-3#km^-1
        self.dTau=np.einsum('ijzxy,z->izxy',extp3d,self.dz*1e-3)   
    def setup_zTau0top(self,):
        '''
        Setting-up tau vertical axis. Both 0-top and o-bottom
        self.dTau must be computed before callling this.
        '''
        self.zTau=np.cumsum(self.dTau,axis=1)
        self.zTau0top=np.zeros_like(self.zTau,dtype=float)
        for i in np.arange(self.x.size):
            for j in np.arange(self.y.size):
                for wv in np.arange(3):
                    self.zTau0top[wv,:,i,j]=self.zTau[wv,:,i,j].max()-self.zTau[wv,:,i,j]
        print('self.zTau, self.zTau0top created.')
    def get_qe_avg(self,mie_name,lesCname,mie_path=None,lesBinEd=None):
        self.readMie(mie_name,mie_path)
        c_mie=self.c_mie
        qe=c_mie.qe
        al=c_mie.alb
        r=c_mie.d/2
        if lesCname==None and lesBinEd==None:
            print('Assign either lesBinEd or lesCname')
        elif not(lesBinEd==None):
            r_l=lesBinEd
        elif not(lesCname==None):
            print('Make sure to give correct lesCname.')
            try:
                r_l=self.get_pDefr_l(lesCname)
            except KeyError:
                print('KeyError: Do not have predefined radii bin edges for this LES case yet!!')
                print('An example to get binedges is shown in self.get_pDefr_l function help.')
        #Bext= sum(Qe(r_i)*pi*r_i^2*dN_i)
        qe_avg=np.zeros((qe.shape[1],r_l.size-1),dtype=float)
        alb_av=np.zeros_like(qe_avg,dtype=float)
        for i in np.arange(0,r_l.size-1,1):
            lt=r_l[i]
            rt=r_l[i+1]
            qe_avg[:,i]=np.einsum('ij,i->j'        ,qe[lt:rt,:],r[lt:rt]**2)/np.sum(r[lt:rt]**2)  
            alb_av[:,i]=np.einsum('ij,i,ij->j'     ,qe[lt:rt,:],r[lt:rt]**2,al[lt:rt,:],)/np.einsum('i,ij->j',r[lt:rt]**2,qe[lt:rt])
            
        self.qe_avg=qe_avg; self.alb_av=alb_av
        print('self.qe_avg, self.alb_av created')
        return qe_avg,c_mie
    def find_reVW(self,mie_name,lesCname,dgSZA,dgVZA,a=1,b=0,lesBinEd=None,mie_path=None,band=None):
        '''
        Two parameter bidirectional weighting function implementation for both re and ve
        Platnick et. al. (2000), Alexandrov et. al. (2012), Zhang et. al. (2017), Miller et. al. (2017)
        -----------------------------------------------------------------------------------------------
        mie_name: Mie file name without extenction
        lesCname: see self.get_pDefr_l for lesCname/lesBinEd
        dgSZA,dgVZA: Solar zenith and viewing zenith in degrees respectively
        a,b: W=tau**b*exp(-a*B*tau)
        '''
        start=time.time()
        if lesCname is not None:
            self.setup_dTau(mie_name,lesCname,mie_path=mie_path)
        elif lesBinEd is not None:
            self.setup_dTau(mie_name,lesBinEd=lesBinEd,mie_path=mie_path)
        if band is None:
            band=1
            print('Deafault wavelength %0.3f was selected'%(self.c_mie.wvl[band]))
            print('Available wavelengths:'+str(self.c_mie.wvl))
        self.setup_zTau0top()
        tau=self.zTau0top[band,:]
        sza=np.deg2rad(dgSZA)
        vza=np.deg2rad(dgVZA)
        B=1/np.cos(sza)+1/np.cos(vza)
        c=1/np.einsum('zxy,zxy->xy',tau**b*np.exp(-a*B*tau),self.dTau[band,:])
        #c=1/np.trapz(tau**b*np.exp(-a*B*tau),tau,axis=0)
        w_tauDIVc=tau**b*np.exp(-a*B*tau)
        print('Calculating re and ve ...')
        #calculating re_tau (Hansen & Travis 1974 2.53)
        num=np.einsum('r,rzxy->zxy',self.qe_avg[band,:]*self.alb_av[band,:]*self.DHARMA.r_drops**3,self.DHARMA.dN_drops)
        #num=np.trapz(num,self.DHARMA.r_drops,axis=0) dN_drops is dN, not n(r)
        den=np.einsum('r,rzxy->zxy',self.qe_avg[band,:]*self.alb_av[band,:]*self.DHARMA.r_drops**2,self.DHARMA.dN_drops)
        #den=np.trapz(den,self.DHARMA.r_drops,axis=0) dN_drops is dN, not n(r)
        re_tau=num/den
        Re_vw=c*np.einsum('zxy,zxy->xy',re_tau*w_tauDIVc,self.dTau[band,:])
        #Re_vw=c*np.trapz(re_tau*w_tauDIVc,tau,axis=0)
        #calculating ve_tau (Hanse & Travis 1974 2.54)
        rMinre=np.zeros_like(self.DHARMA.dN_drops,dtype=float)
        for i in np.arange(self.x.size):
            for j in np.arange(self.y.size):
                for k in np.arange(self.DHARMA.z.size):
                    rMinre[:,k,i,j]=self.DHARMA.r_drops-re_tau[k,i,j]
        num=np.einsum('r,rzxy->zxy',self.qe_avg[band,:]*self.alb_av[band,:]*self.DHARMA.r_drops**2,self.DHARMA.dN_drops*rMinre**2)
        #num=np.trapz(num,self.DHARMA.r_drops,axis=0)
        ve_tau=num/den/re_tau**2
        Ve_vw=c*np.einsum('zxy,zxy->xy',ve_tau*w_tauDIVc,self.dTau[band,:])
        #Ve_vw=c*np.trapz(ve_tau*w_tauDIVc,tau,axis=0)
        w_tau=np.zeros_like(w_tauDIVc,dtype=float)
        for i in np.arange(self.x.size):
            for j in np.arange(self.y.size):
                w_tau[:,i,j]=w_tauDIVc[:,i,j]*c[i,j]
        #computing vertically weighted number concentration
        dN_vw = -np.einsum('izxy,zxy,zxy->ixy',self.DHARMA.dN_drops,w_tau,self.dTau[band,:])
        #computing re_dN_vw, re from vertically weighted DSDs (Hansen & Travis 1974 2.53, Miller et. al. 2016)
        num = np.einsum('r,rxy->xy',self.qe_avg[band,:]*self.alb_av[band,:]*self.DHARMA.r_drops**3,dN_vw)
        #num = np.trapz(num,self.DHARMA.r_drops,axis=0)
        den = np.einsum('r,rxy->xy',self.qe_avg[band,:]*self.alb_av[band,:]*self.DHARMA.r_drops**2,dN_vw)
        #den = np.trapz(den,self.DHARMA.r_drops,axis=0)
        re_dN_vw = num/den
        #computing ve_dN_vw, ve from vertically weighted DSDs (Hnasen & Travis 1974 2.53, Miller et. al. 2016) 
        ve_dN_vw=np.zeros_like(re_dN_vw,dtype=float)
        for i in np.arange(self.x.size):
            for j in np.arange(self.y.size):
                rMinre = self.DHARMA.r_drops-re_dN_vw[i,j]
                num    = self.qe_avg[band,:]*self.alb_av[band,:]*self.DHARMA.r_drops**2*dN_vw[:,i,j]*rMinre**2
    #            num    = np.trapz(num,obj.DHARMA.r_drops)
                ve_dN_vw[i,j] = 1/re_dN_vw[i,j]**2*num.sum()/den[i,j]
        
        #calculating vertically weighted cloud optical thickness ?? (:D LOL)
        end=time.time()
        print('%0.2f mins elapsed!'%((end-start)/60))
        return Re_vw,Ve_vw,dN_vw,re_tau,ve_tau,w_tau,tau,re_dN_vw,ve_dN_vw
    def setup_reVW(self,mie_name,lesCname,dgSZA,dgVZA,a=1,b=0,lesBinEd=None,mie_path=None,band=None,fpath=None,\
                   replace=None):
        '''
        fpath: Path to saved psudo retrieval hdf5 (default:/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_simulations/)
        See DHARMA_onmp.find_reVW for I/O help.
        '''
        from vertical_weighting_psudo_ret import LES_psudo_rets
        VW=LES_psudo_rets(self,mie_name,mie_path,dgSZA,dgVZA,a,b,band=band,replace=replace)
        if VW.replace=='1':
            VW.Re,VW.Ve,VW.dN,VW.re_tau,VW.ve_tau,VW.w_tau,VW.tau,VW.Re_dN_vw,VW.Ve_dN_vw = self.find_reVW(mie_name,lesCname,dgSZA,dgVZA,a=a,b=b,mie_path=mie_path,band=band)  
            VW.Tau=VW.tau[0,:,:]
            VW.saveVW()
        else:
            VW.loadVW()
        self.VW=VW
    def findLWP(self,obj=None):
        '''
        Same as Dharma_LES_read.findLWP()
        Find LWP by integrating alog the each bin
        obj: Dharma_LES_read.DHARMA object
        return:
            lwc,lwp
        '''
        if obj==None:
            ATEX=self.DHARMA
        else:
            ATEX=obj
        rhol=1e6#g/m^3
        r=ATEX.r_drops#um 
        z=self.dz
        nr=ATEX.dN_drops
        l1=np.einsum('i,ijkl->jkl',r**3,nr)
        lwc=l1*rhol*4.0/3.0*np.pi*1e-12
        lwp=np.einsum('jkl,j->kl',l1,z)*rhol*np.pi*4.0/3*1e-12
        return lwc,lwp
    def find_cldTop(self,obj=None,th=1):
        '''
        Define cloud top based on the liquid water path threshold starting from the top
        th: lwp threshold in g/m^2. Default value is 1 g/m^2
        '''
        if obj==None:
            ATEX=self.DHARMA
        else:
            ATEX=obj
        rhol=1e6#g/m^3
        r=ATEX.r_drops#um 
        z=self.dz
        nr=ATEX.dN_drops
        l1=np.einsum('i,ijkl->jkl',r**3,nr)
        dlwp=np.einsum('jkl,j->jkl',l1,z)*rhol*np.pi*4.0/3*1e-12
        self.ctop=np.zeros((dlwp.shape[1],dlwp.shape[2]),dtype=float)
        for k in np.arange(0,dlwp.shape[1]):
            for l in np.arange(0,dlwp.shape[2]):
                for j in np.arange(0,dlwp.shape[0]):
                    if dlwp[:j,k,l].sum()>th:
                        self.ctop[k,l]=ATEX.z[j]
                        break
    def set_cot_top(self,band=1,ctop_tau=0.2):
        '''
        Setup self.cot_top 2D array with cloud top vertical indices of each column.
        Use self.DHARMA.z get the altitude.
        '''
        if hasattr(self,"cot_top"):
            print("cot_top already exist!! ")
        else:
            if hasattr(self,"dTau"):
                self.setup_zTau0top()
                tau = self.zTau0top
                _,z,x,y = tau.shape
                print("%0.3f band selected to compute COT-based cloud top"%self.c_mie.wvl[band])
                cot_top = np.zeros((x,y),dtype=int)
                for xi in range(x):
                    for yi in range(y):
                        min_tau = ctop_tau 
                        topi = 0
                        for zi in range(z):
                            if tau[band,zi,xi,yi] > min_tau:
                                min_tau = 0
                                topi = zi
                        cot_top[xi,yi] = topi        
                max_ctop = self.DHARMA.z[cot_top.max()]/1e3
                print("COT-based domain cloud top: %0.2f km"%max_ctop)
                self.cot_top = cot_top
                print("self.cot_top setup successfully")    
            else:
                print('self.dTau does not exist!!! Use self.setup_dTau() to setup')
        
    def find_Bext(self,obj):
        '''
        Find the Bext ( volime extinction coefficient area/vol)
        Bext[j]=Qe(r_i)*pi*r_i^2dN_ij ; Einstein summation (over i)
        obj: Dharma_LES_read.DHARMA object
        '''
        Qe=2.09
        nr=obj.dN_drops
        r=obj.r_drops
        Bext=np.einsum('i,ijkl->jkl',r**2,nr)*np.pi*Qe*1e-6#m^-1
        return Bext
    def add_cloud_mask(self,cldM):
        '''
        Mask 'cloud free' pixels
        cldM: 2D array. 'True' are cloudy pixels
		lwc,lwp,Bext,tauQ,dTau,Tau will be masked.
        '''
        if self.cloud_mask=='not_done':
            self.lwc[:,np.invert(cldM)]=np.nan
            self.lwp[np.invert(cldM)]=np.nan
            self.BextQe2p09[:,np.invert(cldM)]=np.nan
            self.tauQe2p09[np.invert(cldM)]=np.nan
            self.dTau[:,:,np.invert(cldM)]=np.nan
            self.Tau[:,np.invert(cldM)]=np.nan
            self.cloud_mask='done'
        else:
            print('Already applied a cloud mask!')
    

def iqu(fig,ax,ds,iqu,v_the,ttl=None,xy_extent=None):
#    xy_extent:[0,6.4,0,6.4] actual xy dimensions
    if ttl==None:
        ttl=['I','Q','U']
        ax.set_title(ttl[iqu])
    else:
        ax.set_title(ttl)
    if v_the==0:
        ctf=ax.imshow(np.mean(ds.MeanPRad[ds.VZA==v_the,:,:,iqu],axis=0),cmap='jet',origin='lower')
    else:
        ctf=ax.imshow(np.squeeze(ds.MeanPRad[ds.VZA==v_the,:,:,iqu]),cmap='jet',origin='lower')
    cb=fig.colorbar(ctf,ax=ax,orientation='horizontal',extend='both',pad=0.2)
    if xy_extent!=None:
        ctf.set_extent(xy_extent)
    return ctf,cb
def ip(fig,ax,ds,v_the,ttl=None,xy_extent=None):
    if ttl==None:
        ttl=r'$\sqrt(I^2+Q^2)$'
        ax.set_title(ttl)
    else:
        ax.set_title(ttl)
    if v_the==0:
        Q=np.mean(ds.MeanPRad[ds.VZA==v_the,:,:,1],axis=0)
        U=np.mean(ds.MeanPRad[ds.VZA==v_the,:,:,2],axis=0)
        Ip=np.sqrt(Q**2+U**2)
        ctf=ax.imshow(Ip,cmap='jet',origin='lower')
    else:
        Q=np.squeeze(ds.MeanPRad[ds.VZA==v_the,:,:,1])
        U=np.squeeze(ds.MeanPRad[ds.VZA==v_the,:,:,2])
        Ip=np.sqrt(Q**2+U**2)
        ctf=ax.imshow(Ip,cmap='jet',origin='lower')
    cb=fig.colorbar(ctf,ax=ax,orientation='horizontal',extend='both',pad=0.2)
    if xy_extent!=None:
        ctf.set_extent(xy_extent)
    return ctf,cb
def iquRMS(fig,ax,ds,iqu,v_the):
    if v_the==0:
        ctf=ax.imshow(np.mean(ds.RMSEPRad[ds.VZA==v_the,:,:,iqu],axis=0),cmap='jet',origin='lower')
    else:
        ctf=ax.imshow(np.squeeze(ds.RMSEPRad[ds.VZA==v_the,:,:,iqu]),cmap='jet',origin='lower')
    cb=fig.colorbar(ctf,ax=ax,orientation='horizontal',extend='both')
    return ctf,cb

def pltIQU(ds,v_the,sv=False):
    '''
    I think this is for LES cases. Check next time.
    ds: cpnLES_MSCARTlib.POLCARTdset
    v_the: Viewing Zenith Angle
    '''
    fig1,ax1=plt.subplots(2,3,figsize=(8,7))
    fig1.suptitle("\n".join(wrap(ds.fname+'(%0.0f mins)'%ds.MeanTiming+\
                ' * %0.0E'%ds.Nbat+'_vthe%d'%v_the,50)))
    fig1.subplots_adjust(top=0.85)
    #fig1.subplots_adjust(hspace=0.5)

    iqurow(fig1,ax1,0,ds,v_the)
    iqurow(fig1,ax1,1,ds,v_the,rms=True)
    
    fig1.show()

    if sv:
        fig1.savefig('figures/'+ds.fname+'_IQU_vthe%d'%v_the+'.png',format='png',dpi=200)
        print('figures/'+ds.fname+'_IQU_vthe%d'%v_the+'.png'+' SAVED!')
        
def iqurow(fig,ax,row,ds,VZA,rms=False):
    if rms:
        ctf,cb=iquRMS(fig,ax[row,0],ds,0,VZA)
        ctf.set_clim(0.001,0.003)
        cb.set_ticks(np.linspace(0.001,0.003,2))
        ctf,cb=iquRMS(fig,ax[row,1],ds,1,VZA)
        ctf.set_clim(0.0,0.001)
        cb.set_ticks(np.linspace(0,0.001,2))
        ctf,cb=iquRMS(fig,ax[row,2],ds,2,VZA)
        ctf.set_clim(0,0.001)
        cb.set_ticks(np.linspace(0,0.001,2))
    else:
        ctf,cb=iqu(fig,ax[row,0],ds,0,VZA)
        ctf.set_clim(0,0.8)
        cb.set_ticks(np.arange(0,1.1,0.25))
        ctf,cb=iqu(fig,ax[row,1],ds,1,VZA)
        ctf.set_clim(-0.020,-0.001)
        cb.set_ticks(np.arange(-0.015,0.001,0.010))
        ctf,cb=iqu(fig,ax[row,2],ds,2,VZA)
        ctf.set_clim(-0.02,0.001)
        cb.set_ticks(np.arange(-0.015,0.001,0.010))
    ax[row,0].scatter([90],[100],c='w')
    ax[row,0].scatter([40],[40],c='k')
    ax[row,1].scatter([90],[100],c='w')
    ax[row,1].scatter([40],[40],c='k')
    ax[row,2].scatter([90],[100],c='w')
    ax[row,2].scatter([40],[40],c='k')

def tick_fun(VZA,SZA,SAA):
    V=scat_ang(SZA,VZA,SAA)
    return ["%d"%z for z in V]
