#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Sat Sep 16 08:55:26 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Retrieval codes to do NJK and polarimetric retrievals
01/21/2017:
    Color ratio method (Jethva 2013) CRLUT added(not finalized)
"""
import numpy as np
import netCDF4, pickle
import matplotlib.pyplot as plt
import h5py, time, os
import scipy.interpolate as itp
from scipy.optimize import curve_fit
import cpnCommonlib as cpn
from cpnMielib import MieSet, bulk_Mie, PSDs

cpn.setup_figures(plt)
def interp(x,xp,yp):
    if xp[0]>xp[-1]:
        xp=np.flipud(xp)
        yp=np.flipud(yp)
    y=np.interp(x,xp,yp)
    
    return y 

#To read these type of Dan's LUT data.
#fname='MODIS_LUT_extended_SZA%03d_RAA000.nc'%(180-int(sza))
class Bispec_LUT(object):
    def __init__(self,fdpath,fname):
        self.fdpath=fdpath
        self.fname=fname
        self.band=[]
        self.re=[]
        self.ve=[]
        self.tau=[]
        self.mu0=[]
        self.phi=[]
        self.mu=[]
        self.scat=[]
        self.I=[]
        self.Q=[]
        self.U=[]
    def readLUT(self,):
        data=netCDF4.Dataset(self.fdpath+self.fname,'r')
        self.band=data.variables['band'][:]
        self.re=data.variables['re'][:]
        self.ve=data.variables['ve'][:]
        self.tau=data.variables['tau'][:]
        self.mu0=data.variables['mu0'][:]
        self.phi=data.variables['phi'][:]
        self.mu=data.variables['mu'][:]
        self.scat=data.variables['scat'][:]
        self.I=data.variables['I'][:]
        self.Q=data.variables['Q'][:]
        self.U=data.variables['U'][:]
    def find_mu_ix(self,VZA,pro_type='original'):
        '''
        Return appropriate VZA index of the objects data set.
        VZA: Viewing zenith angle in degrees
        pro_type: Product type
            - 'original' implies the first set of files
            - 'DISORT' implies the second set of LUT generations that I did after the defense.
            - 'MSCART' implies the second set of LUT generations that I did after the defense.
        '''
        if VZA==-20.0 and pro_type=='original':
            return 107
        elif VZA==0.0 and pro_type=='original':
            return 149
        elif VZA==0.0 and (pro_type=='DISORT' or pro_type='MSCART'):
            return 0
        else:
            mu_VZA=np.cos(np.deg2rad(VZA))
            mu_ix=np.where(np.isclose(self.mu,mu_VZA,atol=0.01))
            print('%%%%VZA indices for LUT%%%%')
            print('VZA indices:'+str(mu_ix))
            print('VZA:'+str(np.rad2deg(np.arccos(self.mu[mu_ix]))))
            mu_ix=input('Select VZA index:')
            return mu_ix
    def plotLUT(self,VZA,re_lines=None,ve=0.02,figAx=None,gcolor='grey',pro_type='original'):
        '''
        re_line: Array of re values that are needed to be shown on the LUT plot.
                Make sure all the values are available in self.re
        figAx: (fig1,ax1)
        pro_type: Product type
            - 'original' implies the first set of files
            - 'DISORT' implies the second set of LUT generations that I did after the defense.
        '''
        #plot LUT using 0.865 and 2.13
        mu_ix=self.find_mu_ix(VZA,pro_type=pro_type)
        LUT_VNIR=self.I[0,:,np.where(self.ve==ve),:,mu_ix].T
        LUT_SWIR=self.I[1,:,np.where(self.ve==ve),:,mu_ix].T
        if re_lines is None:
            fig_LUT,axLUT=plotLUT(LUT_VNIR,LUT_SWIR,self.re[:],self.tau,figAx=figAx,gcolor=gcolor)
        else:        
            fig_LUT,axLUT=plotLUT(LUT_VNIR,LUT_SWIR,self.re,self.tau,re_lines=re_lines,figAx=figAx,gcolor=gcolor)
        axLUT.set_title(r'$\mu$=%0.3f : $\mu_0$=%0.3f : $\phi$=%0.1f : $v_e$=%0.2f'%(self.mu[mu_ix],self.mu0,self.phi,ve))
        return fig_LUT,axLUT


def plotLUT(VNIR_lut,SWIR_lut,reff_lut,tau_lut,re_lines=None,figAx=None,gcolor='grey'):
    '''
    To plot LUT. fig and ax instances will be returned to over plot data.
    re_lines:  Array of re values that are needed to be shown on the LUT plot.
                Make sure all the values are available in self.re
    '''
    if figAx is None:
        fig_LUT,axLUT=plt.subplots()
    else:
        fig_LUT,axLUT=figAx[0],figAx[1]
    if re_lines is None:
        for i in np.arange(0,reff_lut.size,10):
            x=np.squeeze(VNIR_lut[i,:]);y=np.squeeze(SWIR_lut[i,:])
            axLUT.plot(x,y,'k-',color=gcolor)
            axLUT.text(x.max(),y.max(),' %d'%reff_lut[i],color=gcolor)
        re_min_ix=np.argmin(reff_lut)
    else:
        p=np.repeat(re_lines,reff_lut.size).reshape(re_lines.size,reff_lut.size)
        q=np.repeat(reff_lut,re_lines.size).reshape(reff_lut.size,re_lines.size)
        re_ln_ix=np.where((p-q.T)==0)[1]
        for i in re_ln_ix:
            x=np.squeeze(VNIR_lut[i,:]);y=np.squeeze(SWIR_lut[i,:])
            axLUT.plot(x,y,'k-',color=gcolor)
            axLUT.text(x.max(),y.max(),' %d'%reff_lut[i],color=gcolor)
        re_min_ix=np.squeeze(np.argwhere(reff_lut==re_lines.min()))
    for j in np.arange(0,tau_lut.size,10):
        x=np.squeeze(VNIR_lut[re_min_ix:,j]);y=np.squeeze(SWIR_lut[re_min_ix:,j])
        axLUT.plot(x,y,'k--',color=gcolor)
        axLUT.text(x.min(),y.min()-0.02,r'%d'%tau_lut[j],color=gcolor)
        axLUT.set_xlabel('VNIR Reflectance')
        axLUT.set_ylabel('SWIR Reflectance')
    axLUT.set_xlim(0,VNIR_lut.max()+0.2)
    axLUT.arrow(VNIR_lut.max()+.1,0.5,0.03,-0.4,linestyle='-',head_width=0.02,\
              length_includes_head=True,color=gcolor)
    axLUT.annotate(r' $r_e$($\mu m$)',xy=(VNIR_lut.max()+0.12,0.3),color=gcolor,rotation=90)
    axLUT.arrow(0.2,0.05,VNIR_lut.max()-0.3,0,linestyle='-',head_width=0.02,\
              length_includes_head=True,color=gcolor)
    axLUT.annotate(r' $\tau$',xy=(VNIR_lut.max()/2,0),color=gcolor)
    fig_LUT.show()
    return fig_LUT,axLUT


#An object which can keep NJK retrieval data and also can do NJK retrievals by 
#using my "retrieve_NJK" function.
class NJK_retrievals(object):
    def __init__(self,RT865_file,RT2p13_file,Physics_file,LUT_file):
        #Data files related to this retrievals. (Only for information purposes)
        self.RT865_file=RT865_file
        self.RT2p13_file=RT2p13_file
        self.Physics_file=Physics_file
        self.LUT_file=LUT_file
    def set_NJK_inputs(self,):
        self.VNIR_lut=[]#[n_reff_lut,n_tau_lut]
        self.SWIR_lut=[]#[n_reff_lut,n_tau_lut]
        self.reff_lut=[]
        self.tau_lut=[]
        self.obs_VNIR_set=[]
        self.obs_SWIR_set=[]
    def set_case_name(self,):
        '''
        To give a name to the retrieval case.
        '''
        les_case='_'.join(self.Physics_file.split('/')[-1].split('_')[:-1])
        bands=self.RT2p13_file.split('/')[-1].split('_')[3]+'_'+\
            self.RT865_file.split('/')[-1].split('_')[3]
        apnd='_'.join(self.RT2p13_file.split('/')[-1].split('_')[4:]).split('.')[0]
        self.NJK_case_name='NJK_retrievals_'+les_case+'_'+bands+'_'+apnd    

    def NJK_Python(self,dim='1D'):
        '''
        NJK retrievals (Python code)
        ------------------------------
        dim=1:1 or 2 for 1D or 2D obs_VNIR_set array dimensions 
        '''
        if dim=='1D':
            print('Executing NJK retrievals for the given 1D field')
            self.NJK_Python_1D()
        if dim=='2D':
            print('Executing NJK retrievals for the given 2D field')
            st=time.time()
            self.NJK_Python_2D()
            ed=time.time()
            print('%d Elapsed!'%(ed-st))
        
    def NJK_Python_1D(self,):
#        from cpnretrieve_NJK import retrieve_NJK
        if not(hasattr(self,'VNIR_lut')):
            print('Set all inputs as in set_NJK_inputs first !')
        else:
            re=np.zeros((self.obs_VNIR_set.size),dtype=float)
            tau=np.zeros_like(re)
            flag=np.zeros_like(re)
            print('Running NJK_Python ... ')
            for i in np.arange(0,self.obs_VNIR_set.size):
                re[i],tau[i],flag[i]=retrieve_NJK(self.obs_VNIR_set[i],self.obs_SWIR_set[i],\
                self.VNIR_lut,self.SWIR_lut,self.reff_lut,self.tau_lut)
            self.re=re
            self.tau=tau
            self.flag=flag
            print('Bispectral retrievals done!')
    def NJK_Python_2D(self,):
#        from cpnretrieve_NJK import retrieve_NJK
        if not(hasattr(self,'VNIR_lut')):
            print('Set all inputs as in set_NJK_inputs first !')
        else:
            re=np.zeros_like(self.obs_VNIR_set,dtype=float)
            tau=np.zeros_like(re,dtype=float)
            flag=np.zeros_like(re,dtype=int)
            print('Running NJK_Python ... ')
            for i in np.arange(0,self.obs_VNIR_set.shape[0]):
                for j in np.arange(0,self.obs_VNIR_set.shape[1]):
                    re[i,j],tau[i,j],flag[i,j]=retrieve_NJK(self.obs_VNIR_set[i,j],self.obs_SWIR_set[i,j],\
                    self.VNIR_lut,self.SWIR_lut,self.reff_lut,self.tau_lut)
            self.re=re
            self.tau=tau
            self.flag=flag
            print('Bispectral retrievals done!')

def doNJK_LES(LES_case_VNIR,LES_case_SWIR,sza,VZA,RTdim='3D',check_DB=True):
    '''
    To do NJK retrievals of LES cases. 1D or 3D
    ---------------------------------------------------------
    LES_case_VNIR:cpnLES_MSCARTlib.LES_case object for VNIR
    LES_case_SWIR:cpnLES_MSCARTlib.LES_case object for SWIR
    sza: Source zenith angle string ('140')
    VZA: Viewing zenith angle float
    RTdim: '3D' or '1D' RT simulation instance (RT_1D). 
        For 1D case, both LES object must have RT1D instance.
    check_DB: True check data directory for previous retrieved data files
        (/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Retrievals/NJK_retrievals/data/)
    return: 
        cpnRetrievalslib.NJK_retrievals object
    '''
    fdpath='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_retrievals/LUTs/'
    fname='MODIS_LUT_extended_SZA%03d_RAA000.nc'%(180-int(sza))
    LUT=Bispec_LUT(fdpath,fname)
    LUT.readLUT()
    mu_ix=LUT.find_mu_ix(VZA)
    LUT_VNIR=LUT.I[0,:,1,:,mu_ix].T
    LUT_SWIR=LUT.I[1,:,1,:,mu_ix].T

    obVZA_ix=LES_case_VNIR.RT.find_obVZA_ix(VZA)
    if RTdim=='3D':
        RT865_file=LES_case_VNIR.get_file_name(LES_case_VNIR.name+'_3D',res=LES_case_VNIR.res)
        RT2p13_file=LES_case_SWIR.get_file_name(LES_case_SWIR.name+'_3D',res=LES_case_SWIR.res)
        obs_SWIR_set=LES_case_SWIR.RT.MeanPRad[obVZA_ix,:,:,0]
        obs_VNIR_set=LES_case_VNIR.RT.MeanPRad[obVZA_ix,:,:,0]
    elif RTdim=='1D':
        RT865_file=LES_case_VNIR.get_file_name(LES_case_VNIR.name+'_1D',res=LES_case_VNIR.res)
        RT2p13_file=LES_case_SWIR.get_file_name(LES_case_SWIR.name+'_1D',res=LES_case_SWIR.res)
        obs_SWIR_set=LES_case_SWIR.RT1D.MeanPRad[obVZA_ix,:,:,0]
        obs_VNIR_set=LES_case_VNIR.RT1D.MeanPRad[obVZA_ix,:,:,0]
    Physics_file=LES_case_VNIR.RT_field.dpath+LES_case_VNIR.RT_field.fname
    LUT_file='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_retrievals/LUTs/'+'MODIS_LUT_extended_SZA%03d_RAA000.nc'%(180-int(sza))
    
    NJK_ret=NJK_retrievals(RT865_file,RT2p13_file,Physics_file,LUT_file)
    NJK_ret.VNIR_lut=LUT_VNIR
    NJK_ret.SWIR_lut=LUT_SWIR
    NJK_ret.reff_lut=LUT.re
    NJK_ret.tau_lut=LUT.tau
    NJK_ret.obs_SWIR_set=obs_SWIR_set
    NJK_ret.obs_VNIR_set=obs_VNIR_set
    NJK_ret.set_case_name()
    if check_DB and os.path.isfile('/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Retrievals/NJK_retrievals/data/'+\
                   NJK_ret.NJK_case_name+'.pkl'):
        print('Retrievals already exist!: '+'/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Retrievals/NJK_retrievals/data/'+\
                   NJK_ret.NJK_case_name+'.pkl')
        logic=input('Do you want to continue? y/n:')
    else:
        logic='y'
    if logic=='y':
        NJK_ret.NJK_Python('2D')
    return NJK_ret
                
class NJKobj(object):
    '''
    d3D=NJKobj('NJK_retrievals/',"fractal_cld_re12ve05_x40km_MSCART_SZA"+sza+"_VZA"+NJK_vza+".pkl")
    කියවන්ට NJK retrievals ඔක්කොම ඕනවෙන දේවල් එක්කම.
    '''
    def __init__(self,ret_path,fname,resl=None):
        self.ret_path=ret_path
        self.fname=fname
        if resl==None:
            self.read()
        else:
            self.res=resl
    def read(self,):
        with open(self.ret_path+self.fname, "rb") as f:
            self.NJK_ret=pickle.load(f)
            self.d865SZAx=pickle.load(f)
            self.d2p13SZAx=pickle.load(f)
            self.dstep=pickle.load(f)
        
def read_agg_NJK_ret(d_native,res):
    #res: resolution 10,20,100,...
    #NJK_native: NJKobj with native resolution
    #return: NJKobj
    fname=d_native.ret_path+d_native.fname.split('.',1)[0]+'_%dm.pkl'%res
    with open(fname,"rv") as f:
        NJK_ret=pickle.load(f)
    return NJK_ret

def retrieve_NJK(VNIR,SWIR,VNIR_lut,SWIR_lut,reff_lut,tau_lut):
    '''
    Bispectra retrieval code
    Note that LUT lines have to be monotonic.

    VNIR_lut: shorter wavelength(0.865) non-absorbing band
    SWIR_lut: longer wavelength(2.13) absorbing band

    flag=2 : Reflectances are out from the LUT space
    flag=3 : re and tau guesses are out from the LUT space
    flag=4 : re and tau retrievals are out from the LUT space
    flag=5 : last two gusses have > 1e-7 difference
    '''
    reff_lut=reff_lut.reshape(1,reff_lut.size)
    tau_lut=tau_lut.reshape(1,tau_lut.size)
    #If we are outside the LUT
    if VNIR < VNIR_lut.min() or VNIR > VNIR_lut.max() or  SWIR<SWIR_lut.min() or SWIR>SWIR_lut.max():
        flag=2
        tau_guess=np.nan
        reff_guess=np.nan
    else:
        tau_matrix=np.repeat(tau_lut,reff_lut.size,0)
        reff_matrix=np.repeat(reff_lut,tau_lut.size,0).T
    
        #cost fnction
        cost_function=(VNIR_lut-VNIR)**2/(VNIR**2)+(SWIR_lut-SWIR)**2/(SWIR**2)
        p0=np.where(cost_function==cost_function.min())
        if np.size(p0)>2: #When two minimums found
            p0=(p0[0][0],p0[1][0])
        tau_guess=tau_matrix[p0]
        reff_guess=reff_matrix[p0]
        
        for i in np.arange(0,184):
#            print('re_guess:'+str(reff_guess))
            dummy_tau=tau_guess;dummy_reff=reff_guess
            int_refl=np.squeeze(np.zeros_like(tau_lut,dtype=float))
            for j in np.arange(0,tau_lut.size):
                int_refl[j]=interp(reff_guess,np.squeeze(reff_lut),VNIR_lut[:,j])
            tau_guess=interp(VNIR,int_refl,np.squeeze(tau_lut))
            
            int_refl=np.squeeze(np.zeros_like(reff_lut,dtype=float))
            for k in np.arange(0,reff_lut.size):
                int_refl[k]=interp(tau_guess,np.squeeze(tau_lut),SWIR_lut[k,:])
            reff_guess=interp(SWIR,int_refl,np.squeeze(reff_lut))
#            if SWIR>0.25:
#                print('here')
            
            if tau_guess < np.nanmin(tau_lut) or tau_guess > np.nanmax(tau_lut) or \
                    reff_guess < np.nanmin(reff_lut) or reff_guess > np.nanmax(reff_lut):
                        flag=3 ; break
            if abs(dummy_tau-tau_guess) < 1e-7 and abs(dummy_reff-reff_guess) < 1e-7:
                        break
        if dummy_tau-tau_guess > 1e-7 or dummy_reff-reff_guess > 1e-7:
            flag=5 #last two guesses has > 1e-7 difference
                    
    ret_tau=tau_guess;ret_re=reff_guess
    if ret_tau < np.nanmin(tau_lut) or ret_tau > np.nanmax(tau_lut) or \
        ret_re < np.nanmin(reff_lut) or ret_re > np.nanmax(reff_lut) :
            ret_tau=np.nan
            ret_re=np.nan
            flag=4
    elif ~np.isnan(ret_tau) and ~np.isnan(ret_re):
        flag=1

    return ret_re,ret_tau,flag

class CRLUT(object):
    '''
    Color ratio method LUT object
    '''
    def __init__(self,LUT_name):
        self.fname=LUT_name
    def setup_LUT(self,AOT,COT,bnd,SZA,LUT,CR,VZA):
        '''
        setting up LUT elements AOT,COT,bnd,SZA,
        LUT: [AOT,COT,band,SZA,VZA/ScatA]
        CR: [AOT,COT,SZA,VZA/ScatA]
        '''
        self.AOT=AOT
        self.COT=COT
        self.bnd=bnd
        self.SZA=SZA
        self.LUT=LUT
        self.CR=CR
        self.VZA=VZA
    def save_LUT(self):
        '''
        Saving as a hdf5 file
        '''
        f=h5py.File(self.fname,'w')
        f.attrs['object']='CRLUT'
        f.attrs['Action']='color_ratio_method_Jethva_2013'
        
        PCentry=f.create_dataset('LUT',data=self.LUT)
        PCentry.dims[0].label='AOT'
        PCentry.dims[1].label='COT'
        PCentry.dims[2].label='bnd'
        PCentry.dims[3].label='SZA'
        PCentry.dims[4].label='VZA'
        PCentry.attrs['units']='Reflectance'
        PCentry.attrs["long_name"]='Reflectance_LUT'
        
        PCentry=f.create_dataset('CR',data=self.CR)
        PCentry.dims[0].label='AOT'
        PCentry.dims[1].label='COT'
        PCentry.dims[2].label='SZA'
        PCentry.dims[3].label='VZA_ScatA'
        PCentry.attrs['units']='none'
        PCentry.attrs["long_name"]='Color_ratio_470_div_860'
        
        PCentry=f.create_dataset('AOT',data=self.AOT)
        PCentry.attrs['units']='none'
        PCentry.attrs['long_name']='Aerosol_Optical_depth_at_860'
        
        PCentry=f.create_dataset('COT',data=self.COT)
        PCentry.attrs['units']='none'
        PCentry.attrs['long_name']='Cloud_Optical_depth_at_860'
        
        PCentry=f.create_dataset('bnd',data=self.bnd)
        PCentry.attrs['units']='mu'
        PCentry.attrs['long_name']='band'
        
        PCentry=f.create_dataset('SZA',data=self.SZA)
        PCentry.attrs['units']='degrees'
        PCentry.attrs['long_name']='Solar_Zenith_Angle'
        
        PCentry=f.create_dataset('VZA',data=self.VZA)
        PCentry.attrs['units']='degrees'
        PCentry.attrs['long_name']='Viewing_Zenith_Angle'

        f.close()
        print(self.fname+' SAVED!')
    def read_LUT(self,fdpath=None):
        '''
        Reading hdf5 file
        '''
        if fdpath==None:
            f=h5py.File(self.fname,'r')
        else:
            f=h5py.File(fdpath+self.fname,'r')
        self.LUT=f['LUT'][:]
        self.CR=f['CR'][:]
        self.AOT=f['AOT'][:]
        self.COT=f['COT'][:]
        self.bnd=f['bnd'][:]
        self.SZA=f['SZA'][:]
        self.VZA=f['VZA'][:]

    def plt_LUT(self,vza,sza):
        '''
        Plot the LUT of given vza and sza
        return:
            fig1,ax1,fig1_ttl
        '''
        fig1,ax1=plt.subplots()
        fig1_ttl=self.fname.split('.',1)[0]+'_LUT_SZA%03d_VZA%03d'%(sza,vza)
        fig1.suptitle(fig1_ttl)
        for i in np.arange(0,self.COT.size):
            x=np.mean(self.LUT[:,i,self.bnd==0.860,self.SZA==sza,self.VZA==vza],axis=1)
            y=np.mean(self.CR[:,i,self.SZA==sza,self.VZA==vza],axis=1)
            if i==0:
                ax1.plot(x,y,'k.-',label='COT')
            else:
                ax1.plot(x,y,'k.-',label='_nolegend_')
            ax1.text(np.nanmax(x),np.nanmax(y),'%d'%self.COT[i],size=12)
        for i in np.arange(0,self.AOT.size):
            x=np.mean(self.LUT[i,:,self.bnd==0.860,self.SZA==sza,self.VZA==vza],axis=0)
            y=np.mean(self.CR[i,:,self.SZA==sza,self.VZA==vza],axis=0)
            if i==0:
                ax1.plot(x,y,'k.--',label='AOT')
            else:
                ax1.plot(x,y,'k.--',label='_nolegend_')
            ax1.text(np.nanmax(x),np.nanmin(y),'%0.2f'%self.AOT[i])
        ax1.set_xlabel(r'$R_{0.860}$')
        ax1.set_ylabel(r'$R_{0.470}/R_{0.860}$')
        ax1.legend(loc='best')
        fig1.show()
        return fig1,ax1,fig1_ttl

def retrieve_CR(VNIR,CR,VNIR_lut,CR_lut,AOT_lut,COT_lut):
    '''
    Color ratio method of Jethva et al. 2013. 
    Same as NJK. So uses retrieve_NJK()
    
    b865: VNIR reflectance (float)
    CR: color ration (float)
    AOT_lut: AOT array of the LUT space (1D array)
    COT_lut: COT array of the LUT space (1D array)
    b865_lut: VNIR (0.865) band LUT reflectance array (AOT_lut.size,COT_lut.size)
    CR_lut: CR(.470/.860) LUT array (AOT_lut.size,COT_lut.size)
    ----------------------------------------------------------------------
    return
    ret_AOT:
    ret_COT:
    flats:
        flag=1 : Successful retrievals
        flag=2 : Reflectances are out from the LUT space
        flag=3 : AOT and COT guesses are out from the LUT space
        flag=4 : AOT and COT retrievals are out from the LUT space
        flag=5 : last two gusses have > 1e-7 difference
    '''
    VNIR=VNIR
    SWIR=CR
    VNIR_lut=VNIR_lut
    SWIR_lut=CR_lut
    reff_lut=AOT_lut
    tau_lut=COT_lut
    ret_AOT,ret_COT,flag=retrieve_NJK(VNIR,SWIR,VNIR_lut,SWIR_lut,reff_lut,tau_lut)
    return ret_AOT,ret_COT,flag
    

def retrieve_CR2(b865,CR,b865_lut,CR_lut,AOT_lut,COT_lut):
    '''
    Don't use this. Almost same as retrieve_CR()
    Color ratio method of Jethva et al. 2013.
    b865: VNIR reflectance (float)
    CR: color ration (float)
    AOT_lut: AOT array of the LUT space (1D array)
    COT_lut: COT array of the LUT space (1D array)
    b865_lut: VNIR (0.865) band LUT reflectance array (AOT_lut.size,COT_lut.size)
    CR_lut: CR(.470/.860) LUT array (AOT_lut.size,COT_lut.size)
    ----------------------------------------------------------------------
    return
    ret_AOT:
    ret_COT:
    flag:
        #flag=1 : Sucessful retreival
        #flag=2 : Reflectances are out from the LUT space
        #flag=3 : re and tau guesses are out from the LUT space
        #flag=4 : re and tau retrievals are out from the LUT space
    '''
    AOT_lut=AOT_lut.reshape(1,AOT_lut.size)
    COT_lut=COT_lut.reshape(1,COT_lut.size)
    #If we are outside the LUT
    if b865 < b865_lut.min() or b865 > b865_lut.max() or  CR<CR_lut.min() or CR>CR_lut.max():
        flag=2
        COT_guess=np.nan
        AOT_guess=np.nan
    else:
        tau_matrix=np.repeat(COT_lut,AOT_lut.size,0)
        reff_matrix=np.repeat(AOT_lut,COT_lut.size,0).T
    
        #cost fnction
        cost_function=(b865_lut-b865)**2/(b865**2)+(CR_lut-CR)**2/(CR**2)
        p0=np.where(cost_function==cost_function.min())
        COT_guess=tau_matrix[p0]
        AOT_guess=reff_matrix[p0]
        
        for i in np.arange(0,184):
#            print('re_guess:'+str(AOT_guess))
            dummy_tau=COT_guess;dummy_reff=AOT_guess
            int_refl=np.squeeze(np.zeros_like(COT_lut,dtype=float))
            for j in np.arange(0,COT_lut.size):
                int_refl[j]=interp(AOT_guess,np.squeeze(AOT_lut),b865_lut[:,j])
            COT_guess=interp(b865,int_refl,np.squeeze(COT_lut))
            
            int_refl=np.squeeze(np.zeros_like(AOT_lut,dtype=float))
            for k in np.arange(0,AOT_lut.size):
                int_refl[k]=interp(COT_guess,np.squeeze(COT_lut),CR_lut[k,:])
            AOT_guess=interp(CR,int_refl,np.squeeze(AOT_lut))
#            if CR>0.25:
#                print('here')
            
            if COT_guess < np.nanmin(COT_lut) or COT_guess > np.nanmax(COT_lut) or \
                    AOT_guess < np.nanmin(AOT_lut) or AOT_guess > np.nanmax(AOT_lut):
                        flag=3 ; break
            if abs(dummy_tau-COT_guess) < 1e-7 and abs(dummy_reff-AOT_guess) < 1e-7:
                        break
                    
    ret_COT=COT_guess;ret_AOT=AOT_guess
    if ret_COT < np.nanmin(COT_lut) or ret_COT > np.nanmax(COT_lut) or \
        ret_AOT < np.nanmin(AOT_lut) or ret_AOT > np.nanmax(AOT_lut):
            ret_COT=np.nan
            ret_AOT=np.nan
            flag=4
    else:
        flag=1

    return ret_AOT,ret_COT,flag

class Pol_ret_P12Lib(object):
    '''
    Polarimetric cloud retrievals P12 (theoretical) object
    '''
    def __init__(self,fname='Pol_ret_PM_library_0p860_2p13.hdf5',path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Retrievals/Pol_retrievals/'):
        '''
        fname:
            Pol_ret_PM_library_0p860_2p13_V2.hdf5 LUT ve has special features. Finer resolutions from 0.01-0.05.
                After that 0.025 steps upto 0.3.
        '''
        if not(path==None):
            self.fname=path+fname
        else:
            self.fname=fname
    def generateLib(self,re=np.arange(2,30,0.25),ve=np.arange(0.01,0.3,0.01)):
        '''
        To generate LUT library for polarimetric cloud property retrievals
        '''
        cloud_mie=MieSet('Cloud_mie_470_860_2p13',path='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_retrievals/polarimetric_cloud/')
        cloud_mie.readMie()
        psd=PSDs('mod_gamma',np.asarray(cloud_mie.d,dtype=float))
        self.re=re
        self.ve=ve
        self.avP12={'0p860':[],'2p13':[]}
        wli={'0p860':1,'2p13':2}
        for wl in self.avP12.keys():
            avP12_x=np.zeros((re.size,ve.size,cloud_mie.ang.size),dtype=float)    
            for i in np.arange(0,re.size):
                for j in np.arange(0,ve.size):
                    psd.set_mod_gamma_norm(re[i],ve[j])
                    bulk=bulk_Mie('bulk_Mie_470_860_2p13',psd=psd,Mie=cloud_mie)
                    bulk.cal_bulk_Pmat()
                    avP12_x[i,j,:]=bulk.bulk_P12[:,wli[wl]]
            print(wl+':%0.2f'%(bulk.Mie.wvl[wli[wl]]))
            self.avP12[wl]=avP12_x
        self.bulk_Mie_ang=bulk.Mie.ang
        self.bulk_Mie=bulk
    def saveP12Lib(self,):
        if os.path.isfile(self.fname):
            inp=raw_input('File aready exist. Continue? y/n:')
        else:
            inp='y'
        if inp=='y':
            f=h5py.File(self.fname,'w')
            entry=f.create_dataset('P12_0p860',data=self.avP12['0p860'])
            entry.dims[0].label='re'
            entry.dims[1].label='ve'
            entry.dims[2].label='angle'
            entry.attrs['long_name']='Polarized Phase function (bulk)'
        
            entry=f.create_dataset('P12_2p13',data=self.avP12['2p13'])
            entry.dims[0].label='re'
            entry.dims[1].label='ve'
            entry.dims[2].label='angle'
            entry.attrs['long_name']='Polarized Phase function (bulk)'
            
            entry=f.create_dataset('re',data=self.re)
            entry.attrs['units']='mu'
            entry.attrs['long_name']='Cloue Effective Radius'
            entry=f.create_dataset('ve',data=self.ve)
            entry.attrs['units']='None'
            entry.attrs['long_name']='Cloud effective variance'
            entry=f.create_dataset('angle',data=self.bulk_Mie_ang)
            entry.attrs['units']='degrees'
            entry.attrs['long_name']='Scattering angle'
            f.close()
            print(self.fname+' SAVED!')
    def loadP12Lib(self,):
        f=h5py.File(self.fname,'r')
        self.avP12={'op860':[],'2p13':[]}
        self.avP12['0p860']=f['P12_0p860'][:]
        self.avP12['2p13'] =f['P12_2p13' ][:]
        self.bulk_Mie_ang   =f['angle'][:]
        self.re = f['re'][:]
        self.ve = f['ve'][:]
        f.close()

def retrieve_pol_veFixd(P12_in,re_in,MieA,Q_in,ScatA,VZA,SZA):
    '''
    Breon et al. (1998) polarimetric retrievals (in my way)
    Exploits only the pattern of the polarize reflection AFTER normalizing to the clodbow value.
        (It seems normalization to the cloudbow value is not important. Only viewing gemetry correction would be enough)
    Uses correlation coefficients.
    Seems nice!
    -----------------------------------------------------
    P12_in: 2D Array [re_in,MieA] P12 from Mie calculations
    re_in: 1D Array
    MieA: 1D array scattering angle
    Q_in: 3D/2D array (single column/2D image) [ScatA,nxdim,nydim] or [ScatA]
    VZA: Viewing Zenith
    SZA: Solar Zenith
    return:
        re_ret: if Q_in[ScatA,nxdim,nydim] -> [nxdim,nydim]
    '''
    muS=np.cos(np.deg2rad(SZA))
    muV=abs(np.cos(np.deg2rad(VZA)))
    p12_a1=np.argwhere(MieA==130)
    p12_a2=np.argwhere(MieA==160)
    Q_a1=np.argmin(abs(ScatA-130))
    Q_a2=np.argmin(abs(ScatA-160))
    if Q_a1>Q_a2:
        Q_a1,Q_a2=Q_a2,Q_a1
    #P12 normalization
    cbow_P12=np.max(-P12_in[:,p12_a1:p12_a2],axis=1)#maximum -P12 at the cloud bow
    mNo_avP12=np.einsum('ij,i->ij',-P12_in,1/cbow_P12)#avP12 normalized by the maximum at the cloudbow
    
    #interpolation
    theP12=MieA[p12_a1-1:p12_a2+2]
    theQ=ScatA[Q_a1:Q_a2]
    mNo_P12=mNo_avP12[:,p12_a1-1:p12_a2+2]
    mNo_P12_itp=np.zeros((re_in.size,theQ.size),dtype='float')
    for i in np.arange(0,re_in.size):
        itpF=itp.interp1d(theP12,mNo_P12[i,:])#interpolated avP12 to match with observed scattering angles
        mNo_P12_itp[i,:]=itpF(theQ)
    
    #Q_in normalization
    if Q_in.ndim==3:
        Qdom=np.einsum('ijk,i->ijk',-Q_in,4*(muS+muV))    
        nom_Q=Qdom[Q_a1:Q_a2,:,:]/np.max(Qdom[Q_a1:Q_a2,:,:],axis=0)
        ret_re=np.zeros((nom_Q.shape[1],nom_Q.shape[2]),dtype=float)
        for i in np.arange(nom_Q.shape[1]):
            for j in np.arange(nom_Q.shape[2]):
                ret_re[i,j]=re_in[findbestPatveFixd(mNo_P12_itp,nom_Q[:,i,j],re_in)]
    elif Q_in.ndim==1:
        Qdom=-Q_in*4*(muS+muV)
        nom_Q=Qdom[Q_a1:Q_a2]/np.max(Qdom[Q_a1:Q_a2])
        ret_re=re_in[findbestPatveFixd(mNo_P12_itp,nom_Q,re_in)]
    return ret_re
def findbestPatveFixd(mNo_P12_itp,nom_Q,re_in):
    cr=np.corrcoef(mNo_P12_itp,nom_Q)
    reix=np.argmax(cr[re_in.size,:-1])
    return reix
def retrieve_pol(P12_in,re_in,ve_in,MieA,Q_in,ScatA,VZA,SZA):
    '''
    Breon et al. (1998) polarimetric retrievals (in my way)
    Exploits only the pattern of the polarize reflection AFTER normalizing to the clodbow value.
        (It seems normalization to the cloudbow value is not important. Only viewing gemetry correction would be enough)
    Uses correlation coefficients.
    02/05/2018:
    Good for smaller library. When have more re and ve values in the library does not work.
    Should invetigate furthur to confirm whether a bug or not.
    -----------------------------------------------------
    P12_in: 2D Array [re_in,ve_in,MieA] P12 from Mie calculations
    re_in: 1D Array
    ve_in: 1D Array
    MieA: 1D array scattering angle
    Q_in: 3D/2D array (single column/2D image) [ScatA,nxdim,nydim] or [ScatA]
    VZA: Viewing Zenith
    SZA: Solar Zenith
    return:
        re_ret: if Q_in[ScatA,nxdim,nydim] -> [nxdim,nydim]
    '''
    muS=np.cos(np.deg2rad(SZA))
    muV=abs(np.cos(np.deg2rad(VZA)))
    p12_a1=np.argwhere(MieA==130)
    p12_a2=np.argwhere(MieA==160)
    Q_a1=np.argmin(abs(ScatA-130))
    Q_a2=np.argmin(abs(ScatA-160))
    if Q_a1>Q_a2:
        Q_a1,Q_a2=Q_a2,Q_a1
    #P12 normalization
    cbow_P12=np.max(-P12_in[:,:,p12_a1:p12_a2],axis=2)#maximum -P12 at the cloud bow
    cbow_P12=cbow_P12*0+1 #to temperolly remove normalization
    mNo_avP12=np.einsum('ijk,ij->ijk',-P12_in,1/cbow_P12)#avP12 normalized by the maximum at the cloudbow
    
    #interpolation
    theP12=MieA[p12_a1-1:p12_a2+2]
    theQ=ScatA[Q_a1:Q_a2]
    mNo_P12=mNo_avP12[:,:,p12_a1-1:p12_a2+2]
    mNo_P12_itp=np.zeros((re_in.size,ve_in.size,theQ.size),dtype='float')
    for i in np.arange(0,re_in.size):
        for j in np.arange(0,ve_in.size):
            itpF=itp.interp1d(theP12,mNo_P12[i,j,:])#interpolated avP12 to match with observed scattering angles
            mNo_P12_itp[i,j,:]=itpF(theQ)
    
    #Q_in normalization
    if Q_in.ndim==3:
        Qdom=np.einsum('ijk,i->ijk',-Q_in,4*(muS+muV))    
        nom_Q=Qdom[Q_a1:Q_a2,:,:]/1#np.max(Qdom[Q_a1:Q_a2,:,:],axis=0) normalization temporally removed
        ret_re=np.zeros((nom_Q.shape[1],nom_Q.shape[2]),dtype=float)
        ret_ve=np.zeros((nom_Q.shape[1],nom_Q.shape[2]),dtype=float)
        for i in np.arange(nom_Q.shape[1]):
            for j in np.arange(nom_Q.shape[2]):
#                print('%0.2f'%(100*float(i*1)*(j*1)/ret_re.size))
                re_ix,ve_ix=findbestPat2(mNo_P12_itp,nom_Q[:,i,j],re_in,ve_in)
                ret_re[i,j],ret_ve[i,j]=re_in[re_ix],ve_in[ve_ix]
    elif Q_in.ndim==1:
        Qdom=-Q_in*4*(muS+muV)
        nom_Q=Qdom[Q_a1:Q_a2]/np.max(Qdom[Q_a1:Q_a2])
        re_ix,ve_ix=findbestPat2(mNo_P12_itp,nom_Q,re_in,ve_in)
        ret_re,ret_ve=re_in[re_ix],ve_in[ve_ix]
    return ret_re,ret_ve
def findbestPat2(mNo_P12_itp,nom_Q,re_in,ve_in):
    cr=np.corrcoef(mNo_P12_itp.reshape(re_in.size*ve_in.size,nom_Q.size),nom_Q)
    ix=np.argmax(cr[re_in.size*ve_in.size,:-1])
    re_ix=ix/re_in.size
    ve_ix=np.remainder(ix,ve_in.size)
    return re_ix,ve_ix


'''
Parametric curvefitting method
==============================
'''
class Pmat(object):
    '''
    To adjust the original Polarized phase matrix (P12) to the observation geometry
    Only useful to do polarimetric retrievals based on parametric curve fitting method.
    2019/02/11:
        Different fitting functions can be used by changing method
    method: 
        Breon   :Breon et. al. 2005 
        BreonMod:Breon et. al. 2005 albeit F=aP+b*cos()**2+c
    Qa1Qa2, SZA:
        Qa1Qa2, SZA = None,None => uses self.findQ12a2_2() to get Q_a1, Q_a2
        Qa1Qa2, SZA = <any>,'160' => uses SZA to get predefined Q_a1, Q_a2
        Qa1Qa2, SZA = [55,85],None => uses given Qa1Qa2    
    '''
    def __init__(self,Re,Ve,MieA,MieP12,ObsA,method='Breon',primaryBow=False,Qa1Qa2=None,SZA=None):
        print('Pol. ret. method: '+method+' primaryBow: '+str(primaryBow))
        preQa1Qa2 = {160:[55,85],140:[35,65],120:[15,45]}
        self.method=method
        if primaryBow:
            slt=135 #left scattering angle
            srt=165 #right scattering angle
        else:
            slt=145 #left scattering angle
            srt=165 #right scattering angle
        self.Re=Re
        self.Ve=Ve
        #Selecting matching observation scattering angles to the library
        p12_a1=int(np.argwhere(MieA==slt))
        p12_a2=int(np.argwhere(MieA==srt))
        if SZA is not None:
            Qa1Qa2 = preQa1Qa2[SZA]
        elif Qa1Qa2 is None:
            Qa1Qa2 = np.array(self.findQa1a2_2(ObsA,slt,srt))
        print('Qa1Qa2: '+str(Qa1Qa2))
        Q_a1,Q_a2=Qa1Qa2[0],Qa1Qa2[1]
        if Q_a1>Q_a2:
            Q_a1,Q_a2=Q_a2,Q_a1
        theQ=ObsA[Q_a1:Q_a2]
        theP12=MieA[p12_a1-1:p12_a2+2]
        MieP12=MieP12[:,:,p12_a1-1:p12_a2+2]
        MieP12_itp=np.zeros((Re.size,Ve.size,theQ.size),dtype='float')
        for i in np.arange(0,Re.size):
            for j in np.arange(0,Ve.size):
                itpF=itp.interp1d(theP12,MieP12[i,j,:])#interpolated avP12 to match with observed scattering angles
                MieP12_itp[i,j,:]=itpF(theQ)
        self.MieP12_itp=MieP12_itp#Interpolated P12 library
        self.ScatA=theQ#Scattering angles of ther self.MieP12_itp
        self.Q_a1=Q_a1
        self.Q_a2=Q_a2
    def findQa1a2_1(self,ObsA,slt,srt):
        Q_a1=np.argmin(abs(ObsA-slt))
        Q_a2=np.argmin(abs(ObsA-srt))
        return Q_a1,Q_a2
    def findQa1a2_2(self,ObsA,slt,srt):
        q_a1=np.argwhere(abs(ObsA-slt)<0.01)
        q_a2=np.argwhere(abs(ObsA-srt)<0.01)
        if q_a1.size==1 and q_a2.size==1:
            print('only 1 occurrence of each angle.')
            Q_a1,Q_a2=self.findQa1a2_1(ObsA,slt,srt)
        elif q_a1.size==1 or q_a2.size==1:
            print('either one angle has only one occurrence')
            argsDiff=abs(q_a1-q_a2)
            if argsDiff.size==q_a1.size:
                Q_a1=q_a1[np.where(argsDiff==(srt-slt))]
                Q_a2=q_a2
            elif argsDiff.size==q_a2.size:
                Q_a2=q_a2[np.where(argsDiff==(srt-slt))]
                Q_a1=q_a1
        elif q_a1.size==q_a2.size:
            print('both angles has same number of occurrences')
            argsDiff=abs(q_a1-q_a2)
            Q_a1=q_a1[np.where(argsDiff==(srt-slt))]
            Q_a2=q_a2[np.where(argsDiff==(srt-slt))]
        else:
            print('!!! Check whether required scattering angle range is available in the ObsA !!!!')
        return np.squeeze(Q_a1),np.squeeze(Q_a2)
    def set_reve(self,re,ve):
        '''
        To set re and ve values. 
        '''
        self.re=re;self.ve=ve
        self.re_i=np.argmin(np.abs(self.re-self.Re))
        self.ve_i=np.argmin(np.abs(self.ve-self.Ve))
    def getP(self,Theta):
        '''
        Theta MUST be same as ScatA
        '''
        return self.MieP12_itp[self.re_i,self.ve_i,:]
    def imitateF(self,Theta,a,b,c):
        '''
        To imitate F() in do_fitting to reproduce fitted curve later
        '''
        if self.method=='Breon':
            y=a*self.getP(Theta)+b*np.deg2rad(Theta)+c
        elif self.method=='BreonMod':
            y=a*self.getP(Theta)+b*np.cos(np.deg2rad(Theta))**2+c
        return y

def do_fitting(x,y,P,ygabc=None):
    '''
    Polarimetric cloud property retrievals from parametric curve fitting
    P:PMat object
    x: Observed scattering angles
    y: Corrected polarized reflectance observations
    method:
        Breon   : Breon et. al. 2005
        BreonMod: Breon et. al. 2005 F=aP()+b*cos()**2+c
    Re,Ve,abc
    Qsqr: Quality index (Breon et. al. 2005)
    '''
    if P.method=='Breon':
        return fit_Breon(x,y,P,ygabc)
    elif P.method=='BreonMod':
        return fit_BreonMod(x,y,P,ygabc)
def fit_Breon(x,y,P,ygabc):
    '''
    F=a*P.getP(Theta)+b*np.deg2rad(Theta)+c
    '''
    def F(Theta,a,b,c):
        return a*P.getP(Theta)+b*np.deg2rad(Theta)+c
    return fit(x,y,P,F,ygabc)
def fit_BreonMod(x,y,P,ygabc):
    '''
    F=a*P.getP(Theta)+b*np.cos(np.deg2rad(Theta))**2+c
    '''
    def F(Theta,a,b,c):
        return a*P.getP(Theta)+b*np.cos(np.deg2rad(Theta))**2+c
    return fit(x,y,P,F,ygabc)    
def fit(x,y,P,F,ygabc):
    R_sq=0.0
    for i in P.Re:
        for j in P.Ve:
            P.set_reve(i,j)
            if ygabc is None:
                popt, pcov = curve_fit(F, x, y)
            else:
                popt, pcov = curve_fit(F, x, y,ygabc)
            #Finding R^2 value
            residuals = y- F(x, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y-np.mean(y))**2)
            r_sq = 1 - (ss_res / ss_tot)
            fl=open('Rsqs.dat',"a")
            fl.write("\n%0.6f\t%0.2f\t%0.2f"%(r_sq,i,j))
            fl.close()

            if r_sq>R_sq:
                R_sq=r_sq
                abc=popt
                Re=i
                Ve=j
                fl=open('Rsqs_max.dat',"a")
                fl.write("\n%0.6f\t%0.2f\t%0.2f"%(R_sq,Re,Ve))
                fl.close()
    Qsqr=abc[0]**2*(np.mean(P.getP(x)**2)-np.mean(P.getP(x))**2)/np.mean((y-F(x,*abc))**2)
    
    return Re,Ve,abc,Qsqr,R_sq


def fitBreon(x,y,P,ygabc=None):
    def F(Theta,a,b,c):
        return a*P.getP(Theta)+b*np.deg2rad(Theta)+c
    def fit(x,y,P,F,ygabc):
        R_sq=0.0
        for i in P.Re:
            for j in P.Ve:
                P.set_reve(i,j)
                popt, pcov = curve_fit(F, x, y,ygabc)
                #Finding R^2 value
                residuals = y- F(x, *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_sq = 1 - (ss_res / ss_tot)
                fl=open('Rsqs.dat',"a")
                fl.write("\n%0.6f\t%0.2f\t%0.2f"%(r_sq,i,j))
                fl.close()
    
                if r_sq>R_sq:
                    R_sq=r_sq
                    abc=popt
                    Re=i
                    Ve=j
                    fl=open('Rsqs_max.dat',"a")
                    fl.write("\n%0.6f\t%0.2f\t%0.2f"%(R_sq,Re,Ve))
                    fl.close()
        Qsqr=abc[0]**2*(np.mean(P.getP(x)**2)-np.mean(P.getP(x))**2)/np.mean((y-F(x,*abc))**2)
        
        return Re,Ve,abc,Qsqr,R_sq
    return fit(x,y,P,F,ygabc)    

def fitBreon_noRsq(x,y,P,ygabc=None):
    '''
    Flags:
        0: Successful retrievals
        1: Optimal parameters not found (curve_fit). maxfev=800 reached"
            o Check how many these flags were given in each case. Seems not a consistent issue.
    '''
    def F(Theta,a,b,c):
        return a*P.getP(Theta)+b*np.deg2rad(Theta)+c
    def fit(x,y,P,F,ygabc):
        R_sq=0.0
        flag=0.0
        abc=np.array([1,1,1])*np.nan;Re,Ve=np.nan,np.nan #For failed retrievals
        for i in P.Re:
            for j in P.Ve:
                P.set_reve(i,j)
                try:
                    popt, pcov = curve_fit(F, x, y,p0=ygabc)
                except RuntimeError:
                    print("Optimal parameters not found: Number of calls to function has reached maxfev = 800")
                    flag=1.0
                #Finding R^2 value
                residuals = y- F(x, *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_sq = 1 - (ss_res / ss_tot)
#                fl=open('Rsqs.dat',"a")
#                fl.write("\n%0.6f\t%0.2f\t%0.2f"%(r_sq,i,j))
#                fl.close()
    
                if r_sq>R_sq:
                    R_sq=r_sq
                    abc=popt
                    Re=i
                    Ve=j
#                    fl=open('Rsqs_max.dat',"a")
#                    fl.write("\n%0.6f\t%0.2f\t%0.2f"%(R_sq,Re,Ve))
#                    fl.close()
        Qsqr=abc[0]**2*(np.mean(P.getP(x)**2)-np.mean(P.getP(x))**2)/np.mean((y-F(x,*abc))**2)
        
        return Re,Ve,abc,Qsqr,R_sq,flag
    return fit(x,y,P,F,ygabc)    
def fitBreonMod(x,y,P,ygabc=None):
    '''
    F=a*P.getP(Theta)+b*np.cos(np.deg2rad(Theta))**2+c
    Save *.dat files that contain Rsq values.
    '''
    def F(Theta,a,b,c):
        return a*P.getP(Theta)+b*np.cos(np.deg2rad(Theta))**2+c
    def fit(x,y,P,F,ygabc):
        R_sq=0.0
        for i in P.Re:
            for j in P.Ve:
                P.set_reve(i,j)
                if ygabc is None:
                    popt, pcov = curve_fit(F, x, y)
                else:
                    popt, pcov = curve_fit(F, x, y,ygabc)
                #Finding R^2 value
                residuals = y- F(x, *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_sq = 1 - (ss_res / ss_tot)
                fl=open('Rsqs.dat',"a")
                fl.write("\n%0.6f\t%0.2f\t%0.2f"%(r_sq,i,j))
                fl.close()
    
                if r_sq>R_sq:
                    R_sq=r_sq
                    abc=popt
                    Re=i
                    Ve=j
                    fl=open('Rsqs_max.dat',"a")
                    fl.write("\n%0.6f\t%0.2f\t%0.2f"%(R_sq,Re,Ve))
                    fl.close()
        Qsqr=abc[0]**2*(np.mean(P.getP(x)**2)-np.mean(P.getP(x))**2)/np.mean((y-F(x,*abc))**2)
        
        return Re,Ve,abc,Qsqr,R_sq
    return fit(x,y,P,F,ygabc)    
def fitBreonMod_noRsq(x,y,P,ygabc=None):
    '''
    Same as fitBreonMod but does not save Rsq values.
    '''
    def F(Theta,a,b,c):
        return a*P.getP(Theta)+b*np.cos(np.deg2rad(Theta))**2+c
    def fit(x,y,P,F,ygabc):
        R_sq=0.0
        for i in P.Re:
            for j in P.Ve:
                P.set_reve(i,j)
                if ygabc is None:
                    popt, pcov = curve_fit(F, x, y)
                else:
                    popt, pcov = curve_fit(F, x, y,ygabc)
                #Finding R^2 value
                residuals = y- F(x, *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_sq = 1 - (ss_res / ss_tot)
#                fl=open('Rsqs.dat',"a")
#                fl.write("\n%0.6f\t%0.2f\t%0.2f"%(r_sq,i,j))
#                fl.close()
    
                if r_sq>R_sq:
                    R_sq=r_sq
                    abc=popt
                    Re=i
                    Ve=j
#                    fl=open('Rsqs_max.dat',"a")
#                    fl.write("\n%0.6f\t%0.2f\t%0.2f"%(R_sq,Re,Ve))
#                    fl.close()
        Qsqr=abc[0]**2*(np.mean(P.getP(x)**2)-np.mean(P.getP(x))**2)/np.mean((y-F(x,*abc))**2)
        
        return Re,Ve,abc,Qsqr,R_sq
    return fit(x,y,P,F,ygabc)    
'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
----------- VERTICALLY WEIGHTED PSEUDO RETRIEVALS------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
class LES_pseudo_rets(object):
    '''
    This object will be used in DHARMA_onmp.setup_reVW to compute and setup vertically weighted psudo retrievals.
    Use DHARMA_onmp.find_reVW if you don't want save or reload results.
    -----------------------------------------------------------------------------
    les_onmp: DHARMA_onmp object
    fpath: path to search previously generated psudo retrieval results.
    '''
    def __init__(self,les_onmp,mie_name,mie_path,dgSZA,dgVZA,a,b,band=None,fpath=None,replace=None):
        self.replace=replace# '0' to load from existing files. '1' to replace previous results.
        les_onmp.readMie(mie_name,mie_path)
        if band is None:
            band=1
            print('Deafault wavelength %0.3f being selected'%(les_onmp.c_mie.wvl[band]))
            print('Available wavelengths:'+str(les_onmp.c_mie.wvl))        
        fname=les_onmp.fname.split('.',1)[0]+'_'+mie_name+'_SZA%03d_VZA%03d_b'%(dgSZA,dgVZA)+\
            ("%0.3f"%les_onmp.c_mie.wvl[band]+'a%0.1f_b%0.1f'%(a,b)).replace('.','p')
        self.fname=fname.replace('.','p')
        self.mie_path=mie_path
        self.les_onmp=les_onmp
        if fpath==None:
            self.fpath='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_simulations/'
        else:
            self.fpath=fpath
        print('File name: '+self.fname)
        test=os.path.isfile(self.fpath+self.fname+'.hdf5')
        if test:
            print('File already exist in: '+self.fpath)
            if self.replace is None:
                self.replace=input('Replace/Reload 1/0:')
        if self.replace=='1' and test:
            print('Previous results will be replaced!')
        elif self.replace=='0' and test:
            print('Previous results will be loaded!')
            self.loadVW()
        elif self.replace=='0' and not(test):
            self.replace='1'
        elif replace is None:
            self.replace='1'
    def saveVW(self,):
        out_name=self.fpath+self.fname+'.hdf5'
        f=h5py.File(out_name,'w')
        f.attrs['DHARMA_LES_PATH']=self.fpath
        f.attrs['MIE_PATH']=self.mie_path
        
        PC=f.create_dataset('Re',data=self.Re)
        PC.dims[0].label='x'
        PC.dims[1].label='y'
        PC.attrs['units']='Microns'
        PC.attrs["long_name"]='Vertically_Weighted_Cloud_Effective_Radius'

        PC=f.create_dataset('Ve',data=self.Ve)
        PC.dims[0].label='x'
        PC.dims[1].label='y'
        PC.attrs['units']='None'
        PC.attrs["long_name"]='Vertically_Weighted_Cloud_Effective_Variance'
        
        PC=f.create_dataset('Re_dN_vw',data=self.Re_dN_vw)
        PC.dims[0].label='x'
        PC.dims[1].label='y'
        PC.attrs['units']='Microns'
        PC.attrs["long_name"]='Cloud_Effective_Radius_from_Vertically_weighted_DSD'
    
        PC=f.create_dataset('Ve_dN_vw',data=self.Ve_dN_vw)
        PC.dims[0].label='x'
        PC.dims[1].label='y'
        PC.attrs['units']='None'
        PC.attrs["long_name"]='Cloud_Effective_Variance_from_Vertically_weighted_DSD'

        PC=f.create_dataset('Tau',data=self.Tau)
        PC.dims[0].label='x'
        PC.dims[1].label='y'
        PC.attrs['units']='None'
        PC.attrs["long_name"]='Vertically_Weighted_Cloud_Optical_Thickness'

        PC=f.create_dataset('dN',data=self.dN)
        PC.dims[0].label='x'
        PC.dims[1].label='y'
        PC.attrs['units']='counts/cm^3'
        PC.attrs["long_name"]='Vertically_Weighted_Number_Concentration'
        
        PC=f.create_dataset('Re_3D',data=self.re_tau)
        PC.dims[0].label='tau'
        PC.dims[1].label='x'
        PC.dims[2].label='y'
        PC.attrs['units']='Microns'
        PC.attrs["long_name"]='Cloud_Effective_Radius_3D'
    
        PC=f.create_dataset('Ve_3D',data=self.ve_tau)
        PC.dims[0].label='tau'
        PC.dims[1].label='x'
        PC.dims[2].label='y'
        PC.attrs['units']='None'
        PC.attrs["long_name"]='Cloud_Effective_Variance_3D'
    
        PC=f.create_dataset('W',data=self.w_tau)
        PC.dims[0].label='tau'
        PC.dims[1].label='x'
        PC.dims[2].label='y'
        PC.attrs['units']='None'
        PC.attrs["long_name"]='Vertical_weighting_function'
        
        PC=f.create_dataset('tau',data=self.tau)
        PC.dims[0].label='z'
        PC.dims[1].label='x'
        PC.dims[2].label='y'
        PC.attrs['units']='None'
        PC.attrs['long_name']='Cloud_Optical_Thickness_3D'    
    
        PC=f.create_dataset('x',data=self.les_onmp.x/1e3)
        PC.attrs['units']='km'
        PC.attrs['long_name']='X distance'
    
        PC=f.create_dataset('y',data=self.les_onmp.y/1e3)
        PC.attrs['units']='km'
        PC.attrs['long_name']='Y distance'
        
        PC=f.create_dataset('z',data=self.les_onmp.DHARMA.z/1e3)
        PC.attrs['units']='km'
        PC.attrs['long_name']='Z distance'
        f.close()
        print(out_name+' SAVED!')

    def loadVW(self,):
        filename=self.fpath+self.fname+'.hdf5'
        f=h5py.File(filename,'r')
        self.fpath    = f.attrs['DHARMA_LES_PATH']
        self.mie_path = f.attrs['MIE_PATH']
        self.Re       = f['Re'][:] 
        self.Ve       = f['Ve'][:]    
        self.Re_dN_vw = f['Re_dN_vw'][:]
        self.Ve_dN_vw = f['Ve_dN_vw'][:]
        self.Tau      = f['Tau'][:]    
        self.dN       = f['dN'][:]        
        self.re_tau   = f['Re_3D'][:]    
        self.ve_tau   = f['Ve_3D'][:]    
        self.w_tau    = f['W'][:]        
        self.tau      = f['tau'][:]    
        self.x        = f['x'][:]   
        self.y        = f['y'][:]        
        self.z        = f['z'][:]
        f.close()
