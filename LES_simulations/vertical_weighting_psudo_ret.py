#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
********************************************
Created on Wed Feb 27 09:35:39 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************

Implementation of the bidirectional weghting retrievals on LES simulations
Two variable parametric vertical weighting (check my 07/30/2018 notes)
Platnick(2000), Alexandrove (2012), Zhang(2017), Miller(2017)
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from textwrap import wrap
from matplotlib.animation import FuncAnimation
import os, h5py
  
from cpnLES_MSCARTlib import DHARMA_onmp,LES_field
from cpnMielib import MieSet
import cpnCommonlib as cpn

def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'figures/')

class LES_psudo_rets(object):
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
            ("%0.3f"%les_onmp.c_mie.wvl[band]).replace('.','p')+'a%d_b%d'%(a,b)
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
        self.fpath=f.attrs['DHARMA_LES_PATH']
        self.mie_path=f.attrs['MIE_PATH']
        self.Re=f['Re'][:]   
        self.Ve=f['Ve'][:]    
        self.Tau=f['Tau'][:]    
        self.dN=f['dN'][:]        
        self.re_tau=f['Re_3D'][:]    
        self.ve_tau=f['Ve_3D'][:]    
        self.w_tau=f['W'][:]        
        self.tau=f['tau'][:]    
        self.x=f['x'][:]   
        self.y=f['y'][:]        
        self.z=f['z'][:]
        f.close()
        
if __name__=='__main__':
    cpn.setup_figures(plt)
    #DYCOMS2------------------------------------------------------------------------------------------
    DYCOMS2_path  ='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/DYCOMS/Ackerman_DYCOMS2_25bins/'
    DYCOMS2=DHARMA_onmp('dharma_008036.cdf',DYCOMS2_path );DYCOMS2.readDHARMA()
    mie_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/'
    mie_name='DYCOM2_dharma_008036_mie_470_860_2p13'
    lesCname='DYCOMS2';lesBinEd=None
    a=1;b=0;SZA=40;VZA=0
    band=1
    wvl={1:'0p860',2:'2p13'}
#    Re_vw,Ve_vw,dN_vw,re_tau,ve_tau,w_tau,tau=DYCOMS2.find_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)        
    DYCOMS2.setup_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)
    #DYCOMS2 end ---------------------------------------------------------------------------------------
    #ATEXc ----------------------------------------------------------------------------------------------
    ATEXc_path  ='/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_ATEX_highres/CCN40/'  
    ATEXc=DHARMA_onmp('dharma_007877.cdf',ATEXc_path );ATEXc.readDHARMA()
    mie_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/'
    mie_name='DYCOM2_dharma_008036_mie_470_860_2p13'
    lesCname='ATEXc'
    a=1;b=0;SZA=40;VZA=0
    band=1
    wvl={1:'0p860',2:'2p13'}
#    Re_vw,Ve_vw,dN_vw,re_tau,ve_tau,w_tau,tau=ATEXc.find_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)        
    ATEXc.setup_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)
    #ATEXc end---------------------------------------------------------------------------------------------
    #ATEXp -------------------------------------------------------------------------------------------------------------------
    ATEXp_path  ='/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_ATEX_highres/CCN600/'  
    ATEXp=DHARMA_onmp('dharma_013067.cdf',ATEXp_path );ATEXp.readDHARMA()
    mie_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/'
    mie_name='DYCOM2_dharma_008036_mie_470_860_2p13'
    lesCname='ATEXp'
    a=1;b=0;SZA=40;VZA=0;band=1;wvl={1:'0p860',2:'2p13'}
#    Re_vw,Ve_vw,dN_vw,re_tau,ve_tau,w_tau,tau=ATEXp.find_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)        
    ATEXp.setup_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)
    #ATEXp end ------------------------------------------------------------------------------------------------------------------------
    #RICO ------------------------------------------------------------------------------------------------------------------------------
    RICO_path='/umbc/xfs1/zzbatmos/users/mild1/Research/LES/LES_Satellite_Simulator_DJM/LES_Simulations/Ackerman_RICO_NEW/'
    RICO=DHARMA_onmp('dharma_005044.cdf',RICO_path );RICO.readDHARMA()
    mie_path='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_simulations/'
    mie_name='RICO_dharma_005044_mie_470_860_2p13'
    lesCname='RICO'
    a=1;b=0;SZA=40;VZA=0;band=1;wvl={1:'0p860',2:'2p13'}
#    Re_vw,Ve_vw,dN_vw,re_tau,ve_tau,w_tau,tau=RICO.find_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)        
    RICO.setup_reVW(mie_name,lesCname,SZA,VZA,a=a,b=b,mie_path=mie_path,band=band)
    #RICO end---------------------------------------------------------------------------------------------------------------------------
    case=DYCOMS2
    #Tau, Re_vw and Ve_vw plot
    fig2,ax2=plt.subplots(1,3,figsize=(14,5),subplot_kw={'aspect':'equal'})
    fig2_ttl=case.fname.split('.',1)[0]+'_'+mie_name+'_a%db%dSZA%dVZA%dwvl'%(a,b,SZA,VZA)+wvl[band]+'tau_CER_CEV_vw'
    fig2.suptitle("\n".join(wrap(fig2_ttl,80)))
    ctf3=ax2[0].contourf(case.x/1e3,case.y/1e3,case.VW.tau[0,:],np.linspace(0,40,50),cmap=plt.cm.jet)
    ctf1=ax2[1].contourf(case.x/1e3,case.y/1e3,case.VW.Re,np.linspace(0,25,100),cmap=plt.cm.jet)   
    ctf2=ax2[2].contourf(case.x/1e3,case.y/1e3,case.VW.Ve,np.linspace(0,0.3,50),cmap=plt.cm.jet)
    cpn.add_cb(fig2,ctf3,ax2[0],ticks=np.arange(0,40.1,10),orientation='vertical',label=r'$COT$')
    cpn.add_cb(fig2,ctf1,ax2[1],ticks=np.arange(0,25.1,5),orientation='vertical',label=r'$CER_{VW}$ ($\mu m$)')
    cpn.add_cb(fig2,ctf2,ax2[2],ticks=np.arange(0,0.31,0.1),orientation='vertical',label=r'$CEV_{VW}$')
    for n, x in enumerate(ax2.flat):
        x.set_xlabel('X (km)')
        x.set_ylabel('Y (km)')
    fig2.tight_layout(rect=[0,0,1,0.99])
    fig2.show()
    
    #COT
    array = tau.swapaxes(0,2)
    x,y,z = case.x/1e3,case.y/1e3,case.DHARMA.z/1e3
    xslice_i=64
    fig3 = plt.figure(figsize=(16,6))
    fig3_ttl = case.fname.split('.',1)[0]+mie_name+'_COT'
    ax3p1 = fig3.add_subplot(131, aspect='equal')
    ax3p2 = fig3.add_subplot(132, projection='3d')
    ax3p3 = fig3.add_subplot(133, projection='3d')
    mnmx={'vmin':0,'vmax':40}
    cpn.array3D_slice('x',array,x,y,z,64,fig3,ax3p2,mnmx,cb=False)
    cpn.array3D_slice('y',array,x,y,z,64,fig3,ax3p3,mnmx,cb=False)
    ctf=ax3p1.contourf(case.x*1e-3,case.y*1e-3,tau[0,:],np.linspace(0,40,50),cmap=plt.cm.jet)
    cpn.add_cb(fig3,ctf,ax3p1,ticks=np.arange(0,40.1,10),orientation='vertical',label='COT')
    ax3p1.set_xlabel('X (km)')
    ax3p1.set_ylabel('Y (km)')
    fig3.suptitle(fig3_ttl)
    fig3.tight_layout(rect=[0,0,1,0.98])
    fig3.show()

    #re_tau
    array = re_tau.swapaxes(0,2)
    x,y,z = case.x/1e3,case.y/1e3,case.DHARMA.z/1e3
    xslice_i=64
    fig4 = plt.figure(figsize=(16,6))
    fig4_ttl = case.fname.split('.',1)[0]+mie_name+'_re_tau'
    ax4p1 = fig4.add_subplot(131, projection='3d')
    ax4p2 = fig4.add_subplot(132, projection='3d')
    ax4p3 = fig4.add_subplot(133, projection='3d')
    mnmx={'vmin':0,'vmax':50}
    cpn.array3D_slice('z',array,x,y,z,90,fig4,ax4p1,mnmx={'vmin':0,'vmax':10},cb_label=r'$r_e(\tau)$ $(\mu m)$')
    cpn.array3D_slice('x',array,x,y,z,64,fig4,ax4p2,mnmx=mnmx,cb_label=r'$r_e(\tau)$ $(\mu m)$')
    cpn.array3D_slice('y',array,x,y,z,64,fig4,ax4p3,mnmx=mnmx,cb_label=r'$r_e(\tau)$ $(\mu m)$')
    fig4.suptitle(fig4_ttl)
#    fig4.tight_layout(rect=[0,0,1,0.98])
    fig4.show()

    #re_tau and re_vw and ve_vw gif
    fig5,ax5=plt.subplots(1,2,subplot_kw={'aspect':'equal'})
    fig5_ttl=case.fname.split('.',1)[0]+mie_name+'re_tau_re_vw'    
    xcens=case.x*1e-3#km
    ycens=case.x*1e-3#km
    zi=0
    v=np.linspace(0,25,50)
    Re_vw_tau=re_tau*w_tau
    fig5.suptitle('Z = %0.2f km'%(case.DHARMA.z[zi]/1e3))
    ctf1=[ax5[0].contourf(xcens,ycens,re_tau[zi,:],np.linspace(0,250,50),cmap=plt.cm.jet)]
    ctf2=[ax5[1].contourf(xcens,ycens,np.trapz(Re_vw_tau[:zi,:],tau[:zi,:],axis=0),v,cmap=plt.cm.jet)]
    cpn.add_cb(fig5,ctf1[0],ax5[0],pad=0.3,ticks=np.arange(0,250.1,50),label=r'$r_e(Z)$')
    cpn.add_cb(fig5,ctf2[0],ax5[1],pad=0.3,ticks=np.arange(0,25.1,5),label=r'$r_e^{vw}(Z)$')
    fig5.tight_layout()
    fig5.show()
    def update(zi):
        fig5.suptitle('Z = %0.2f km'%(case.DHARMA.z[zi]/1e3))
        for coll in ctf1[0].collections: 
            coll.remove() 
        ctf1[0]=ax5[0].contourf(xcens,ycens,re_tau[zi,:],np.linspace(0,250,50),cmap=plt.cm.jet)    
        ctf2[0]=ax5[1].contourf(xcens,ycens,np.trapz(Re_vw_tau[:zi,:],tau[:zi,:],axis=0),v,cmap=plt.cm.jet)
        plt.draw()
        return
    fig5.show()
    anim = FuncAnimation(fig5, update, frames=np.arange(0, 96), interval=100)
    anim.save(fig5_ttl+'.gif', dpi=200, writer='imagemagick')
    
    #re_tau,ve_tau and re_vw,ve_vw with altitude gif
    fig6=plt.figure(figsize=(7,6))
    fig6_ttl=case.VW.fname.split('.',1)[0]+'LES_re_ve_summary'  
    fig6.suptitle("\n".join(wrap(fig6_ttl,60)),size=12)
    ax6p1=plt.subplot2grid((2,7),(0,0),rowspan=2,)
    ax6p2=plt.subplot2grid((2,7),(0,1),colspan=3,aspect='equal')
    ax6p3=plt.subplot2grid((2,7),(0,4),colspan=3,aspect='equal')
    ax6p4=plt.subplot2grid((2,7),(1,1),colspan=3,aspect='equal')
    ax6p5=plt.subplot2grid((2,7),(1,4),colspan=3,aspect='equal')
    xcens=case.x*1e-3#km
    ycens=case.x*1e-3#km
    ax6p1.plot(case.VW.tau.sum(axis=(1,2)),case.DHARMA.z/1e3)
    ax6p1.set_xscale('log')
    zi=95
    Re_vw_tau=case.VW.re_tau*case.VW.w_tau
    Ve_vw_tau=case.VW.ve_tau*case.VW.w_tau
    line,=ax6p1.plot(np.array([0,1]),case.DHARMA.z[zi]/1e3*np.array([1,1]),'k-',linewidth=3)
    ax6p1.set_ylim(0,case.DHARMA.z[-1]/1e3)
    ctf1=[ax6p2.contourf(xcens,ycens,case.VW.re_tau[zi,:],np.linspace(0,250,50),cmap=plt.cm.jet)]
    ctf2=[ax6p3.contourf(xcens,ycens,np.trapz(Re_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),np.linspace(0,25,50),cmap=plt.cm.jet)]
    ctf3=[ax6p4.contourf(xcens,ycens,case.VW.ve_tau[zi,:],np.linspace(0,0.5,50),cmap=plt.cm.jet)]
    ctf4=[ax6p5.contourf(xcens,ycens,np.trapz(Ve_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),np.linspace(0,0.2,50),cmap=plt.cm.jet)]
    cpn.add_cb(fig6,ctf1[0],ax6p2,pad=0.3,ticks=np.arange(0,250.1,125),label=r'$r_e(Z)$ $\mu m$')
    cpn.add_cb(fig6,ctf2[0],ax6p3,pad=0.3,ticks=np.arange(0,25.1,5),label=r'$r_e^{vw}(Z)$ $\mu m$')
    cpn.add_cb(fig6,ctf3[0],ax6p4,pad=0.3,ticks=np.arange(0,0.51,0.25),label=r'$v_e(Z)$')
    cpn.add_cb(fig6,ctf4[0],ax6p5,pad=0.3,ticks=np.arange(0,0.21,0.1),label=r'$v_e^{vw}(Z)$')
    ax6p1.set_ylabel('Altitude (km)')
    ax6p1.set_xticklabels([0,np.log(case.VW.tau.sum(axis=(1,2)).max())])
    fig6.tight_layout(rect=[0,0,1,0.96])
    fig6.show()
    def updatefig6(zi):
        line.set_ydata(case.DHARMA.z[zi]/1e3*np.array([1,1]))
        ctf1[0]=ax6p2.contourf(xcens,ycens,case.VW.re_tau[zi,:],np.linspace(0,250,50),cmap=plt.cm.jet)
        ctf2[0]=ax6p3.contourf(xcens,ycens,np.trapz(Re_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),np.linspace(0,25,50),cmap=plt.cm.jet)
        ctf3[0]=ax6p4.contourf(xcens,ycens,case.VW.ve_tau[zi,:],np.linspace(0,0.5,50),cmap=plt.cm.jet)
        ctf4[0]=ax6p5.contourf(xcens,ycens,np.trapz(Ve_vw_tau[:zi,:],case.VW.tau[:zi,:],axis=0),np.linspace(0,0.2,50),cmap=plt.cm.jet)
        plt.draw()
    anim = FuncAnimation(fig6, updatefig6, frames=np.arange(0, 96), interval=100)
    anim.save('figures/'+fig6_ttl+'dpi200.gif', dpi=200, writer='imagemagick')
    
    
    #plotting optical thickness
    xcens=case.x*1e-3#km
    ycens=case.x*1e-3#km
#    tau3=np.einsum('izxy->ixy',case.dTau)
    tau3=case.zTau[:,-1,:]
    fig1,ax1=plt.subplots(subplot_kw={'aspect':'equal'})
    fig1_ttl=case.fname.split('.',1)[0]+mie_name+'_COT'
    fig1.suptitle(fig1_ttl)
    ctf=ax1.contourf(xcens,ycens,tau3[1],np.linspace(0,40,50),cmap=plt.cm.jet)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    fig1.colorbar(ctf,cax=cax,ticks=np.arange(0,40.1,10),label=r'$\tau$')
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    fig1.show()
    
    #bidirctional weighting
    MSCART_field=LES_field('case_dharma_008036_b0p860.nc',dpath='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/')
    MSCART_field.readLES_field()
    fig2,fig2_ttl=MSCART_field.plot_Tau(v=np.linspace(0,40,50))
    
