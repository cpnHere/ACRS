#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Wed Jun  6 10:44:59 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Polarimetric cloud retrievals (mine)
Priliminary tests.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as itp
from textwrap import wrap
import time, os, sys

from cpnMielib import MieSet, bulk_Mie, PSDs
from cpnLES_MSCARTlib import LES_case, DHARMA_onmp, POLCARTdset
import cpnCommonlib as cpn
from cpnRetrievalslib import  Pol_ret_P12Lib,Pmat
def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'figures/')
def ret_re1(P12_in,re_in,MieA,Q_in,ScatA,VZA,SZA):
    '''
    Exploits only the pattern of the polarize reflection.
    Uses correlation coefficients.
    Seems not that good!
    P12_in: 2D Array[re_in,MieA],MieA,Q_in,ScatA,VZA,SZA
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
    #P12 normalization
    cbow_P12=np.max(-P12_in[:,p12_a1:p12_a2],axis=1)#maximum -P12 at the cloud bow
    mNo_avP12=np.einsum('ij,i->ij',-P12_in,1/cbow_P12)#avP12 normalized by the maximum at the cloudbow
    
    #interpolation
    theP12=MieA[p12_a1-1:p12_a2+2]
    theQ=ScatA[Q_a2:Q_a1]
    mNo_P12=mNo_avP12[:,p12_a1-1:p12_a2+2]
    mNo_P12_itp=np.zeros((re_in.size,theQ.size),dtype='float')
    for i in np.arange(0,re_in.size):
        itpF=itp.interp1d(theP12,mNo_P12[i,:])#interpolated avP12 to match with observed scattering angles
        mNo_P12_itp[i,:]=itpF(theQ)
    
    #Q_in normalization
    if Q_in.ndim==3:
        Qdom=np.einsum('ijk,i->ijk',-Q_in,4*(muS+muV))    
        nom_Q=Qdom[Q_a2:Q_a1,:,:]/np.max(Qdom[Q_a2:Q_a1,:,:],axis=0)
        ret_re=np.zeros((nom_Q.shape[1],nom_Q.shape[2]),dtype=float)
        for i in np.arange(nom_Q.shape[1]):
            for j in np.arange(nom_Q.shape[2]):
                ret_re[i,j]=re[findbestPat(mNo_P12_itp,nom_Q[:,i,j])]
    elif Q_in.ndim==2:
        ret_re=findbestPat(mNo_P12_itp,nom_Q)
    return ret_re
def findbestPat(mNo_P12_itp,nom_Q):
    cr=np.corrcoef(mNo_P12_itp,nom_Q)
    reix=np.argmax(cr[re.size,:-1])
    return reix
def demonstrating_pol_ret():
    #runing Mie code to generate LUTs
    cloud_mie=MieSet('Cloud_mie_470_860_2p13')

#    refs={'lam':[0.443,0.470,0.532,0.550,0.660,0.860,0.865,0.960,1.064,1.24,1.63,2.10,2.13,3.75],\
#          'nr':[1.344590,1.341951,1.337115,1.335943,1.330267,1.324481,1.324372,1.322450,1.320509,1.317240,1.308836,1.291837,1.290108,1.351865],\
#          'nc':[0.8905030E-09,0.7285308E-09,0.1818175E-08,0.2461249E-08,0.1920720E-07,0.3380264E-06,0.3547385E-06,\
#        0.3367559E-05,0.1279514E-05,0.1134782E-04,0.8084244E-04,0.4615043E-03,0.3941974E-03,0.3402400E-02]}
#    lam=[refs['lam'][1],refs['lam'][5],refs['lam'][12]]
#    nr =[refs[ 'nr'][1],refs[ 'nr'][5],refs[ 'nr'][12]]
#    nc =[refs[ 'nc'][1],refs[ 'nc'][5],refs[ 'nc'][12]]
#
#    D=np.linspace(0.1,60,1000)#diameter um
##    ---------------------------------------
#    na=1801
#    ang=np.linspace(0,180,na)
#    ang=np.hstack((np.array([na]),ang))
#    np.savetxt('ang.dat',ang,delimiter='\n',fmt='%.2f')
#    print('Running Mie code')
#    cloud_mie.runMie(lam,nr,nc,D)

    
    cloud_mie.readMie()
    psd=PSDs('mod_gamma',np.asarray(cloud_mie.d,dtype=float))
    #--------------------------------------------------------------------------
    #P12 LUT generation
#    re=np.array([9,10,12,15,20,30])
    re=np.arange(10,20,1)
    ve=0.02
    avP12_860=np.zeros((re.size,cloud_mie.ang.size),dtype=float)
    for i in np.arange(0,re.size):
        psd.set_mod_gamma_norm(re[i],ve)
        bulk=bulk_Mie('bulk_Mie_470_860_2p13',psd=psd,Mie=cloud_mie)
        bulk.cal_bulk_Pmat()
        avP12_860[i,:]=bulk.bulk_P12[:,1]
    p12_a1=np.argwhere(bulk.Mie.ang==130)
    p12_a2=np.argwhere(bulk.Mie.ang==160)

    # Q observations from RT simulations
#    ixx=find_peaks_cwt(data,np.arange(1,50))
    muS=np.cos(np.deg2rad(180.0-dycyx0p860_140.RT.SZA))
    muV=abs(np.cos(np.deg2rad(dycyx0p860_140.RT.VZA)))

    Q_a1=np.argmin(abs(dycyx0p860_140.RT.ScatA-130))
    Q_a2=np.argmin(abs(dycyx0p860_140.RT.ScatA-160))
    Qdom=np.average(np.einsum('ijk,i->ijk',dycyx0p860_140.RT.MeanPRad[:,:,:,1],4*(muS+muV)).reshape(122,128*128),axis=1)#domain averaged Q    
    
    fig1,ax1=plt.subplots()
    fig1_ttl=dycyx0p860_140.RT.fname.split('.',1)[0]+'_Q_obs_P12'
    fig1.suptitle("\n".join(wrap(fig1_ttl,60)))
    ax1.plot(dycyx0p860_140.RT.ScatA[Q_a2:Q_a1],-Qdom[Q_a2:Q_a1],'k--',label=r'$-4Q(\mu_s+\mu_v)(\Theta)$')
    for i in np.arange(0,re.size):
        ax1.plot(bulk.Mie.ang[p12_a1:p12_a2],-avP12_860[i,p12_a1:p12_a2],label=r'$-P_{12}(\Theta);r_e=$%.1f'%re[i])
    ax1.legend(loc='best',prop={'size':10})
    ax1.set_xlabel(r'Scattering Angle ($^{o}$)')
    fig1.show()
    #--------------------------------------------------------------------------
    #-P12 normalization w.r.t. cloudbow maximum
    lang=np.argwhere(bulk.Mie.ang==135)
    rang=np.argwhere(bulk.Mie.ang==145)
    cbow_P12=np.max(-avP12_860[:,lang:rang],axis=1)#maximum -P12 at the cloud bow
    mNo_avP12=np.einsum('ij,i->ij',-avP12_860,1/cbow_P12)#avP12 normalized by the maximum at the cloudbow
    
    fig2,ax2=plt.subplots(1,2,figsize=(10,5))
    fig2_ttl='-P12norm_of_clouds_ve%.2f'%(ve)
    fig2.suptitle(fig2_ttl)
    for i in np.arange(0,re.size):
        ax2[0].plot(bulk.Mie.ang[p12_a1:p12_a2],-avP12_860[i,p12_a1:p12_a2],label=r'$r_e=$%.1f'%re[i])
        ax2[1].plot(bulk.Mie.ang[p12_a1:p12_a2],mNo_avP12[i,p12_a1:p12_a2],label=r'$r_e=$%.1f'%re[i])
    ax2[0].set_title(r'$-P_{12}$',size=10)
    ax2[0].legend(loc='best',prop={'size':10})
    ax2[0].set_xlabel(r'Scattering Angle ($^{o}$)')
    ax2[1].set_title(r'$-P_{12}(\Theta)/-P_{12}(\Theta_{cb})$',size=10)
    ax2[1].set_xlabel(r'Scattering Angle ($^{o}$)')
    fig2.tight_layout()
    fig2.show()

    #--------------------------------------------------------------------------
    #Q domain average normalized w.r.t. the cloudbow value
    lang=np.argwhere(dycyx0p860_140.RT.ScatA==135)
    rang=np.argwhere(dycyx0p860_140.RT.ScatA==145)
    cbow_Qdom=np.max(-Qdom[Q_a2:Q_a1])#maximum -Q at the cloud bow
    mNo_Qdom=-Qdom/cbow_Qdom#Qdom normalized by the maximum at the cloudbow
       

    fig3,ax3=plt.subplots()
    fig3_ttl=dycyx0p860_140.RT.fname.split('.',1)[0]+'_Qweighted_observations'
    fig3.suptitle("\n".join(wrap(fig3_ttl,60)))
    l1,=ax3.plot(dycyx0p860_140.RT.ScatA[Q_a2:Q_a1],-Qdom[Q_a2:Q_a1])
    ax3_2=ax3.twinx()
    l2,=ax3_2.plot(dycyx0p860_140.RT.ScatA[Q_a2:Q_a1],mNo_Qdom[Q_a2:Q_a1],'k-')
    ax3.legend((l1,l2),(r'$Q^{\mu_s \mu_v}(\Theta)$',r'$Q^{\mu_s \mu_v}(\Theta)$/$Q^{\mu_s \mu_v}(\Theta_{cb})$'),loc='best',prop={'size':10})
    ax3.set_xlabel(r'Scattering Angle ($^{o}$)')
    fig3.show()
    #--------------------------------------------------------------------------
    # Domain avg Q and theoretical P12 comparison
    theP12=bulk.Mie.ang[p12_a1-1:p12_a2+2]
    theQ=dycyx0p860_140.RT.ScatA[Q_a2:Q_a1]
    mNo_P12=mNo_avP12[:,p12_a1-1:p12_a2+2]
    mNo_Q=mNo_Qdom[Q_a2:Q_a1]
    mNo_P12_itp=np.zeros((re.size,theQ.size),dtype='float')
    for i in np.arange(0,re.size):
        itpF=itp.interp1d(theP12,mNo_P12[i,:])#interpolated avP12 to match with observed scattering angles
        mNo_P12_itp[i,:]=itpF(theQ)

    fig4,ax4=plt.subplots()
    fig4_ttl=dycyx0p860_140.RT.fname.split('.',1)[0]+'_Cloudbow_max_normalized_-Q_and_-P12_obs_scat_interp'
    fig4.suptitle("\n".join(wrap(fig4_ttl,65)))
    ax4.plot(theQ,mNo_Q,'k-',linewidth=2.0,label=r'$Q^{\mu_s \mu_v}$')
    for i in np.arange(0,re.size):
        ax4.plot(theQ,mNo_P12_itp[i,:],label=r'$r_e=$%.1f'%re[i])
    ax4.legend(loc='best',prop={'size':10})
    ax4.set_xlabel(r'Scattering Angle ($^{o}$)')
    fig4.show()
    #--------------------------------------------------------------------------
    #Domain avg effective radius.    
    cr=np.corrcoef(mNo_P12_itp,mNo_Q)
    re_i=np.argmax(cr[re.size,:-1])#retrieved re index
    re[re_i]
def Example01():
    '''
    Polarimetric retrievals of DYCOMS2 LES RT simulations by using my method.
    '''
    dycyx2p13_140=LES_case('DYCOMS2_dharma_008036_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH2e6.hdf5','/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/')
    dyc2_dharma=DHARMA_onmp('dharma_008036.cdf','/umbc/xfs1/zzbatmos/users/charaj1/LES_simulations/DYCOMS/Ackerman_DYCOMS2_25bins/');dyc2_dharma.readDHARMA()    

    cloud_mie=MieSet('Cloud_mie_470_860_2p13')
    cloud_mie.readMie()
    psd=PSDs('mod_gamma',np.asarray(cloud_mie.d,dtype=float))
    #--------------------------------------------------------------------------
    #P12 LUT generation
#    re=np.array([9,10,12,15,20,30])
    re=np.arange(5,24,0.25)
    ve=0.02
    avP12_860=np.zeros((re.size,cloud_mie.ang.size),dtype=float)
    for i in np.arange(0,re.size):
        psd.set_mod_gamma_norm(re[i],ve)
        bulk=bulk_Mie('bulk_Mie_470_860_2p13',psd=psd,Mie=cloud_mie)
        bulk.cal_bulk_Pmat()
        avP12_860[i,:]=bulk.bulk_P12[:,2]

    # Q observations from RT simulations
#    ixx=find_peaks_cwt(data,np.arange(1,50))
    P12_in=avP12_860
    MieA=bulk.Mie.ang
    Q_in=dycyx2p13_140.RT.MeanPRad[:,:,:,1]
    ScatA=dycyx2p13_140.RT.ScatA
    SZA=180-dycyx2p13_140.RT.SZA
    VZA=dycyx2p13_140.RT.VZA
    re_in=re

    ret_re=retrieve_pol(P12_in,re_in,MieA,Q_in,ScatA,VZA,SZA)
            
    fig5,ax5=plt.subplots()
    fig5_ttl=dycyx2p13_140.RT.fname.split('.',1)[0]+'_whole_domain_pol_ret'
    fig5.suptitle("\n".join(wrap(fig5_ttl,65)))
    ctf=ax5.contourf(dycyx2p13_140.xcens,dycyx2p13_140.ycens,ret_re,np.linspace(10,20,100),extend='both')
    fig5.colorbar(ctf,ax=ax5,ticks=np.arange(10,20.1,2),label=r'$r_e$ ($\mu m$)')
    ax5.set_xlabel('km')
    ax5.set_ylabel('km')
    fig5.show()
def getGuess(cname,method,rTyp,tail):
    '''
    To read abc for domain mean retrievals from previous results. 
    '''
    filename=cname+'_'+method+'_mean'+tail
    abc=None
    if os.path.isfile(filename+'.pkl'):
        print('Domain mean retrievals exist: '+filename)
        if rTyp=='mean':
            print('Since \'mean\' retrievals being executed, an initial guess will not be used.')
        else:
            data=cpn.load_obj(filename)
            abc   =data['abc']
    else:
        print('Domain mean retrievals not found!')
    if abc is not None:
        print('Domain mean retrievals will be used as the initial guess.')
    elif not(rTyp == 'mean'):
        print('To increase the efficienty, do domain mean retrievals first to get an initial guess.')
    
    return abc
    
if __name__=='__main__':
    cpn.setup_figures(plt)
    band=str(sys.argv[1])   
    sza=int(sys.argv[2])    #Source Zenith Angle
    dim=str(sys.argv[3])    #RT vector/scalar
    method=str(sys.argv[4]) #Retrieval technique theory
    rTyp=str(sys.argv[5])   # 'mean'/'full'/'pixel'
    tail=str(sys.argv[6])   #Tail to the output file 
    spRT=False               #Use special RT runs with higher angular resolution
    mpi=True and (rTyp=='full')               #To use MPI version
    case=str(sys.argv[7])
    pBow=True               #Use primary bow

#    band='0p860'
#    sza=140 
#    dim='3D' 
#    method='Breon'
#    rTyp='full' # 'mean'/'full'/'pixel'
    X,Y=120,120 # Give pixel location if rTyp='pixel'
#    tail='_pol_retX%dY%d'%(X,Y)    
#    tail='_pol_ret'    
#    spRT=False
    #--------------------------------------------------------------------------
    #P12 LUT generation    
    P12Lib=Pol_ret_P12Lib(fname='Pol_ret_PM_library_0p860_2p13_V3.hdf5')
    P12Lib.loadP12Lib()
    
#    P12Lib.generateLib(re=np.arange(2,30,0.25),ve=np.array([0.01 , 0.02 , 0.03 , 0.04 , 0.05 ,0.06, 0.07, 0.1,0.11,0.125, 0.15,0.175, 0.2  , 0.225, 0.25 , 0.275]))
#    P12Lib.saveP12Lib()
    print(P12Lib.fname+' is being used')
    #--------------------------------------------------------------------------
    if not(spRT):
        # Q observations from RT simulations
        if case=='DY':
            LEScase=LES_case('DYCOMS2_dharma_008036_b'+band+'_MSCART_SZA'+str(sza)+'_SAA000_VAA000plus_NPH2e6.hdf5',\
                         '/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/',\
                         RT1Dname='DYCOMS2_dharma_008036_b'+band+'_MSCART_1D_bins_SZA'+str(sza)+'_SAA000_VAA000plus_NPH1e5.hdf5')
        elif case=='RC':
            LEScase=LES_case('RICO_dharma_005044_b'+band+'_MSCART_SZA'+str(sza)+'_SAA000_VAA000plus_NPH1e6.hdf5',\
                         '/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/RICO/',\
                         RT1Dname='RICO_dharma_005044_b'+band+'_MSCART_1D_bins_SZA'+str(sza)+'_SAA000_VAA000plus_NPH1e5.hdf5')
        elif case=='AC':
            LEScase=LES_case('ATEXc_dharma_007877_b'+band+'_MSCART_SZA'+str(sza)+'_SAA000_VAA000plus_NPH1e6.hdf5',\
                         '/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXc/',\
                         RT1Dname='ATEXc_dharma_007877_b'+band+'_MSCART_1D_bins_SZA'+str(sza)+'_SAA000_VAA000plus_NPH1e5.hdf5')
        elif case=='AP':
            LEScase=LES_case('ATEXp_dharma_013067_b'+band+'_MSCART_SZA'+str(sza)+'_SAA000_VAA000plus_NPH1e6.hdf5',\
                         '/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXp/',\
                         RT1Dname='ATEXp_dharma_013067_b'+band+'_MSCART_1D_bins_SZA'+str(sza)+'_SAA000_VAA000plus_NPH1e5.hdf5')
            
                
    '''
    Corelation coefficient method (yet to be investigated)
    ======================================================
    MieA=P12Lib.bulk_Mie_ang
    Q_in=LEScase.RT.MeanPRad[:,:,:,1]
    ScatA=LEScase.RT.ScatA
    SZA=180-LEScase.RT.SZA
    VZA=LEScase.RT.VZA
    re_in=P12Lib.re
    ve_in=P12Lib.ve
    P12_in=P12Lib.avP12['0p860']
    start=time.time()
    ret_re,ret_ve=retrieve_pol(P12_in,re_in,ve_in,MieA,Q_in,ScatA,VZA,SZA)
    end=time.time()
    print('%d mins elapsed'%(end/60-start/60))
    
    fig5,ax5=plt.subplots()
    fig5_ttl=LEScase.RT.fname.split('.',1)[0]+'_whole_domain_pol_ret'
    fig5.suptitle("\n".join(wrap(fig5_ttl,65)))
    ctf=ax5.contourf(LEScase.xcens,LEScase.ycens,ret_re,np.linspace(0,25,50),extend='both')
    fig5.colorbar(ctf,ax=ax5,ticks=np.arange(0,25.1,5),label=r'$r_e$ ($\mu m$)')
    ax5.set_xlabel('km')
    ax5.set_ylabel('km')
    fig5.show()
    '''
    
    '''
    Parametric curve fitting method
    ===============================
    '''
    if not(spRT):
        #Selecting 1D/3D RT
        if dim=='1D':
            RT=LEScase.RT1D
        elif dim=='3D':
            RT=LEScase.RT
    else:
        RT=POLCARTdset('LES','/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/For_Pols/DYCOMS2/')
        RT.readPOLCARThdf5('DYCOMS2_dharma_008036_b0p860_MSCART_SZA140_SAA000_VAA180Pol_NPH1e6.hdf5',dpath='/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/For_Pols/DYCOMS2/'+'results/b'+band+'/')
    #Parametric curvefitting method
    if rTyp=='mean':
        Q_in2= RT.MeanPRad[:,:,:,1].mean(axis=1).mean(axis=1)#Averaging whole domain
    elif rTyp=='pixel':
        Q_in2= RT.MeanPRad[:,X,Y,1]
        tail="X%dY%d"%(X,Y)+tail
    elif rTyp=='full':
        Q_in2= RT.MeanPRad[:,:,:,1]
    else:
        print('Error! Give valid rTyp!')
    try:
        cname=RT.fname.split('.',1)[0]
    except AttributeError:
        cname=(RT.fname[0].astype(str)).split('.',1)[0]
    savename=cname+'_'+method+'_'+rTyp+tail
    print(savename+' will be saved')

    muV=np.cos(np.deg2rad(RT.VZA))
    muS=np.cos(np.deg2rad(180-RT.SZA))
    obsSca=RT.ScatA
    P=Pmat(P12Lib.re,P12Lib.ve,P12Lib.bulk_Mie_ang,P12Lib.avP12['0p860'],obsSca,method=method,primaryBow=pBow)
    #gemet=4*(muS+muV)*0+1#No geometric correction
    gemet=4*(muS+muV)
    ygabc=getGuess(cname,method,rTyp,tail)   
    x=obsSca[P.Q_a1:P.Q_a2]
    start=time.time()
    if Q_in2.ndim==3:
        y=np.einsum('ijk,i->ijk',-Q_in2,4*gemet) 
        if mpi:
            mpi_data={'x':x,'y':y,'ygabc':ygabc,'P':P,'savename':savename}
            cpn.save_obj(mpi_data,'mpi_data_'+savename,rp=True)
            print('Running with MPI...')
            os.system('mpirun -n 16 python pol_ret_mpi.py '+'mpi_data_'+savename)
        else:   
            if P.method=='Breon':
                from cpnRetrievalslib import fitBreon_noRsq as do_fitting
            elif P.method=='BreonMod':
                from cpnRetrievalslib import fit_BreonMod as do_fitting
            ret_Re=np.zeros((Q_in2.shape[1],Q_in2.shape[2]),dtype=float)
            ret_Ve=np.zeros((Q_in2.shape[1],Q_in2.shape[2]),dtype=float)
            Qls   =np.zeros((Q_in2.shape[1],Q_in2.shape[2]),dtype=float)
            Rsq   =np.zeros((Q_in2.shape[1],Q_in2.shape[2]),dtype=float)
            abc   =np.zeros((Q_in2.shape[1],Q_in2.shape[2],3),dtype=float)
            yAll  =np.zeros((Q_in2.shape[1],Q_in2.shape[2],x.size))
            c=0
#            os.system("rm Rsqs_"+savename+".dat")
#            os.system("rm Rsqs_max_"+savename+".dat")
#            os.system("echo \"Rsqs\t Re\t Ve\" > Rsqs.dat")
#            os.system("echo \"Rsqs_max\t Re\t Ve\" > Rsqs_max.dat")
            for i in np.arange(Q_in2.shape[1]):
                for j in np.arange(Q_in2.shape[2]):
                    ret_Re[i,j],ret_Ve[i,j],abc[i,j,:],Qls[i,j],Rsq[i,j]=do_fitting(x,np.squeeze(y[P.Q_a1:P.Q_a2,i,j]),P,ygabc)
                    yAll[i,j,:]=y[P.Q_a1:P.Q_a2,i,j]
                    c+=1
                    pc=c/Q_in2.shape[1]/Q_in2.shape[2]*100.0
                    tm=time.time()/60/60-start/60/60
                    print("\r%0.2f%% %0.2f hours remaining ..."%(pc,tm/pc*100-tm),end=" ")#
    elif Q_in2.ndim==1:
        if P.method=='Breon':
            from cpnRetrievalslib import fitBreon as do_fitting
        elif P.method=='BreonMod':
            from cpnRetrievalslib import fit_BreonMod as do_fitting
        y=(-Q_in2*gemet)[P.Q_a1:P.Q_a2]
        os.system("rm Rsqs.dat")
        os.system("rm Rsqs_max.dat")
        os.system("echo \"Rsqs\t Re\t Ve\" > Rsqs.dat")
        os.system("echo \"Rsqs_max\t Re\t Ve\" > Rsqs_max.dat")
        ret_Re,ret_Ve,abc,Qls,Rsq=do_fitting(x,y,P,ygabc)
        os.system("mv Rsqs.dat Rsqs_"+savename+".dat")
        os.system("mv Rsqs_max.dat Rsqs_max_"+savename+".dat")
        P.set_reve(ret_Re,ret_Ve)
        yAll=y
        fig2,ax2=plt.subplots()
        ax2.plot(x,y,'k.-')
        ax2.plot(x,P.imitateF(x,*abc))
        fig2.show()
    end=time.time()
    print('%0.2f mins elapsed'%(end/60-start/60))
    
    if not(mpi):
        data={'ret_Re':ret_Re,'ret_Ve':ret_Ve,'abc':abc,'Qls':Qls,'Rsq':Rsq,'yAll':yAll,'x':x}
        print(abc)
        cpn.save_obj(data,savename,rp=True)    



    '''
    ***************************************************************************
    Plots
    '''
    '''
    #Checking the fit
    P.set_reve(ret_Re,ret_Ve)
    fig6,ax6=plt.subplots(figsize=(8,4))
    fig6_ttl=savename+'_checking_fit'
    fig6.suptitle("\n".join(wrap(fig6_ttl,70)),size=12)
    ax6.plot(RT.ScatA[P.Q_a1:P.Q_a2],P.imitateF(RT.ScatA[P.Q_a1:P.Q_a2],*abc),'g--',label='fit')
    ax6.plot(RT.ScatA[P.Q_a1:P.Q_a2],(-Q_in2*4*(muV+muS))[P.Q_a1:P.Q_a2],'k.-')
#    ax6.plot(RT.ScatA[P.Q_a1:P.Q_a2],-P.getP(RT.ScatA[P.Q_a1:P.Q_a2]),'r--')
    ax6.text(0.5,0.6,r'$r_e$ %0.2f, $v_e$ %0.2f, $R^2$ %0.2f'%(ret_Re,ret_Ve,Rsq),transform=ax6.transAxes)
    ax6.set_xlabel('Scattering Angle')
    ax6.set_ylabel(r'$R_p^*$')
    ax6.legend(frameon=False,loc='best')
    fig6.tight_layout(rect=[0,0,1,0.96])
    fig6.show()

    fig5,ax5=plt.subplots()
    fig5_ttl=RT.fname[0].astype(str).split('.',1)[0]+'_whole_domain_pol_ret_curve_fitting'
    fig5.suptitle("\n".join(wrap(fig5_ttl,65)))
    ctf=ax5.contourf(LEScase.xcens,LEScase.ycens,ret_Re,np.linspace(10,20,100),extend='both')
    fig5.colorbar(ctf,ax=ax5,ticks=np.arange(10,20.1,2),label=r'$r_e$ ($\mu m$)')
    ax5.set_xlabel('km')
    ax5.set_ylabel('km')
    fig5.show()
    '''
