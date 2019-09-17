#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Sep  2 11:41:34 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Polarimetric retrievals for Toy clouds
"""
import numpy as np
import matplotlib.pyplot as plt
import time, os
import cpnCommonlib as cpn
from cpnRetrievalslib import  Pol_ret_P12Lib,Pmat
from cpnLES_MSCARTlib import POLCARTdset,LES_field
from pol_cloud_ret import getGuess

def read_step_cloud(band,sza):
    '''
    To read 3D and 1D RT simulations of a step cloud case
    '''
#    band='865'
#    sza='120'
#    mainDpath='/home/cpnhere/spot3DPaper/Data/'   #To work on local machine
    mainDpath='/umbc/xfs1/zzbatmos/users/charaj1/'#To work on Taki
    nmlpath=mainDpath+'Toy_clouds_MSCART/Step_Cloud/'
    dsp120b865=POLCARTdset('dsp120b865',nmlpath)
    dsp120b865.readMSCARTplus('step_cld_b'+band+'re10ve02_x15km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e6.nc',\
                         fdpath=nmlpath+'results/stepb'+band+'/',step=True)
    #column runs (1D)
    cfree=POLCARTdset('cfree',nmlpath+'1Druns/')
    cfree.readMSCARTplus('step_cld_b'+band+'re10ve02_x1Dex0p1_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc',\
                           fdpath=nmlpath+'1Druns/results/stepb'+band+'/',clm=True,step=True)

    cstep=POLCARTdset('cstep',nmlpath+'1Druns/')
    cstep.readMSCARTplus('step_cld_b'+band+'re10ve02_x1Dex10_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc',\
                           fdpath=nmlpath+'1Druns/results/stepb'+band+'/',clm=True,step=True)

    step1Dto3D=POLCARTdset('step1Dto3D','generated_using_'+cfree.fdpath+cfree.fname)
    step1Dto3D.fname=cfree.fname.split('.',1)[0]+'_1Dto3D.nc'
    step1Dto3D.fdpath=cfree.fdpath
    step1Dto3D.VZA=cfree.VZA
    step1Dto3D.ScatA=cfree.ScatA
    step1Dto3D.MeanPRad=np.zeros_like(dsp120b865.MeanPRad)
    step1Dto3D.SZA=cfree.SZA
    for i in range(0,200):
        step1Dto3D.MeanPRad[:,i,:]=cfree.MeanPRad
    for i in range(200,1200):
        step1Dto3D.MeanPRad[:,i,:]=cstep.MeanPRad
    for i in range(1200,1500):
        step1Dto3D.MeanPRad[:,i,:]=cfree.MeanPRad
    #field file
    dstep=LES_field(dsp120b865.fname.split('_MSCART',2)[0]+'.nc',dpath=nmlpath)#DYCOM field file for MSCART
    dstep.readLES_field()
    return dsp120b865,step1Dto3D,dstep
class Fractal_CLD(object):
    def __init__(self,sza):
        '''
        sza: SZA string
        '''
        #mainDpath='/home/cpnhere/spot3DPaper/Data/'
        mainDpath='/umbc/xfs1/zzbatmos/users/charaj1/' #To work on Taki

        self.d865SZAx=POLCARTdset('d865SZAx',mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/')
        self.d865SZAx.readPOLCARThdf5('fractal_cld_b865re12ve05_x40km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                                      dpath=mainDpath+'/Toy_clouds_MSCART/Fractal_Cloud/results/fracb865/')
        self.d2p13SZAx=POLCARTdset('d2p13SZAx',mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/')
        self.d2p13SZAx.readPOLCARThdf5('fractal_cld_b2p13re12ve05_x40km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                                       dpath=mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/results/fracb2p13/')
        
        self.dfrac865=LES_field(self.d865SZAx.fname.split('_MSCART',2)[0]+'.nc',dpath=mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/')#DYCOM field file for MSCART
        self.dfrac865.readLES_field()
        self.dfrac2p13=LES_field(self.d2p13SZAx.fname.split('_MSCART',2)[0]+'.nc',dpath=mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/')#DYCOM field file for MSCART
        self.dfrac2p13.readLES_field()
        
        self.frac1D865=POLCARTdset('frac1D865',mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb0p865_bins/')
        self.frac1D865.readPOLCARThdf5('fractal_cld_b865re12ve05_x40km_1D_bins_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                                           dpath=mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/1Druns/results/fracb0p865_bins/')
        self.frac1D2p13=POLCARTdset('frac1D2p13',mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb2p13_bins/')
        self.frac1D2p13.readPOLCARThdf5('fractal_cld_b2p13re12ve05_x40km_1D_bins_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5',\
                                      dpath=mainDpath+'Toy_clouds_MSCART/Fractal_Cloud/1Druns/results/fracb2p13_bins/')
        self.xcens=(self.dfrac865.xgrd[1:]+self.dfrac865.xgrd[0:-1])/2        
    def find_3DimpFac(self,VZAi,R2p13_th,LUT_SWIRre=None,LUT_VNIRre=None):
        '''
        Find 3D impact factor
            f=(3D-1D)/1D
        VZAi: viewing zenith angle index (usually 61 for zero degrees)
        R2p13_th: R2p13 threshold 
        --------------------------
        Give LUT_SWIRre and LUT_VNIRre to calculate separate F_0p860
        LUT_SWIRre,LUT_VNIRre is the values of constant re line on the NJK LUT space
        '''
        self.f_R2p13 =(self.d2p13SZAx.MeanPRad[VZAi,:,0]-self.frac1D2p13.MeanPRad[VZAi,:,0])/self.frac1D2p13.MeanPRad[VZAi,:,0]
        self.f_R0p860=(self.d865SZAx.MeanPRad[VZAi,:,0]-self.frac1D865.MeanPRad[VZAi,:,0])/self.frac1D865.MeanPRad[VZAi,:,0]         
        self.F_R2p13 =(self.d2p13SZAx.MeanPRad[VZAi,:,0]-R2p13_th)/R2p13_th
        if LUT_SWIRre is not None:
            VNIRbias=np.sign(self.F_R2p13)*np.interp(abs(self.F_R2p13)*R2p13_th,LUT_SWIRre,LUT_VNIRre)
            self.F_R0p860=VNIRbias/(self.d865SZAx.MeanPRad[VZAi,:,0]-VNIRbias)
            if (self.d865SZAx.MeanPRad[VZAi,:,0]-VNIRbias).any(0):
                print('Zero 1D reflectance !!!!')


if __name__=='__main__':
    cpn.setup_figures(plt)
    '''
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
    '''
    #Input specifications: See run slurm for details
    band='865'
    sza='120'
    dim='3D'    #RT vector/scalar
    rTyp='full' # 'mean'/'full'/'pixel'
    tail=''     #Tail to the output file
    X=40        #Pixel location if rTyp='pixel'
    method='Breon' #Retrieval technique theory
    pBow=True      #Use primary bow
    case='step' #step/fractal
    mpi=False
    #--------------------------------------------------------------------------
    #P12 LUT generation    
    P12Lib=Pol_ret_P12Lib(fname='Pol_ret_PM_library_0p860_2p13_V3.hdf5')
    P12Lib.loadP12Lib()
    
#    P12Lib.generateLib(re=np.arange(2,30,0.25),ve=np.array([0.01 , 0.02 , 0.03 , 0.04 , 0.05 ,0.06, 0.07, 0.1,0.11,0.125, 0.15,0.175, 0.2  , 0.225, 0.25 , 0.275]))
#    P12Lib.saveP12Lib()
    print(P12Lib.fname+' is being used')
    #--------------------------------------------------------------------------
    if case=='step':
        dsp,dsp1D,dstep=read_step_cloud(band,sza)
    elif case=='fractal':
        fc=Fractal_CLD(sza)
        if band=='865':
            dsp,dsp1D,dstep=fc.d865SZAx,fc.frac1D865,fc.dfrac865
            
    if dim=='1D':
        RT=dsp1D
    elif dim=='3D':
        RT=dsp
    
    if rTyp=='mean':
        Q_in2= RT.MeanPRad[:,:,1].mean(axis=1)#Averaging whole domain
    elif rTyp=='pixel':
        Q_in2= RT.MeanPRad[:,X,1]
        tail="X%d"%(X)+tail
    elif rTyp=='full':
        Q_in2= RT.MeanPRad[:,:,1]
    else:
        print('Error! Give valid rTyp!')

    try:
        cname=RT.fname.split('.',1)[0]
    except AttributeError:
        cname=(RT.fname[0].astype(str)).split('.',1)[0]
    savename=cname+'_'+method+'_'+rTyp+tail
    muV=np.cos(np.deg2rad(RT.VZA))
    muS=np.cos(np.deg2rad(180-RT.SZA))
    obsSca=RT.ScatA
    P=Pmat(P12Lib.re,P12Lib.ve,P12Lib.bulk_Mie_ang,P12Lib.avP12['0p860'],obsSca,method=method,primaryBow=pBow)
    #gemet=4*(muS+muV)*0+1#No geometric correction
    gemet=4*(muS+muV)
    ygabc=getGuess(cname,method,rTyp,tail)   
    x=obsSca[P.Q_a1:P.Q_a2]
    start=time.time()

    if P.method=='Breon':
        from cpnRetrievalslib import fitBreon as do_fitting
    if Q_in2.ndim==1:
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
    elif Q_in2.ndim==2:
        y=np.einsum('ij,i->ij',-Q_in2,4*gemet) 
        ret_Re=np.zeros((Q_in2.shape[1]),dtype=float)
        ret_Ve=np.zeros((Q_in2.shape[1]),dtype=float)
        Qls   =np.zeros((Q_in2.shape[1]),dtype=float)
        Rsq   =np.zeros((Q_in2.shape[1]),dtype=float)
        abc   =np.zeros((Q_in2.shape[1],3),dtype=float)
        yAll  =np.zeros((Q_in2.shape[1],x.size))
        c=0
        for i in np.arange(Q_in2.shape[1]):
            ret_Re[i],ret_Ve[i],abc[i,:],Qls[i],Rsq[i]=do_fitting(x,np.squeeze(y[P.Q_a1:P.Q_a2,i]),P,ygabc)
            yAll[i,:]=y[P.Q_a1:P.Q_a2,i]
            c+=1
            pc=c/Q_in2.shape[1]*100.0
            tm=time.time()/60/60-start/60/60
            print("\r%0.2f%% %0.2f hours remaining ..."%(pc,tm/pc*100-tm),end=" ")#
    end=time.time()
    print('%0.2f mins elapsed'%(end/60-start/60))


    if not(mpi):
        data={'ret_Re':ret_Re,'ret_Ve':ret_Ve,'abc':abc,'Qls':Qls,'Rsq':Rsq,'yAll':yAll,'x':x}
        print(abc)
        cpn.save_obj(data,savename,rp=False)    
