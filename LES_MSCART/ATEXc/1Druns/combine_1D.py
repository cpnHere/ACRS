#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Sat Sep  2 13:16:42 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To combine fractal clouds 1-D MSCART RT simulations
12/18/2017:
    Harigiye nathi ewa ganan karanawa dan.
03/29/2017:
    Second "missing data" handling exception implemented to look for invalid high radiances.
"""
import cpnLES_MSCARTlib as MSlib
import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

t1=time.time()

def iqu_cb(fig,ctf,ax,ticks=None,orientation='horizontal',label='label',pad=0.2):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=pad)
    if ticks==None:    
        fig.colorbar(ctf, cax=cax,orientation=orientation,label=label)
    else:
        fig.colorbar(ctf, cax=cax,ticks=ticks,orientation=orientation,label=label)

cname=sys.argv[1]
pc_format=sys.argv[2]
action=sys.argv[3]
sza=cname.split('_SZA',1)[1].split('_',1)[0]
band=cname.split('_b',1)[1].split('_',1)[0]
nmlpath='LESb'+band+'_bins/'#*.nml path for 1D runs
fdpath='results/LESb'+band+'_bins/'

dsNew=MSlib.POLCARTdset('dsNew','/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/LES_MSCART/ATEXc/1Druns/LESb'+band+'_bins')
dsNew.readPOLCARThdf5(cname.split('NPH',1)[0]+'NPH1e6.hdf5','../results/b'+band+'/')

d3D=MSlib.LES_field('ATEXc_dharma_007877_b'+band+'.nc',dpath='../')
d3D.readLES_field()

dsNew.MeanPRad[:]=np.nan
have_all_files=True
second_test=False
missing_count=0
for i in np.arange(0,dsNew.MeanPRad.shape[1],1):
    for j in np.arange(0,dsNew.MeanPRad.shape[2],1):
        fn1=cname.split('_MSCART',1)[0]
        fn2='_MSCART'+cname.split('_MSCART',1)[1]
        fname=fn1+'y'+str(j)+'x'+str(i)+fn2+'.nc'
        ds1D=MSlib.POLCARTdset('ds1D',nmlpath)
        print(fname)
        try:
            ds1D.readMSCARTplus(fname,fdpath=fdpath,clm=True,prnt=False,step=True)
            if j==0 and i==0:
                dsNew.SZA=ds1D.SZA
                dsNew.SAA=ds1D.SAA
                dsNew.NPH=ds1D.NPH
                dsNew.Nbat=ds1D.Nbat
                dsNew.MeanTiming=ds1D.MeanTiming
                dsNew.VZA=ds1D.VZA
                dsNew.VAA=ds1D.VAA
                dsNew.cc3D=ds1D.cc3D
                dsNew.nmldpath=ds1D.nmldpath
                dsNew.ScatA=ds1D.ScatA
            dsNew.MeanPRad[:,i,j,:]=ds1D.MeanPRad
            dsNew.RMSEPRad[:,i,j,:]=ds1D.RMSEPRad
            second_test=True
    
        except FileNotFoundError:
            print(fname+' RuntimeError: Possibly does not exist!!!')
            if missing_count==0:
                missing_fnames=np.array(fname)
            else:
                missing_fnames=np.vstack((missing_fnames,fname))
            missing_count=missing_count+1
            if have_all_files:
                have_all_files=False
        
        #Checking for another set of missing file (having very high dummy values as radiances)
        if second_test:
            rad_val=dsNew.MeanPRad[10,i,j,0]
            if rad_val>10:
                print(fname+' Invalid high radiance value found!!')
                if missing_count==0:
                    missing_fnames=np.array(fname)
                else:
                    missing_fnames=np.vstack((missing_fnames,fname))
                missing_count=missing_count+1
                if have_all_files:
                    have_all_files=False
   
NPH=cname.split('NPH',1)[1]
if have_all_files:
    print('No missing files found!!')
    dsNew.fname=dsNew.fname.split('_MSCART',1)[0]+'_MSCART'+'_1D_bins'+dsNew.fname.split('_MSCART',1)[1].split('.',1)[0]
    dsNew.savePOLCARThdf5(dsNew.fname.split('NPH',1)[0]+'NPH'+NPH+'.hdf5',dpath=fdpath,pc_format=pc_format,action=action)

else:
    print(missing_fnames)
    missing_job=np.zeros(missing_fnames.size,dtype=int)
    missing_fnames=missing_fnames.reshape(missing_fnames.size)
    for i in np.arange(0,missing_fnames.size):
        yi=int(missing_fnames[i].split('y',1)[1].split('x',1)[0])
        xj=int(missing_fnames[i].split('y',1)[1].split('x',1)[1].split('_MSCART',1)[0])
        missing_job[i]=yi*dsNew.MeanPRad.shape[1]+xj
#    for x in missing_job: sys.stdout.write(str(x)+' ')
    umissing_job=np.unique(missing_job)
    missing=np.insert(umissing_job,0,umissing_job.size)
    datName=dsNew.fname.split('_MSCART',1)[0]+'_MSCART'+'_1D_bins'+dsNew.fname.split('_MSCART',1)[1].split('.',1)[0]
    datName=datName.split('NPH',1)[0]+'NPH'+NPH
    np.savetxt(datName+'_missingRuns.dat',missing,fmt='%d')
    print('%d files are missing!!'%missing_count)
    print('Missing job numbers saved to: '+datName+'_missingRuns.dat' )
    
def contourf_Iqu(fig,ax,obj,v=None,ticks=None,cmap=None):
    #pb_VZA_ang: primarybow VZA
    #xb_VZA_ang: supernumerary VZA
    xcens=(d3D.xgrd[1:]+d3D.xgrd[0:-1])/2
    ycens=(d3D.ygrd[1:]+d3D.ygrd[0:-1])/2
    if v==None:
        v = np.linspace(-0.05,0.15,num=256)
    else:
        v=v
    ctf=ax[0].contourf(xcens,ycens,obj.MeanPRad[61,:,:,0] ,v,extend='both',cmap=cmap);iqu_cb(fig,ctf,ax[0],label='R_I')
    ctf=ax[1].contourf(xcens,ycens,-obj.MeanPRad[61,:,:,1],v,extend='both',cmap=cmap);iqu_cb(fig,ctf,ax[1],label='-Q_I')
    ctf=ax[2].contourf(xcens,ycens,obj.MeanPRad[61,:,:,2] ,v,extend='both',cmap=cmap);iqu_cb(fig,ctf,ax[2],label='U_I')

t2=time.time()
print('%0.2f mins elapsed'%((t2-t1)/60))
fig2,ax2=plt.subplots(1,3,subplot_kw={'aspect':'equal'})
fig2_ttl=contourf_Iqu(fig2,ax2,dsNew)
fig2.show()
