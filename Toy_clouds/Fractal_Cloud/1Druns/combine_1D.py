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

sza='140'
band ='0p865'
bandx='865'
nmlpath='fracb'+band+'_bins/'#*.nml path for 1D runs
fdpath='results/fracb'+band+'_bins/'

dsNew=MSlib.POLCARTdset('dsNew','/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/Toy_clouds/Fractal_Cloud/1Druns/fracb'+band+'_bins')
dsNew.readPOLCARThdf5('fractal_cld_b'+bandx+'re12ve05_x40km_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5','../results/b'+band+'/')
#dsNew.readPOLCARThdf5('fractal_StCu_ACA_b'+band+'_x40km_H0p5_AOT0p25_H0p5_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.hdf5','../results_aca/b'+band+'/')
#lib.collectColsMSCART('OP_dharma_008036_full_3_26_binsy9x9_MSCART_SZA060_SAA030_VAA000plus_NPH1e5.nc','results/OP_dharma_008036_full_3_26_bins/',dsNew,'MSCARToutTest.hdf5')

d3D=MSlib.LES_field('fractal_cld_b865re12ve05_x40km.nc',dpath='../')
d3D.readLES_field()

dsNew.MeanPRad[:]=np.nan
have_all_files=True
second_test=False
missing_count=0
for j in np.arange(0,dsNew.MeanPRad.shape[1],1):
    fname='fractal_cld_b'+bandx+'re12ve05_x40kmy0x'+str(j)+'_MSCART_SZA'+sza+'_SAA000_VAA000plus_NPH1e5.nc'
    ds1D=MSlib.POLCARTdset('ds1D',nmlpath)
    print(fname)
    try:
        ds1D.readMSCARTplus(fname,fdpath=fdpath,clm=True,prnt=False,step=True)
        if j==0:
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
        dsNew.MeanPRad[:,j,:]=ds1D.MeanPRad
        dsNew.RMSEPRad[:,j,:]=ds1D.RMSEPRad
        second_test=True

    except RuntimeError:
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
        rad_val=dsNew.MeanPRad[10,j,0]
        if rad_val>10:
            print(fname+' Invalid high radiance value found!!')
            if missing_count==0:
                missing_fnames=np.array(fname)
            else:
                missing_fnames=np.vstack((missing_fnames,fname))
            missing_count=missing_count+1
            if have_all_files:
                have_all_files=False
   
    
if have_all_files:
    print('No missing files found!!')
    dsNew.fname=dsNew.fname.split('_MSCART',1)[0]+'_1D_bins_MSCART'+dsNew.fname.split('_MSCART',1)[1].split('.',1)[0]+'.hdf5'
    dsNew.savePOLCARThdf5(dsNew.fname,dpath='results/fracb'+band+'_bins/',pc_format='fractal_1D',action='Spot3D_paper_revision')

else:
    print('%d files are missing!!'%missing_count)
    print(missing_fnames)
    missing_fname_xbin=np.zeros(missing_fnames.size,dtype=int)
    missing_fnames=missing_fnames.reshape(missing_fnames.size)
    for i in np.arange(0,missing_fnames.size):
        missing_fname_xbin[i]=int(missing_fnames[i].split('y0x',1)[1].split('_MSCART',1)[0])
    print(missing_fname_xbin)

xcens=(d3D.xgrd[1:]+d3D.xgrd[0:-1])/2
def contourf_Iqu(fig,ax,iqu,pb_VZA_ang,sb_VZA_ang,obj_VZA,obj_MeanPRad,v=None,ticks=None,cmap=None):
    #pb_VZA_ang: primarybow VZA
    #xb_VZA_ang: supernumerary VZA
    if iqu==0:
        sgn=1;iq_lab='I'
    elif iqu==1:
        sgn=-1;iq_lab='-Q'
    fig_ttl=iq_lab+'_contours'
    if v==None:
        v = np.linspace(-0.05,0.15,num=256)
    else:
        v=v
    ctf=ax.contourf(xcens,obj_VZA,sgn*obj_MeanPRad[:,:,iqu],v,extend='both',cmap=cmap)
    if ticks==None:
        fig.colorbar(ctf,ax=ax,orientation='vertical',ticks=np.arange(-0.05,0.16,0.05))
    else:
        fig.colorbar(ctf,ax=ax,orientation='vertical',ticks=ticks)
    if np.size(pb_VZA_ang)>1:
        for i in range(0,np.size(pb_VZA_ang)):
            ax.plot(xcens,np.zeros_like(xcens)+pb_VZA_ang[i],'k--')
    else:
        ax.plot(xcens,np.zeros_like(xcens)+pb_VZA_ang,'k--')
    if np.size(sb_VZA_ang)>1:
        for i in range(0,np.size(sb_VZA_ang)):
            ax.plot(xcens,np.zeros_like(xcens)+sb_VZA_ang[i],'k--')
    else:
        ax.plot(xcens,np.zeros_like(xcens)+sb_VZA_ang,'k--')
    ax.set_ylabel('VZA',size=10)
    return fig_ttl

fig2,ax2=plt.subplots()
fig2_ttl=contourf_Iqu(fig2,ax2,1,dsNew.VZA[50],dsNew.VZA[60],dsNew.VZA,dsNew.MeanPRad)
ax2.set_xlabel('x (km)',size=10,y=-0.5)
ax2.set_title(r'$-Q_{1D}$')
fig2.show()
