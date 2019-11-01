#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Wed Apr 26 21:32:36 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To convert 3D field files to individual column field files to run MSCART separately.
09/01/2017:
    To handle fractal cloud case. Trying to make it compatible with any LES_field object
05/14/2018:
    For the LES cases. Uses the "actual 1D feature" of MSCART
"""
import numpy as np
import cpnLES_MSCARTlib as lib


band='2p13'
d3D=lib.LES_field('ATEXc_dharma_007877_b2p13.nc',dpath='../')
d3D.readLES_field()
#i1D=lib.LES_field('OP_A5.nc')
#i1D.readLES_field()

def dobins():
    for i in np.arange(0,d3D.ny):
        for j in np.arange(0,d3D.nx):
            fobj=lib.LES_field(d3D.fname.split('.',1)[0]+'y'+str(i)+'x'+str(j)+'.nc')
            fobj.nxp=2
            fobj.nyp=2
            fobj.nzp=d3D.nzp
            fobj.nx=1
            fobj.ny=1
            fobj.nz=d3D.nz
            fobj.nz3=1#d3D.nz3 #Number of horizontally inhomogeneous layers
            fobj.iz3l=1
            fobj.npar1=d3D.npar3
            fobj.npar3=2#d3D.npar3
            fobj.nkd=d3D.nkd
            fobj.nstokes=d3D.nstokes
            fobj.nang=d3D.nang
            fobj.npf=d3D.npf
            fobj.npsfc=d3D.npsfc
            fobj.xgrd=np.array([0,1])
            fobj.ygrd=np.array([0,1])
        
            fobj.zgrd=d3D.zgrd
            fobj.extp1d=d3D.extp3d[:,:,i,j]
            fobj.omgp1d=d3D.omgp3d[:,:,i,j]
            fobj.jpfp1d=d3D.jpfp3d[:,:,i,j]
            fobj.rpfp1d=d3D.rpfp3d[:,:,i,j]
            fobj.extp3d=np.zeros((fobj.npar3,fobj.nz3,1,1),dtype=float)
            fobj.omgp3d=np.ones ((fobj.npar3,fobj.nz3,1,1),dtype=float)
            fobj.jpfp3d=np.zeros((fobj.npar3,fobj.nz3,1,1),dtype=float)
            fobj.rpfp3d=np.zeros((fobj.npar3,fobj.nz3,1,1),dtype=float)
            fobj.absg1d=d3D.absg1d
            fobj.wkd=d3D.wkd
            fobj.scatAng=d3D.scatAng
            fobj.PMatVal=d3D.PMatVal
            fobj.tsfc2d=d3D.tsfc2d[i,j]
            fobj.jsfc2d=d3D.jsfc2d[i,j]
            fobj.psfc2d=d3D.psfc2d[i,j,:]

            fobj.writeLES_field(fobj.fname,dpath='LESb'+band+'_bins/')
	

dobins()
