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
"""
import numpy as np
import cpnLES_MSCARTlib as lib


band='0p865'
d3D=lib.LES_field('fractal_cld_b865re12ve05_x40km.nc',dpath='../')
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
            fobj.nz3=d3D.nz3
            fobj.iz3l=d3D.iz3l
            fobj.npar1=d3D.npar1
            fobj.npar3=d3D.npar3
            fobj.nkd=d3D.nkd
            fobj.nstokes=d3D.nstokes
            fobj.nang=d3D.nang
            fobj.npf=d3D.npf
            fobj.npsfc=d3D.npsfc
            fobj.xgrd=np.array([0,1])
            fobj.ygrd=np.array([0,1])
        
            fobj.zgrd=d3D.zgrd
            fobj.extp1d=d3D.extp1d
            fobj.omgp1d=d3D.omgp1d
            fobj.jpfp1d=d3D.jpfp1d
            fobj.rpfp1d=d3D.rpfp1d
            fobj.extp3d=d3D.extp3d[:,:,i,j]
            fobj.omgp3d=d3D.omgp3d[:,:,i,j]
            fobj.jpfp3d=d3D.jpfp3d[:,:,i,j]
            fobj.rpfp3d=d3D.rpfp3d[:,:,i,j]
            fobj.absg1d=d3D.absg1d
            fobj.wkd=d3D.wkd
            fobj.scatAng=d3D.scatAng
            fobj.PMatVal=d3D.PMatVal
            fobj.tsfc2d=d3D.tsfc2d[i,j]
            fobj.jsfc2d=d3D.jsfc2d[i,j]
            fobj.psfc2d=d3D.psfc2d[i,j,:]

            fobj.writeLES_field(fobj.fname,dpath='fracb'+band+'_bins/')
	

dobins()
