#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri Jul 10 20:58:06 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************

To run DISORT for Kokhanovsky cloud case
"""

import numpy as np
from scipy.special import legendre
import os

if __name__ == "__main__":
    
    Data = np.loadtxt('P11/Zhibo_P11/Kokhanovsky_benchmark_cloud.PDA',skiprows=6)
    re = 1 #um
    ve = 0.3
    Ang = Data[:,0]
    P11 = Data[:,1]
    Mu = np.cos(np.radians(Ang))
    
    # check if the phase function is normalized
    C =-0.5*np.trapz(P11,Mu) # negative sign is because the integration should be from -1 to 1, but Mu is from 1 to -1
    print('normalization',C)
    if C>0.99999:
        print("Phase function is normalized.")
    else:
        print("Warning!!: Phase function is not normalized!!")
    
    print('Computing Legendre polynomials....')
    N=np.arange(0,900)
    g = []
    for n in N:
        f =  legendre(n)
        L = f(Mu)
        g.append(-0.5*np.trapz(P11*L,Mu))
    g = np.array(g)
    
    print('Writing input file...')
    COT = 5.0
    MU0 = np.cos(np.deg2rad(60))
    input_content=np.concatenate([[COT],[MU0],g])
    jobid='Kokhanovsky'+'_ve%0.1fre%0.1fCOT%0.1fMU0%0.4f'%(ve,re,COT,MU0)
    np.savetxt(jobid.replace('.','p')+'_inputFile.dat',input_content)
    print(jobid.replace('.','p')+'_inputFile.dat')
    print("jobid: "+jobid.replace('.','p'))
    
    #Executing DISORT
    os.system('./disort_driver.sh Kokhanovsky_ve0p3re1p0COT5p0MU00p5000')
    
