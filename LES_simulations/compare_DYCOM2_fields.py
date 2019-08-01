#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu May 10 12:46:34 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Comparing new DYCOMS2 MSCART field file with the old one

"""
import numpy as np
import cpnLES_MSCARTlib as MSCARTlib
import cpnCommonlib as cpn

def check_PMs():
    figdyn,ttldyn=dyn.plot_PM('PM_new',show=False)
    figdy ,ttldy =dy.plot_PM('PM_old',show=False)
    for i in np.arange(0,25,1):
        cpn.savefig(figdyn[i],ttldyn[i],'figures/DYCOMS2_comp/')
        cpn.savefig(figdy[i],ttldy[i],'figures/DYCOMS2_comp/')
        
dy=MSCARTlib.LES_field('OP_dharma_008036_full_3_26.nc',dpath='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/')#DYCOM field file for MSCART
dy.readLES_field()

dyn=MSCARTlib.LES_field('DYCOMS2_dharma_008036_b0p860.nc')#DYCOM field file for MSCART
dyn.readLES_field()

