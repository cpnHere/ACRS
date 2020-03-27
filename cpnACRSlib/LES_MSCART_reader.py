#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Fri May 11 06:58:38 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
"""
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap

from cpnLES_MSCARTlib import LES_case, POLCARTdset, LES_field
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cpnCommonlib as cpn

def savefig(fig,fig_ttl):
    cpn.savefig(fig,fig_ttl,'figures/ATEXp/')
    
def iqu_cb(fig,ctf,ax,ticks=None,orientation='horizontal',label='label',pad=0.2):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=pad)
    if ticks is None:    
        fig.colorbar(ctf, cax=cax,orientation=orientation,label=label)
    else:
        fig.colorbar(ctf, cax=cax,ticks=ticks,orientation=orientation,label=label)
def get_rt_vc(vci):
    '''
    For RT radiance 2D histogram plots
    vci: Give index (ix) as the following table
        Contourf color bars and v (VC)
        +++++++++++++++++++++++++++++++++++++++
        DYCOMS2
        +++++++++++++++++++++++++++++++++++++++
        ---------------------------------------
        band   |  SZA  |  VZA  |  ix  |  SAA  |
        ---------------------------------------
        0p860  |  140  |  000  |   1  |  030  |
        2p13   |  140  |  000  |   2  |  030  |
        2p13   |  120  |  000  |   4  |  030  |
        ---------------------------------------
        0p860  |  120  |  000  |   3  |  000  |
        0p860  |  140  |  000  |   6  |  000  |
        0p860  |  160  |  000  |  19  |  000  |
        2p13   |  120  |  000  |  20  |  000  |
        2p13   |  140  |  000  |   5  |  000  |
        2p13   |  160  |  000  |  21  |  000  |
        
        +++++++++++++++++++++++++++++++++++++++
        RICO
        +++++++++++++++++++++++++++++++++++++++
        ---------------------------------------
        0p860  |  120  |  000  |  23  |  000  |
        0p860  |  140  |  000  |   9  |  000  |
        0p860  |  160  |  000  |  23  |  000  |
        2p13   |  120  |  000  |  24  |  000  |
        2p13   |  140  |  000  |  25  |  000  |
        2p13   |  160  |  000  |  26  |  000  |
        
        +++++++++++++++++++++++++++++++++++++++
        ATEXc
        +++++++++++++++++++++++++++++++++++++++
        ---------------------------------------
        0p860  |  120  |  000  |   7  |  000  |
        0p860  |  140  |  000  |  10  |  000  |
        0p860  |  160  |  000  |  12  |  000  |
        2p13   |  120  |  000  |   8  |  000  |
        2p13   |  140  |  000  |  11  |  000  |
        2p13   |  160  |  000  |  22  |  000  |
        
        +++++++++++++++++++++++++++++++++++++++
        ATEXp
        +++++++++++++++++++++++++++++++++++++++
        ---------------------------------------
        0p860  |  120  |  000  |   13  |  000  |
        0p860  |  140  |  000  |   14  |  000  |
        0p860  |  160  |  000  |   15  |  000  |*not finalized
        2p13   |  120  |  000  |   16  |  000  |
        2p13   |  140  |  000  |   17  |  000  |
        2p13   |  160  |  000  |   18  |  000  |    
    '''
    VC={1:{'vIR':np.linspace(0   ,1,50) ,'vQR':np.linspace(-.05,0,50)    ,'vUR':np.linspace(0.0,0.1,50) ,\
           'vIe':np.linspace(0.0,2.0,20),'vQe':np.linspace(0.0,1.00,20)  ,'vUe':np.linspace(0.0,1.0,20),\
           'cIR':np.arange(0,1.1,0.25)  ,'cQR':np.arange(-.05,0.011,0.02),'cUR':np.arange(0.0,0.11,0.05),\
           'cIe':np.arange(0.0,2.1,0.5) ,'cQe':np.arange(0.0,1.1,0.5)    ,'cUe':np.arange(0.0,1.1,0.5)},
        2:{'vIR':np.linspace(0.1 ,0.3,50) ,'vQR':np.linspace(-.02,0,50)    ,'vUR':np.linspace(0.02,0.03,50) ,\
           'vIe':np.linspace(0.75,2.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.1,0.31,0.1)  ,'cQR':np.arange(-.02,0.001,0.01),'cUR':np.arange(0.02,0.031,0.005),\
           'cIe':np.arange(0.75,2.1,0.25) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,2.1,0.5)},\
        3:{'vIR':np.linspace(0.2 ,0.8,50) ,'vQR':np.linspace(-.05,0,50)    ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.05,.001,0.025),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        4:{'vIR':np.linspace(0.1 ,0.3,50) ,'vQR':np.linspace(-.005,0,50)   ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.1,0.31,0.1)  ,'cQR':np.arange(-.005,.001,.002),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,2.1,0.5)},\
        5:{'vIR':np.linspace(0.2 ,0.4,50) ,'vQR':np.linspace(-.05,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)  ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.2,0.41,0.1)  ,'cQR':np.arange(-.05,.01,.02),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)   ,'cUe':np.arange(0,2.1,0.5)},\
        6:{'vIR':np.linspace(0.2 ,0.82,50) ,'vQR':np.linspace(-.1,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)  ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.2,0.82,0.2)  ,'cQR':np.arange(-.1,.01,.05),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)   ,'cUe':np.arange(0,2.1,0.5)},\
        7:{'vIR':np.linspace(0.2 ,0.8,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        8:{'vIR':np.linspace(0.0 ,0.3,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,0.31,0.1)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        9:{'vIR':np.linspace(0.0 ,1.0,50) ,'vQR':np.linspace(-.1,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
           'vIe':np.linspace(0.0,3.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,1.1,0.5)  ,'cQR':np.arange(-.1,0.01,0.05),'cUR':np.arange(-0.01,0.011,0.010),\
           'cIe':np.arange(0,3.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,2.5)},\
        10:{'vIR':np.linspace(0.0 ,1.0,50) ,'vQR':np.linspace(-.1,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
            'vIe':np.linspace(0.0,3.00,20),'vQe':np.linspace(0,0.5,20)     ,'vUe':np.linspace(0.0,5,20),\
            'cIR':np.arange(0.0,1.01,0.5)  ,'cQR':np.arange(-.1,0.01,0.05),'cUR':np.arange(-0.01,0.011,0.010),\
            'cIe':np.arange(0,3.1,1) ,'cQe':np.arange(0,0.51,0.25)      ,'cUe':np.arange(0,5.1,2.5)},\
        11:{'vIR':np.linspace(0.0 ,0.5,50) ,'vQR':np.linspace(-.05,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
            'vIe':np.linspace(0.0,3.00,20),'vQe':np.linspace(0,0.5,20)     ,'vUe':np.linspace(0.0,5,20),\
            'cIR':np.arange(0.0,5.51,0.25)  ,'cQR':np.arange(-.05,0.01,0.025),'cUR':np.arange(-0.01,0.011,0.010),\
            'cIe':np.arange(0,3.1,1) ,'cQe':np.arange(0,0.51,0.25)      ,'cUe':np.arange(0,5.1,2.5)},\
        12:{'vIR':np.linspace(0.2 ,0.8,50) ,'vQR':np.linspace(-.005,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.2,0.81,0.2)  ,'cQR':np.arange(-.005,.001,0.0025),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        13:{'vIR':np.linspace(0.0 ,1.0,50) ,'vQR':np.linspace(-.02,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,1.1,0.5)  ,'cQR':np.arange(-.02,.01,0.01),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        14:{'vIR':np.linspace(0.0 ,1.0,50) ,'vQR':np.linspace(-.05,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,1.1,0.5)  ,'cQR':np.arange(-.05,.01,0.025),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        15:{'vIR':np.linspace(0.0 ,1.0,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,1.1,0.5)  ,'cQR':np.arange(-.01,.01,0.005),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        16:{'vIR':np.linspace(0.0 ,0.6,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,0.61,0.2)  ,'cQR':np.arange(-.01,.01,0.005),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        17:{'vIR':np.linspace(0.0 ,0.6,50) ,'vQR':np.linspace(-.1,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,0.61,0.2)  ,'cQR':np.arange(-.1,.01,0.05),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)},\
        18:{'vIR':np.linspace(0.0 ,0.6,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
           'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,10.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,0.61,0.2)  ,'cQR':np.arange(-.01,.01,0.005),'cUR':np.arange(0.0,0.021,0.010),\
           'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,10.1,2)      ,'cUe':np.arange(0,5.1,1)},\
        19:{'vIR':np.linspace(0.2 ,0.82,50) ,'vQR':np.linspace(-.01,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)  ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.2,0.82,0.2)  ,'cQR':np.arange(-.01,.01,.005),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)   ,'cUe':np.arange(0,2.1,0.5)},\
        20:{'vIR':np.linspace(0.2 ,0.3,50) ,'vQR':np.linspace(-.01,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,1.0,20)  ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.2,0.31,0.05)  ,'cQR':np.arange(-.01,.01,.005),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,1.1,0.5)   ,'cUe':np.arange(0,2.1,0.5)},\
        21:{'vIR':np.linspace(0.2 ,0.3,50) ,'vQR':np.linspace(-.01,0,50) ,'vUR':np.linspace(-0.05,0.05,50) ,\
           'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,10.0,20)  ,'vUe':np.linspace(0.0,2,20),\
           'cIR':np.arange(0.2,0.31,0.05)  ,'cQR':np.arange(-.01,.01,.005),'cUR':np.arange(-0.05,0.051,0.05),\
           'cIe':np.arange(0.0,1.1,0.5) ,'cQe':np.arange(0,10.1,2)   ,'cUe':np.arange(0,2.1,0.5)},\
        22:{'vIR':np.linspace(0.0 ,0.3,50) ,'vQR':np.linspace(-.01,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
            'vIe':np.linspace(0.0,1.00,20),'vQe':np.linspace(0,10.0,20)     ,'vUe':np.linspace(0.0,5,20),\
            'cIR':np.arange(0.0,0.31,0.1)  ,'cQR':np.arange(-.01,0.01,0.005),'cUR':np.arange(-0.01,0.011,0.010),\
            'cIe':np.arange(0,1.1,0.5) ,'cQe':np.arange(0,10.1,2)      ,'cUe':np.arange(0,5.1,2.5)},\
        23:{'vIR':np.linspace(0.0 ,1.0,50) ,'vQR':np.linspace(-.01,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
           'vIe':np.linspace(0.0,3.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,1.1,0.5)  ,'cQR':np.arange(-.01,0.01,0.005),'cUR':np.arange(-0.01,0.011,0.010),\
           'cIe':np.arange(0,3.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,2.5)},\
        24:{'vIR':np.linspace(0.0 ,0.5,50) ,'vQR':np.linspace(-.01,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
           'vIe':np.linspace(0.0,3.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,0.51,0.25)  ,'cQR':np.arange(-.01,0.01,0.005),'cUR':np.arange(-0.01,0.011,0.010),\
           'cIe':np.arange(0,3.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,2.5)},\
        25:{'vIR':np.linspace(0.0 ,0.5,50) ,'vQR':np.linspace(-.1,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
           'vIe':np.linspace(0.0,3.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,0.51,0.25)  ,'cQR':np.arange(-.1,0.01,0.05),'cUR':np.arange(-0.01,0.011,0.010),\
           'cIe':np.arange(0,3.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,2.5)},\
        26:{'vIR':np.linspace(0.0 ,0.3,50) ,'vQR':np.linspace(-.01,0.0,50)    ,'vUR':np.linspace(-0.01,0.01,50) ,\
           'vIe':np.linspace(0.0,3.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
           'cIR':np.arange(0.0,0.31,0.1)  ,'cQR':np.arange(-.01,0.01,0.005),'cUR':np.arange(-0.01,0.011,0.010),\
           'cIe':np.arange(0,3.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,2.5)}
        }      
    return VC[vci]
def iqu_DYCOMS2(case,VZA=0,SZA=140,vci=None,vc=None,RTdim='3D',case_name='DYCOMS2'): 
    '''
    vc: dictionary to define v and colorbars for contoursf
        vc={'vIR':np.linspace(0   ,1,50) ,'vQR':np.linspace(-.05,0,50)    ,'vUR':np.linspace(0.0,0.1,50) ,\
           'vIe':np.linspace(0.0,2.0,20),'vQe':np.linspace(0.0,1.00,20)  ,'vUe':np.linspace(0.0,1.0,20),\
           'cIR':np.arange(0,1.1,0.25)  ,'cQR':np.arange(-.05,0.011,0.02),'cUR':np.arange(0.0,0.11,0.05),\
           'cIe':np.arange(0.0,2.1,0.5) ,'cQe':np.arange(0.0,1.1,0.5)    ,'cUe':np.arange(0.0,1.1,0.5)},
        
    RTdim: '1D' or '3D' RT transfer string
    '''

    if VZA==0:
        VZAi=61
    if vc==None and vci==None:
        print('Give either vc or vci')
    elif not(vci==None):
        vc=get_rt_vc(vci)
    else:
        print('vc given explicitly')
    fig1,ax1=plt.subplots(2,3,figsize=(8,6),subplot_kw={'aspect':'equal'})

    if RTdim=='3D':
        ctfI,ctfQ,ctfU,cteI,cteQ,cteU,fig1_ttl=iqu_LEScase_3D(fig1,case,ax1,VZAi,vc,case_name=case_name)
    elif RTdim=='1D':
        ctfI,ctfQ,ctfU,cteI,cteQ,cteU,fig1_ttl=iqu_LEScase_1D(fig1,case,ax1,VZAi,vc,case_name=case_name)
    ax1[0,1].tick_params(labelleft=False)
    ax1[0,2].tick_params(labelleft=False)
    ax1[1,1].tick_params(labelleft=False)
    ax1[1,2].tick_params(labelleft=False)
    
    '''
    cpn.add_cb(fig1,ctfI,ax1[0,0],orientation='vertical',ticks=vc['cIR'],label='R_I')
    cpn.add_cb(fig1,ctfQ,ax1[0,1],orientation='vertical',ticks=vc['cQR'],label='R_Q')
    cpn.add_cb(fig1,ctfU,ax1[0,2],orientation='vertical',ticks=vc['cUR'],label='R_U')
    cpn.add_cb(fig1,cteI,ax1[1,0],orientation='vertical',ticks=vc['cIe'],label='RMS%')
    cpn.add_cb(fig1,cteQ,ax1[1,1],orientation='vertical',ticks=vc['cQe'],label='RMS%')
    cpn.add_cb(fig1,cteU,ax1[1,2],orientation='vertical',ticks=vc['cUe'],label='RMS%')
    '''
    
    iqu_cb(fig1,ctfI,ax1[0,0],ticks=vc['cIR'],label='$R_I$',pad=0.4)
    iqu_cb(fig1,ctfQ,ax1[0,1],ticks=vc['cQR'],label='$R_Q$',pad=0.4)
    iqu_cb(fig1,ctfU,ax1[0,2],ticks=vc['cUR'],label='$R_U$',pad=0.4)
    iqu_cb(fig1,cteI,ax1[1,0],ticks=vc['cIe'],label='RMS%',pad=0.4)
    iqu_cb(fig1,cteQ,ax1[1,1],ticks=vc['cQe'],label='RMS%',pad=0.4)
    iqu_cb(fig1,cteU,ax1[1,2],ticks=vc['cUe'],label='RMS%',pad=0.4)
    fig1.suptitle("\n".join(wrap(fig1_ttl,60)),size=10)   
    
    return fig1,fig1_ttl,vc
    
def iqu_LEScase_3D(fig1,case,ax1,VZAi,vc,case_name="DYCOMS2"):
    cmap=plt.cm.jet
    cmap.set_bad(color='gray')
    fig1_ttl=case.RT.fname.split('.',1)[0]+'_'+case_name+'_IQU_top_RMSE_bot'    
    ctfI=ax1[0,0].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,0],vc['vIR'],cmap=cmap,extend='both')
    ctfQ=ax1[0,1].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,1],vc['vQR'],cmap=cmap,extend='both')
    ctfU=ax1[0,2].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,2],vc['vUR'],cmap=cmap,extend='both')
    cteI=ax1[1,0].contourf(case.xcens,case.ycens,case.RT.RMSEPRad[VZAi,:,:,0]/case.RT.MeanPRad[VZAi,:,:,0]*100,vc['vIe'],cmap=cmap,extend='both')
    cteQ=ax1[1,1].contourf(case.xcens,case.ycens,case.RT.RMSEPRad[VZAi,:,:,1]/case.RT.MeanPRad[VZAi,:,:,1]*100,vc['vQe'],cmap=cmap,extend='both')
    cteU=ax1[1,2].contourf(case.xcens,case.ycens,case.RT.RMSEPRad[VZAi,:,:,2]/case.RT.MeanPRad[VZAi,:,:,2]*100,vc['vUe'],cmap=cmap,extend='both')
    return ctfI,ctfQ,ctfU,cteI,cteQ,cteU,fig1_ttl
def iqu_LEScase_1D(fig1,case,ax1,VZAi,vc,case_name="DYCOMS2"):    
    cmap=plt.cm.jet
    cmap.set_bad(color='gray')
    fig1_ttl=case.RT1D.fname.split('.',1)[0]+'_'+case_name+'_IQU_top_RMSE_bot'
    ctfI=ax1[0,0].contourf(case.xcens,case.ycens,case.RT1D.MeanPRad[VZAi,:,:,0],vc['vIR'],cmap=cmap,extend='both')
    ctfQ=ax1[0,1].contourf(case.xcens,case.ycens,case.RT1D.MeanPRad[VZAi,:,:,1],vc['vQR'],cmap=cmap,extend='both')
    ctfU=ax1[0,2].contourf(case.xcens,case.ycens,case.RT1D.MeanPRad[VZAi,:,:,2],vc['vUR'],cmap=cmap,extend='both')
    cteI=ax1[1,0].contourf(case.xcens,case.ycens,(case.RT1D.RMSEPRad[VZAi,:,:,0]/case.RT1D.MeanPRad[VZAi,:,:,0]*100),vc['vIe'],cmap=cmap,extend='both')
    cteQ=ax1[1,1].contourf(case.xcens,case.ycens,(case.RT1D.RMSEPRad[VZAi,:,:,1]/case.RT1D.MeanPRad[VZAi,:,:,1]*100),vc['vQe'],cmap=cmap,extend='both')
    cteU=ax1[1,2].contourf(case.xcens,case.ycens,(case.RT1D.RMSEPRad[VZAi,:,:,2]/case.RT1D.MeanPRad[VZAi,:,:,2]*100),vc['vUe'],cmap=cmap,extend='both')
    return ctfI,ctfQ,ctfU,cteI,cteQ,cteU,fig1_ttl
def DYCOMS2_3D_1D_bias(case,VZA,SZA,vci):
    '''
    DYCOMS2 3D and 1D biases plot
    Uses iqu_DYCOMS2()
    '''
    cmap=plt.cm.jet
    cmap.set_bad(color='gray')
    cmap0=plt.cm.RdBu_r
    cmap0.set_bad(color='gray')
    if VZA==0:
        VZAi=61
    _,_,vc=iqu_DYCOMS2(case,VZA,SZA,vci=vci)
    RTbias=POLCARTdset('LES','')
    RTbias.MeanPRad=case.RT.MeanPRad-case.RT1D.MeanPRad.swapaxes(1,2)

    fig7,ax7=plt.subplots(3,2,figsize=(5,7),subplot_kw={'aspect':'equal'})
    fig7_ttl=case.RT.fname.split('.',1)[0]+'_DYCOM2_IQU_3D_1D_biases'
    fig7.suptitle("\n".join(wrap(fig7_ttl,50)),size=10)
    
    ctfI=ax7[0,0].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,0]  ,vc['vIR'],cmap=cmap,extend='both');ax7[0,0].set_title('I_3D',size=10)
    ctfQ=ax7[0,1].contourf(case.xcens,case.ycens,case.RT.MeanPRad[VZAi,:,:,1]  ,vc['vQR'],cmap=cmap,extend='both');ax7[0,1].set_title('Q_3D',size=10)
    ctI1=ax7[1,0].contourf(case.xcens,case.ycens,(case.RT1D.MeanPRad[VZAi,:,:,0]).T,vc['vIR'],cmap=cmap,extend='both');ax7[1,0].set_title('I_1D',size=10)
    ctQ1=ax7[1,1].contourf(case.xcens,case.ycens,(case.RT1D.MeanPRad[VZAi,:,:,1]).T,vc['vQR'],cmap=cmap,extend='both');ax7[1,1].set_title('Q_1D',size=10)
    
    ctfb1=ax7[2,0].contourf(case.xcens,case.ycens,RTbias.MeanPRad[VZAi,:,:,0],np.linspace(-.1,.1,50)  ,cmap=cmap0,extend='both');ax7[2,0].set_title('I Bias',size=10)
    ctfb2=ax7[2,1].contourf(case.xcens,case.ycens,RTbias.MeanPRad[VZAi,:,:,1],np.linspace(-.01,.01,50),cmap=cmap0,extend='both');ax7[2,1].set_title('Q Bias',size=10)
    
    ax7[0,0].set_ylabel('km',size=10)
    ax7[1,0].set_ylabel('km',size=10)
    ax7[2,0].set_ylabel('km',size=10)
    
    ax7[0,1].tick_params(labelleft=False)
    ax7[1,1].tick_params(labelleft=False)
    ax7[2,1].tick_params(labelleft=False)
    
    iqu_cb(fig7,ctfI,ax7[0,0],ticks=vc['cIR'],label='')
    iqu_cb(fig7,ctfQ,ax7[0,1],ticks=vc['cQR'],label='')
    iqu_cb(fig7,ctI1,ax7[1,0],ticks=vc['cIR'],label='')
    iqu_cb(fig7,ctQ1,ax7[1,1],ticks=vc['cQR'],label='')
    iqu_cb(fig7,ctfb1,ax7[2,0],ticks=np.arange(-.1,.11,.1) ,label='')
    iqu_cb(fig7,ctfb2,ax7[2,1],ticks=np.arange(-.01,.011,.01),label='')
    cpn.sub_labels(ax7)
    return fig7,fig7_ttl

def old_vs_new_SZA140(dynyx,dyo):    
    #Comparing old and new runs
    fig1,ax1=plt.subplots(3,3,figsize=(8,8),subplot_kw={'aspect':'equal'})
    fig1_ttl=dynyx.RT.fname.split('.',1)[0]+'_DYCOM2_IQU_old_vs_new_3D'
    fig1.suptitle("\n".join(wrap(fig1_ttl,60)))
    
    ctfI=ax1[0,0].contourf(dynyx.xcens,dynyx.ycens,dynyx.RT.MeanPRad[61,:,:,0],np.linspace(0   ,1,50) ,cmap=plt.cm.jet)
    ctfQ=ax1[0,1].contourf(dynyx.xcens,dynyx.ycens,dynyx.RT.MeanPRad[61,:,:,1],np.linspace(-.05,0,50) ,cmap=plt.cm.jet)
    ctfU=ax1[0,2].contourf(dynyx.xcens,dynyx.ycens,dynyx.RT.MeanPRad[61,:,:,2],np.linspace(0.0,0.1,50),cmap=plt.cm.jet)
    ax1[1,0].contourf(dyo.xcens,dynyx.ycens  ,dyo.RT.MeanPRad[61,:,:,0],np.linspace(0   ,1,50) ,cmap=plt.cm.jet)
    ax1[1,1].contourf(dyo.xcens,dynyx.ycens  ,dyo.RT.MeanPRad[61,:,:,1]  ,np.linspace(-.05,0,50) ,cmap=plt.cm.jet)
    ax1[1,2].contourf(dyo.xcens,dynyx.ycens  ,dyo.RT.MeanPRad[61,:,:,2]  ,np.linspace(0.0,0.1,50),cmap=plt.cm.jet)
    
    b1,b2,b3=dynyx.RT.MeanPRad[61,:,:,0]-dyo.RT.MeanPRad[61,:,:,0],\
             dynyx.RT.MeanPRad[61,:,:,1]-dyo.RT.MeanPRad[61,:,:,1],\
             dynyx.RT.MeanPRad[61,:,:,2]-dyo.RT.MeanPRad[61,:,:,2]
    ctfb1=ax1[2,0].contourf(dyo.xcens,dynyx.ycens,b1,np.linspace(-.003,.003,50)  ,cmap=plt.cm.RdBu_r,extend='both')
    ctfb2=ax1[2,1].contourf(dyo.xcens,dynyx.ycens,b2,np.linspace(-.002,.002,50),cmap=plt.cm.RdBu_r,extend='both')
    ctfb3=ax1[2,2].contourf(dyo.xcens,dynyx.ycens,b3,np.linspace(-.002,.002,50),cmap=plt.cm.RdBu_r,extend='both')
    
    ax1[0,1].set_title('New',size=10)
    ax1[1,1].set_title('Old',size=10)
    ax1[2,1].set_title('Bias',size=10)
    
    ax1[0,0].set_ylabel('km',size=10)
    ax1[1,0].set_ylabel('km',size=10)
    ax1[2,0].set_ylabel('km',size=10)
    
    ax1[0,1].tick_params(labelleft=False)
    ax1[0,2].tick_params(labelleft=False)
    iqu_cb(fig1,ctfI,ax1[0,0],ticks=np.arange(0,1.1,0.25)     ,label='')
    iqu_cb(fig1,ctfQ,ax1[0,1],ticks=np.arange(-.05,0.011,0.02),label='')
    iqu_cb(fig1,ctfU,ax1[0,2],ticks=np.arange(0.0,0.11,0.05)  ,label='')
    iqu_cb(fig1,ctfI,ax1[1,0],ticks=np.arange(0,1.1,0.25)     ,label='')
    iqu_cb(fig1,ctfQ,ax1[1,1],ticks=np.arange(-.05,0.011,0.02),label='')
    iqu_cb(fig1,ctfU,ax1[1,2],ticks=np.arange(0.0,0.11,0.05)  ,label='')
    iqu_cb(fig1,ctfb1,ax1[2,0],np.arange(-.003,.0031,0.003),label='R_I')
    iqu_cb(fig1,ctfb2,ax1[2,1],np.arange(-.002,.0021,0.002),label='R_Q')
    iqu_cb(fig1,ctfb3,ax1[2,2],np.arange(-.002,.0021,0.002),label='R_U')
    
    fig1.show() 
    return fig1,fig1_ttl
def get_binTicks(data,dc):
    bmx=np.round(data.max(),dc)
    bmn=np.round(data.min(),dc)
    return bmn,bmx


if __name__=='__main__':
    cpn.setup_figures(plt)
    '''
    DYCOMS-II (DYCOMS2) 1D 3D RT examples
    ===========================================================================
    les_new_path='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/'
    TSc0p860_120=LES_case('ToweringSc_dharma_001800_X63_83Y35_55_framed_b0p865_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5',les_new_path,band='0p865')
    
    DYC0p860_140=LES_case('DYCOMS2_dharma_008036_b0p860_MSCART_SZA140_SAA000_VAA000plus_NPH2e6.hdf5',les_new_path,\
                          RT1Dname='DYCOMS2_dharma_008036_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5')
    DYC2p13_140 =LES_case('DYCOMS2_dharma_008036_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH2e6.hdf5',les_new_path,\
                            RT1Dname='DYCOMS2_dharma_008036_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5')

    DYC0p860_120=LES_case('DYCOMS2_dharma_008036_b0p860_MSCART_SZA120_SAA000_VAA000plus_NPH2e6.hdf5',les_new_path,\
                            RT1Dname='DYCOMS2_dharma_008036_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5')
    DYC2p13_120 =LES_case('DYCOMS2_dharma_008036_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH2e6.hdf5',les_new_path,\
                            RT1Dname='DYCOMS2_dharma_008036_b2p13_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5')

    DYC0p860_160=LES_case('DYCOMS2_dharma_008036_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH2e6.hdf5',les_new_path,\
                            RT1Dname='DYCOMS2_dharma_008036_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5')
    DYC2p13_160 =LES_case('DYCOMS2_dharma_008036_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH2e6.hdf5',les_new_path,\
                            RT1Dname='DYCOMS2_dharma_008036_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5')
    
    #Reading old 120 3D and 1D cases
#    dyo120   =LES_case('OP_dharma_008036_full_3_26_MSCART_SZA060_SAA030_VAA000plus_NPH2e6.hdf5'  ,'/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/','results/hdf5/',band='0p860')    
#    RT=POLCARTdset('LES','/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/')
#    RT.readPOLCARThdf5('/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART/results/hdf5/DYCOM_MSCART_1D_SZA%03d_SAA030.hdf5'%60)
    
        

# Comparing with old DYCOMS2 MSCART simulations
#    fig1,fig1_ttl=old_vs_new_SZA140(dynyx0p860_140,dyo140)
#    fig2,fig2_ttl=iqu_DYCOMS2(dyo140,0,140)
    
# 1D and 3D iqu plots
    fig3,fig3_ttl,_=iqu_DYCOMS2(DYC0p860_140,0,140,vci=6);fig3.show()
    fig4,fig4_ttl,_=iqu_DYCOMS2(DYC2p13_140 ,0,140,vci=5);fig4.show()
    fig5,fig5_ttl,_=iqu_DYCOMS2(DYC0p860_140,0,140,vci=6,RTdim='1D');fig5.show()
    fig6,fig6_ttl,_=iqu_DYCOMS2(DYC2p13_140 ,0,140,vci=5,RTdim='1D');fig6.show()
    fig7,fig7_ttl=DYCOMS2_3D_1D_bias(DYC0p860_120,0,120,3);fig7.show()
    fig8,fig8_ttl=DYCOMS2_3D_1D_bias(DYC2p13_120 ,0,120,4);fig8.show()
    
    fig9 ,fig9_ttl ,_=iqu_DYCOMS2(DYC0p860_120,0,120,vci=3);fig9.show()
    fig10,fig10_ttl,_=iqu_DYCOMS2(DYC2p13_120 ,0,120,vci=4);fig10.show()
    fig11,fig11_ttl,_=iqu_DYCOMS2(DYC0p860_120,0,120,vci=3,RTdim='1D');fig11.show()
    fig12,fig12_ttl,_=iqu_DYCOMS2(DYC2p13_120 ,0,120,vci=4,RTdim='1D');fig12.show()
    fig13 ,fig13_ttl =DYCOMS2_3D_1D_bias(DYC0p860_140,0,140,6);fig13.show()
    fig14 ,fig14_ttl =DYCOMS2_3D_1D_bias(DYC2p13_140 ,0,140,5);fig14.show()

    fig15,fig15_ttl ,_=iqu_DYCOMS2(DYC0p860_160,0,160,vci=3);fig15.show()
    fig16,fig16_ttl,_=iqu_DYCOMS2(DYC2p13_160 ,0,160,vci=4);fig16.show()
    fig17,fig17_ttl,_=iqu_DYCOMS2(DYC0p860_160,0,160,vci=3,RTdim='1D');fig17.show()
    fig18,fig18_ttl,_=iqu_DYCOMS2(DYC2p13_160 ,0,160,vci=4,RTdim='1D');fig18.show()
    fig19 ,fig19_ttl =DYCOMS2_3D_1D_bias(DYC0p860_160,0,160,3);fig19.show()
    fig20 ,fig20_ttl =DYCOMS2_3D_1D_bias(DYC2p13_160 ,0,160,4);fig20.show()
    '''
    '''
    ATEXc 1D 3D RT examples
    ===========================================================================
    '''    
    '''
    #ATEXc full cases
    base_dir='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/'    
    ATEXc0p860_120=LES_case('ATEXc_dharma_007877_b0p860_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXc2p13_120 =LES_case('ATEXc_dharma_007877_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXc_dharma_007877_b2p13_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXc0p860_140=LES_case('ATEXc_dharma_007877_b0p860_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXc2p13_140 =LES_case('ATEXc_dharma_007877_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXc_dharma_007877_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXc0p860_160=LES_case('ATEXc_dharma_007877_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXc_dharma_007877_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXc2p13_160 =LES_case('ATEXc_dharma_007877_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXc_dharma_007877_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5')
            
    vc={'vIR':np.linspace(0.0 ,0.4,50) ,'vQR':np.linspace(-.01,0,50)    ,'vUR':np.linspace(0.0,0.02,50) ,\
     'vIe':np.linspace(0.0,5.00,20),'vQe':np.linspace(0,1.0,20)     ,'vUe':np.linspace(0.0,5,20),\
     'cIR':np.arange(0.0,0.41,0.1)  ,'cQR':np.arange(-.01,.001,0.005),'cUR':np.arange(0.0,0.021,0.010),\
     'cIe':np.arange(0,5.1,1) ,'cQe':np.arange(0,1.1,0.5)      ,'cUe':np.arange(0,5.1,1)}    
    
    fig21,fig21_ttl,_=iqu_DYCOMS2(ATEXc0p860_120,0,120,vci=7,case_name='ATEXc');fig21.tight_layout(rect=[0,0,1,0.96]);fig21.show()
    fig22,fig22_ttl,_=iqu_DYCOMS2(ATEXc2p13_120 ,0,120,vci=8,case_name='ATEXc');fig22.tight_layout(rect=[0,0,1,0.96]);fig22.show()
    fig23,fig23_ttl,_=iqu_DYCOMS2(ATEXc0p860_120,0,120,vci=7,case_name='ATEXc',RTdim='1D');fig23.tight_layout(rect=[0,0,1,0.96]);fig23.show()
    fig24,fig24_ttl,_=iqu_DYCOMS2(ATEXc2p13_120 ,0,120,vci=8,case_name='ATEXc',RTdim='1D');fig24.tight_layout(rect=[0,0,1,0.96]);fig24.show()    
    fig25,fig25_ttl  =DYCOMS2_3D_1D_bias(ATEXc0p860_120,0,120,7);fig25.show()
    fig26,fig26_ttl  =DYCOMS2_3D_1D_bias(ATEXc2p13_120,0,120,8);fig26.show()
    
    fig27,fig27_ttl,_=iqu_DYCOMS2(ATEXc0p860_140,0,140,vci=7,case_name='ATEXc');fig27.tight_layout(rect=[0,0,1,0.96]);fig27.show()
    fig28,fig28_ttl,_=iqu_DYCOMS2(ATEXc2p13_140 ,0,140,vci=8,case_name='ATEXc');fig28.show()
    fig29,fig29_ttl,_=iqu_DYCOMS2(ATEXc0p860_140,0,140,vci=7,case_name='ATEXc',RTdim='1D');fig29.show()
    fig30,fig30_ttl,_=iqu_DYCOMS2(ATEXc2p13_140 ,0,140,vci=8,case_name='ATEXc',RTdim='1D');fig30.show()
    fig31,fig31_ttl  =DYCOMS2_3D_1D_bias(ATEXc0p860_140,0,140,7);fig31.show()
    fig32,fig32_ttl  =DYCOMS2_3D_1D_bias(ATEXc2p13_140,0,140,8);fig32.show()

    fig33,fig33_ttl,_=iqu_DYCOMS2(ATEXc0p860_160,0,160,vci=7,case_name='ATEXc');fig33.tight_layout(rect=[0,0,1,0.96]);fig33.show()
    fig34,fig34_ttl,_=iqu_DYCOMS2(ATEXc2p13_160 ,0,160,vci=8,case_name='ATEXc');fig34.show()
    fig35,fig35_ttl,_=iqu_DYCOMS2(ATEXc0p860_160,0,160,vci=7,case_name='ATEXc',RTdim='1D');fig35.show()
    fig36,fig36_ttl,_=iqu_DYCOMS2(ATEXc2p13_160 ,0,160,vci=8,case_name='ATEXc',RTdim='1D');fig36.show()
    fig37,fig37_ttl  =DYCOMS2_3D_1D_bias(ATEXc0p860_160,0,160,7);fig37.show()
    fig38,fig38_ttl  =DYCOMS2_3D_1D_bias(ATEXc2p13_160,0,160,8);fig38.show()
    '''
    '''
    ATEXp 1D 3D RT examples
    ===========================================================================
    '''    
    #ATEXp full cases

    base_dir='/umbc/xfs1/zzbatmos/users/charaj1/LES_MSCART_v2/ATEXp/'    

    ATEXp0p860_120=LES_case('ATEXp_dharma_013067_b0p860_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp2p13_120 =LES_case('ATEXp_dharma_013067_b2p13_MSCART_SZA120_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA120_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp0p860_140=LES_case('ATEXp_dharma_013067_b0p860_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp2p13_140 =LES_case('ATEXp_dharma_013067_b2p13_MSCART_SZA140_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA140_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp0p860_160=LES_case('ATEXp_dharma_013067_b0p860_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                          RT1Dname='ATEXp_dharma_013067_b0p860_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5')
    ATEXp2p13_160 =LES_case('ATEXp_dharma_013067_b2p13_MSCART_SZA160_SAA000_VAA000plus_NPH1e6.hdf5',base_dir,\
                            RT1Dname='ATEXp_dharma_013067_b2p13_MSCART_1D_bins_SZA160_SAA000_VAA000plus_NPH1e5.hdf5')

    fig1,fig1_ttl,_=iqu_DYCOMS2(ATEXp0p860_120,0,120,vci=7,case_name='ATEXp');fig1.tight_layout(rect=[0,0,1,0.96]);fig1.show()
    fig2,fig2_ttl,_=iqu_DYCOMS2(ATEXp2p13_120 ,0,120,vci=8,case_name='ATEXp');fig2.tight_layout(rect=[0,0,1,0.96]);fig2.show()
    fig3,fig3_ttl,_=iqu_DYCOMS2(ATEXp0p860_120,0,120,vci=7,case_name='ATEXp',RTdim='1D');fig3.tight_layout(rect=[0,0,1,0.96]);fig3.show()
    fig4,fig4_ttl,_=iqu_DYCOMS2(ATEXp2p13_120 ,0,120,vci=8,case_name='ATEXp',RTdim='1D');fig4.tight_layout(rect=[0,0,1,0.96]);fig4.show()    
    fig5,fig5_ttl  =DYCOMS2_3D_1D_bias(ATEXp0p860_120,0,120,7);fig5.show()
    fig6,fig6_ttl  =DYCOMS2_3D_1D_bias(ATEXp2p13_120,0,120,8);fig6.show()
    
    fig7,fig7_ttl,_=iqu_DYCOMS2(ATEXp0p860_140,0,140,vci=7,case_name='ATEXp');fig7.tight_layout(rect=[0,0,1,0.96]);fig7.show()
    fig8,fig8_ttl,_=iqu_DYCOMS2(ATEXp2p13_140 ,0,140,vci=8,case_name='ATEXp');fig8.show()
    fig9,fig9_ttl,_=iqu_DYCOMS2(ATEXp0p860_140,0,140,vci=7,case_name='ATEXp',RTdim='1D');fig9.show()
    fig10,fig10_ttl,_=iqu_DYCOMS2(ATEXp2p13_140 ,0,140,vci=8,case_name='ATEXp',RTdim='1D');fig10.show()
    fig11,fig11_ttl  =DYCOMS2_3D_1D_bias(ATEXp0p860_140,0,140,7);fig11.show()
    fig12,fig12_ttl  =DYCOMS2_3D_1D_bias(ATEXp2p13_140,0,140,8);fig12.show()

    fig13,fig13_ttl,_=iqu_DYCOMS2(ATEXp0p860_160,0,160,vci=7,case_name='ATEXp');fig13.tight_layout(rect=[0,0,1,0.96]);fig13.show()
    fig14,fig14_ttl,_=iqu_DYCOMS2(ATEXp2p13_160 ,0,160,vci=8,case_name='ATEXp');fig14.show()
    fig15,fig15_ttl,_=iqu_DYCOMS2(ATEXp0p860_160,0,160,vci=7,case_name='ATEXp',RTdim='1D');fig15.show()
    fig16,fig16_ttl,_=iqu_DYCOMS2(ATEXp2p13_160 ,0,160,vci=8,case_name='ATEXp',RTdim='1D');fig16.show()
    fig17,fig17_ttl  =DYCOMS2_3D_1D_bias(ATEXp0p860_160,0,160,7);fig17.show()
    fig18,fig18_ttl  =DYCOMS2_3D_1D_bias(ATEXp2p13_160,0,160,8);fig18.show()
    plt.close('all')
    savefig(fig1,fig1_ttl)
    savefig(fig2,fig2_ttl)
    savefig(fig3,fig3_ttl)
    savefig(fig4,fig4_ttl)
    savefig(fig5,fig5_ttl)
    savefig(fig6,fig6_ttl)
    savefig(fig7,fig7_ttl)
    savefig(fig8,fig8_ttl)
    savefig(fig9,fig9_ttl)
    savefig(fig10,fig10_ttl)
    savefig(fig11,fig11_ttl)
    savefig(fig12,fig12_ttl)
    savefig(fig13,fig13_ttl)
    savefig(fig14,fig14_ttl)
    savefig(fig15,fig15_ttl)
    savefig(fig16,fig16_ttl)
    savefig(fig17,fig17_ttl)
    savefig(fig18,fig18_ttl)

