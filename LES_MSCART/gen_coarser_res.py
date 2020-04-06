#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Mar 16 12:22:59 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Merging finer resolution MSCART RT results to coarser resolutions 
"""
from cpnLES_MSCARTlib import LES_case
def fctr2rsolution(fctr,c_label):
    if not(c_label.split('_',1)[0] == 'RICO'):
        #DYCOMS2/ATEXp/ATEXc 50m native resolution
        lbl = ('_%.1fkm'%(fctr*0.050)).replace('.','p')
    else:
        #RICO
        lbl = ('_%.1fkm'%(fctr*0.10)).replace('.','p')
    return lbl
def aggregate_Rad(fctr,c_label):
    LEScase=LES_case(c_label)
    des_path3D = LEScase.get_file_name(c_label+'_3D').rsplit('/',1)[0]
    des_path1D = LEScase.get_file_name(c_label+'_1D').rsplit('/',1)[0]
    RT,RTrms = LEScase.RT.movingavg_res(fctr)
    RT1D,RT1Drms = LEScase.RT1D.movingavg_res(fctr)
    LEScase.RT.MeanPRad = RT
    LEScase.RT.RMSEPRad = RTrms
    LEScase.RT.fname = LEScase.get_file_name(c_label+'_3D').rsplit('/',1)[-1].split('.',1)[0]+fctr2rsolution(fctr,c_label)
    LEScase.RT1D.MeanPRad = RT1D
    LEScase.RT1D.RMSEPRad = RT1Drms
    LEScase.RT1D.fname = LEScase.get_file_name(c_label+'_1D').rsplit('/',1)[-1].split('.',1)[0]+fctr2rsolution(fctr,c_label)
    LEScase.RT.savePOLCARThdf5(LEScase.RT.fname+'.hdf5', dpath=des_path3D+'/', pc_format='MSCART_3D_aggregated', action='coaser_resolution')
    LEScase.RT1D.savePOLCARThdf5(LEScase.RT1D.fname+'.hdf5', dpath=des_path1D+'/', pc_format='MSCART_1D_aggregated', action='coaser_resolution')
if __name__ == '__main__':
    
    #DYCOMS 50m native resolution
    #fctr = 2 # convolution domain factor
    #RICO factors [5,10,20,50]
    for fctr in [2,10,20,100]:
        for sza in ['160']:
            for band in ['0p860']:
                c_label = 'ATEXp_'+sza+'_b'+band
                print(c_label)
                aggregate_Rad(fctr,c_label)
