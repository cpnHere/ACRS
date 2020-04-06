#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Apr  6 11:25:30 2020
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To generate coaser resolution vertically-weighted results using native resolution products.
"""
from analysis_lib import vw_rets
from cpnCommonlib import movingaverage2D
import os,h5py
import numpy as np
def aggregate_vw_to_coarser(case_name,SZA,VZA):
    #define 2D moving avg factors and corresponding resolutions 
    if case_name == 'RIC':
        lb_list = ['0p5km','1p0km','2p0km','5p0km']; 
        fctr_li = [5,10,20,50]
    else:
        lb_list = ['0p1km','0p5km','1p0km','5p0km'];     
        fctr_li = [2,10,20,100]

    #Reading native VW retrievals 
    vw_case = vw_rets(case_name,SZA,VZA,a=1,b=0)
    fname = vw_case.VW.fname
    path = vw_case.VW.fpath
    co_res_name = path+fname+'_co_res'
    
    #Looking for [x,y] 2D variables
    twoD_vars = []
    obj = vw_case.VW
    variables = vars(obj)
    for key in variables.keys():
        if type(variables[key])==str:
            print('String variable:'+key)
        elif type(variables[key])==np.ndarray:
            if variables[key].ndim == 2 :
                print(key+' shape:'+str(np.shape(variables[key])))
                twoD_vars = np.append(twoD_vars,key)
    #Aggregating data
    co_res_vars = {}
    for i in range(np.size(fctr_li)):
        for key in twoD_vars:
            co_res_vars[key+'_'+lb_list[i]]=movingaverage2D(variables[key],window=fctr_li[i])
    return co_res_vars,co_res_name,obj

if __name__=='__main__':
    for case_name in ['DYC','RIC','ATp','ATc']:
        for SZA in [20,40,60]:
            print('==================================================================')
            print(case_name+' SZA%d'%SZA)
            VZA = 0
            co_res_vars,co_res_name,VW = aggregate_vw_to_coarser(case_name,SZA,VZA)  
        
            out_name = co_res_name+'.hdf5'
            if not(os.path.isfile(out_name)):
                f=h5py.File(out_name,'w')
                f.attrs['native_file']=out_name.replace('_co_res.','.')
                
                PC=f.create_dataset('x',data=VW.x)
                PC.attrs['units']='km'
                PC.attrs['long_name']='X distance'
            
                PC=f.create_dataset('y',data=VW.y)
                PC.attrs['units']='km'
                PC.attrs['long_name']='Y distance'
            
                for key in co_res_vars.keys():
                    if key.split('_',1)[0] == 'Re':
                        PC=f.create_dataset(key,data=co_res_vars[key])
                        PC.dims[0].label='x'
                        PC.dims[1].label='y'
                        PC.attrs['units']='Microns'
                        PC.attrs["long_name"]='See_native_file'
                    else:
                        PC=f.create_dataset(key,data=co_res_vars[key])
                        PC.dims[0].label='x'
                        PC.dims[1].label='y'
                        PC.attrs['units']='None'
                        PC.attrs["long_name"]='See_native_file'
            
                f.close()
                print(out_name+' SAVED!')
            else:
                print('File already exist!! '+out_name)

        
    
     