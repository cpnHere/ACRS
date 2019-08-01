#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Sep 14 10:01:19 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Bispectra retrieval code (based on Franks sbr_les.pro)
Note that LUT lines have to be monotonic.
"""
import numpy as np

#VNIR_lut: shorter wavelength(0.865) non-absorbing band
#SWIR_lut: longer wavelength(2.13) absorbing band

def interp(x,xp,yp):
    if xp[0]>xp[-1]:
        xp=np.flipud(xp)
        yp=np.flipud(yp)
    y=np.interp(x,xp,yp)
    
    return y 

def retrieve_NJK(VNIR,SWIR,VNIR_lut,SWIR_lut,reff_lut,tau_lut):

    #flag=2 : Reflectances are out from the LUT space
    #flag=3 : re and tau guesses are out from the LUT space
    #flag=4 : re and tau retrievals are out from the LUT space
    
    reff_lut=reff_lut.reshape(1,reff_lut.size)
    tau_lut=tau_lut.reshape(1,tau_lut.size)
    #If we are outside the LUT
    if VNIR < VNIR_lut.min() or VNIR > VNIR_lut.max() or  SWIR<SWIR_lut.min() or SWIR>SWIR_lut.max():
        flag=2
        tau_guess=np.nan
        reff_guess=np.nan
    else:
        tau_matrix=np.repeat(tau_lut,reff_lut.size,0)
        reff_matrix=np.repeat(reff_lut,tau_lut.size,0).T
    
        #cost fnction
        cost_function=(VNIR_lut-VNIR)**2/(VNIR**2)+(SWIR_lut-SWIR)**2/(SWIR**2)
        p0=np.where(cost_function==cost_function.min())
        tau_guess=tau_matrix[p0]
        reff_guess=reff_matrix[p0]
        
        for i in np.arange(0,184):
#            print('re_guess:'+str(reff_guess))
            dummy_tau=tau_guess;dummy_reff=reff_guess
            int_refl=np.squeeze(np.zeros_like(tau_lut))
            for j in np.arange(0,tau_lut.size):
                int_refl[j]=interp(reff_guess,np.squeeze(reff_lut),VNIR_lut[:,j])
            tau_guess=interp(VNIR,int_refl,np.squeeze(tau_lut))
            
            int_refl=np.squeeze(np.zeros_like(reff_lut))
            for k in np.arange(0,reff_lut.size):
                int_refl[k]=interp(tau_guess,np.squeeze(tau_lut),SWIR_lut[k,:])
            reff_guess=interp(SWIR,int_refl,np.squeeze(reff_lut))
#            if SWIR>0.25:
#                print('here')
            
            if tau_guess < np.nanmin(tau_lut) or tau_guess > np.nanmax(tau_lut) or \
                    reff_guess < np.nanmin(reff_lut) or reff_guess > np.nanmax(reff_lut):
                        flag=3 ; break
            if abs(dummy_tau-tau_guess) < 1e-7 and abs(dummy_reff-reff_guess) < 1e-7:
                        break
                    
    ret_tau=tau_guess;ret_re=reff_guess
    if ret_tau < np.nanmin(tau_lut) or ret_tau > np.nanmax(tau_lut) or \
        ret_re < np.nanmin(reff_lut) or ret_re > np.nanmax(reff_lut):
            ret_tau=np.nan
            ret_re=np.nan
            flag=4
    else:
        flag=1

    return ret_re,ret_tau,flag
