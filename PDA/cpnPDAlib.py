#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Mar  1 09:07:12 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Use PDA_Util Zhibo's script
"""
import numpy as np
import PDA_Util

class pdaRad(object):
    def __init__(self,rspname):
        '''
        PDA I/O data object
        
        '''
        self.rspname=rspname
    
    def setupPDA(self,par,wvl,atmos,surf,aer,ss_objects):
        '''
        par:PDA_Util.PDA_Parameters
        wvl:PDA_Util.PDA_Wavelength
        atmos:PDA_Util.PDA_Atmos
        surf:PDA_Util.PDA_Surface
        aer:PDA_Util.PDA_Aerosol
        ss_objects: A list of plib.PDA_SS objects
            [ss_aero,ss_cld] Should follow the same order as SS_FILE of aer (PDA_Aerosol object).
        '''
        self.PDA_Parameters=par
        self.PDA_Wavelength=wvl
        self.PDA_Atmos=atmos=atmos
        self.PDA_Surface=surf
        self.PDA_Aerosol=aer
        self.ss_objects=ss_objects
    
    def readrsp(self,rspname=None,rsp_path=None):
        if not(rspname==None):
            self.rspname=rspname
        if not(rsp_path==None):
            self.rsp_path=rsp_path
            maxview,phi0,xmu0,nview,thetav,\
                rv11,rv21,rv31, rsrf11,rsrf21,rsrf31=PDA_Util.read_rsp_ref(rsp_path+self.rspname+'.rsp')
        else:
            maxview,phi0,xmu0,nview,thetav,\
                rv11,rv21,rv31, rsrf11,rsrf21,rsrf31=PDA_Util.read_rsp_ref(self.rspname+'.rsp')
        self.I,self.Q,self.U=rv11,rv21,rv31
        self.SAA=phi0
        self.SZA=np.round(np.rad2deg(np.arccos(xmu0)))
        self.VZA=thetav
        self.RAA=phi0
    def find_DolP(self,iqu=None):
        '''
        Calculating degree of linear polarization
        '''
        self.lP=np.sqrt(self.Q**2+self.U**2)
        self.DolP=self.lP/self.I
        if not(iqu==None):
            lP=np.sqrt(iqu[1]**2+iqu[2]**2)
            DolP=lP/iqu[0]
            return DolP,lP
    def get_maxmin(self,):
        self.Imx=self.I.max()
        self.Qmx=self.Q.max()
        self.Umx=self.U.max()
        self.Imn=self.I.min()
        self.Qmn=self.Q.min()
        self.Umn=self.U.min()
    def get_ticks(self,nticks,iqu):
        if iqu=='i':
            ticks = np.linspace(self.Imn,self.Imx,nticks)
        elif iqu=='q':
            ticks = np.linspace(self.Qmn,self.Qmx,nticks)
        elif iqu=='u':
            ticks = np.linspace(self.Umn,self.Umx,nticks)
        
        return ticks

def setup_Kky_SS_files(wvl):
    '''
    wvl:PDA_Util.PDA_Wavelength 
    Read Kokhanovsky Phase matrix from a text file and return PDA_Util.PDA_SS
    objects for aerosols(ss_aero) and clouds (ss_cld)
    '''
    fl=open('Kokhanovsky_aero/Kokhanovsky_benchmark_aerosol.txt','r')
    data=fl.readlines()
    header=data[0:6]
    values=np.loadtxt(data[6:])
    fl.close()
    ang=values[:,0]
    P11=values[:,1]
    F22=values[:,2]#F22=P22/P11
    F33=values[:,3]#F33=P33/P11
    F44=values[:,4]#F44=P44/P11
    F12=values[:,5]#F12=P12/P11
    F34=values[:,6]#F34=P34/P11
    ss_aero=PDA_Util.PDA_SS(wvl.LAM[0],1.0,\
                    ang,\
                    P11,\
                    F22,\
                    F33,\
                    F44,\
                    F12,\
                    F34)
    
    
    fl=open('Kokhanovsky_cld/Kokhanovsky_benchmark_cloud.txt','r')
    data=fl.readlines()
    header=data[0:6]
    values=np.loadtxt(data[6:])
    fl.close()
    ang=values[:,0]
    P11=values[:,1]
    F22=values[:,2]#F22=P22/P11
    F33=values[:,3]#F33=P33/P11
    F44=values[:,4]#F44=P44/P11
    F12=values[:,5]#F12=P12/P11
    F34=values[:,6]#F34=P34/P11
    ss_cld=PDA_Util.PDA_SS(wvl.LAM[0],1.0,\
                    ang,\
                    P11,\
                    F22,\
                    F33,\
                    F44,\
                    F12,\
                    F34)
    
    return ss_aero,ss_cld

class Rad(object):
    def __init__(self,I,Q,U):
        '''
        I Q U Radiance object
        '''
        self.I=I
        self.Q=Q
        self.U=U
        self.get_maxmin()
    def get_maxmin(self,):
        self.Imx=self.I.max()
        self.Qmx=self.Q.max()
        self.Umx=self.U.max()
        self.Imn=self.I.min()
        self.Qmn=self.Q.min()
        self.Umn=self.U.min()
    def get_ticks(self,nticks,iqu):
        if iqu=='i':
            ticks = np.linspace(self.Imn,self.Imx,nticks)
        elif iqu=='q':
            ticks = np.linspace(self.Qmn,self.Qmx,nticks)
        elif iqu=='u':
            ticks = np.linspace(self.Umn,self.Umx,nticks)
        
        return ticks
