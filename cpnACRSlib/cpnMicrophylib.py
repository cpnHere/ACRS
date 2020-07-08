#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Aug 31 14:10:17 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Library to handle microphysics and single scattering properties. (Can add Mie
code stuff later)
09/29/2017:
    Added fractal cloud generation part.
12/28/2017:
    pkl_clsses() object added: samahara pkl save karaddi ewaye thiyana objects
    mama mehema ඩමී ක්ලාස් ekakata dala thiyanwa. ethakota ewa read karaddi me class eka
    ලෝඩ් wela thiyenta oona.
05/15/2018:
    generate_fractal2D() added: To generate 2D fractal cloud based on Cahalan et. al. (1994)    
"""
import h5py, os
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks_cwt

class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]

class mphicDset(object):
    '''
    #to read microphysics from Dan
    #for this file "MODIS_band_0.865_polarimetric_microphysical_LUT.h5"
    # or djm_blk_SS_MODIS_extended.h5 files
    '''
    def __init__(self,fname):
        self.fname=fname
    def readd(self,bands='single'):
        #nband='single' or 'multiple'
        f=h5py.File(self.fname)
        self.re=np.squeeze(f['re'][:])
        self.ve=np.squeeze(f['ve'][:])
        self.alb=f['al_blk_avg'][:]
        self.P11=f['P11_blk_avg'][:]
        self.P12=f['P12_blk_avg'][:]
        self.P33=f['P33_blk_avg'][:]
        self.P34=f['P34_blk_avg'][:]
        self.ang=np.squeeze(f['ang'][:])
        if bands=='multiple':
            self.band=np.squeeze(f['band_val'][:])
            self.Qe=f['Qe_blk_avg'][:]

class frac_physics(object):
    '''
    To generate 1D and 2D fracta clouds based on Cahalan et. al. (1994)
    fname: filename
    dpath: path
    dims: 1 or 2 (1D or 2D cloud)
    '''
    def __init__(self,fname,dpath=None,dims=1):
        self.dims=dims
        if dpath==None:
            self.fname=fname
        else:
            self.fname=dpath+fname
    def read(self,):
        # Me Dange fractal cloud eka kiyawanta
        f=h5py.File(self.fname,'r')
        self.re=f['re'][0]
        self.tau=f['tau'][0]
        self.ve=f['ve'][0]
        self.x=f['x'][0]
        self.y=f['y'][0]
        f.close()
    def generate_fractal(self,re,lwp,f=0.5,xdist=None,xorder=12):
        '''
        Generate 1D Fractal cloud - Cahalan et. al. (1994)
        re: re in mu
        lwp: LWP in g/m^2 , 90 for stratocumulus cloud
        f: Fractal parameter fn=fc^n
        xdist: horizontal distance of the cloud in km
        xorder: Order of the fractal cloud
        '''
        self.f=f
        self.orderX=xorder
        bc=lwp
        varReduction=0.8#c=2^{-1/3}=0.8
        var=self.f#fractal parameter fn=fc^n;
        def next_step_1d(bc,var):
            l=np.size(bc)
            bc2=np.zeros(l*2)
            sign=np.random.uniform(-1,1,size=(l))
            sign=sign/abs(sign)
            bc2[0:2*l:2]=bc+sign*var*bc
            bc2[1:2*l:2]=bc-sign*var*bc
            return bc2
            
        for i in np.arange(0,self.orderX):
            bc=next_step_1d(bc,var)
            var=var*varReduction
        self.lwp=bc
        self.tau=bc*3.0/2.0/re
        self.y=np.array([1.0])
        self.re=np.repeat(re,len(bc))
        if xdist==None:
            self.x=np.linspace(0,len(bc)/100.0,len(bc)+1)#let's say this is 0.01km resolution
        else:
            self.x=np.linspace(0,xdist,len(bc)+1)
        self.x=(self.x[0:-1]+self.x[1:])/2#assigning the middle
    def generate_fractal2D(self,re,lwp,f=0.5,xydist=None,xyorder=12):
        '''
        Generate 2D Fractal cloud. Extending same theory in generate_fractal()
        re: re in mu
        lwp: LWP in g/m^2 , 90 for stratocumulus cloud
        f: Fractal parameter fn=fc^n
        xydist: horizontal distance of the cloud in km
        xyorder: Order of the fractal cloud
        '''
        self.f=f
        self.orderXY=xyorder
        bc=lwp
        varReduction=0.8#c=2^{-1/3}=0.8
        var=self.f#fractal parameter fn=fc^n;
        def next_step_2D(bc,var):
            l=np.size(bc)
            xy=int(np.sqrt(l))
            bc2=np.zeros((xy*2,xy*2),dtype=float)
            bc3=np.zeros((xy*2,xy*2),dtype=float)
            #4 random signs
            sign1=np.random.uniform(-1,1,size=(xy,xy))
            sign1=sign1/abs(sign1)
            sign2=np.random.uniform(-1,1,size=(xy,xy))
            sign2=sign2/abs(sign2)
            sign3=np.random.uniform(-1,1,size=(xy,xy))
            sign3=sign3/abs(sign3)
            sign4=np.random.uniform(-1,1,size=(xy,xy))
            sign4=sign4/abs(sign4)
            
            sign4expand=np.zeros((xy,xy),dtype=int)
            for i in np.arange(0,xy):
                for j in np.arange(0,xy):
                    sign4expand[i,j]=sign4[i/2,j/2]

            bc2[0:2*l:2,0:2*l:2]=bc+sign1*var*bc+sign2*var*var*bc
            bc2[0:2*l:2,1:2*l:2]=bc+sign1*var*bc-sign2*var*var*bc
            bc2[1:2*l:2,0:2*l:2]=bc-sign1*var*bc+sign3*var*var*bc
            bc2[1:2*l:2,1:2*l:2]=bc-sign1*var*bc-sign3*var*var*bc

            bc3[0:2*l:2,0:2*l:2]=bc+sign1*var*bc+sign2*var*var*bc
            bc3[1:2*l:2,0:2*l:2]=bc+sign1*var*bc-sign2*var*var*bc
            bc3[0:2*l:2,1:2*l:2]=bc-sign1*var*bc+sign3*var*var*bc
            bc3[1:2*l:2,1:2*l:2]=bc-sign1*var*bc-sign3*var*var*bc

            bc2[np.where(sign4expand==-1)]=bc3[np.where(sign4expand==-1)]

            return bc2        
        for i in np.arange(0,self.orderXY):
            bc=next_step_2D(bc,var)
            var=var*varReduction
        self.lwp=bc
        self.tau=bc*3.0/2.0/re
        self.re=np.repeat(re,len(bc))
        if xydist==None:
            self.x=np.linspace(0,len(bc)/100.0,len(bc)+1)#let's say this is 0.01km resolution
        else:
            self.x=np.linspace(0,xydist,len(bc)+1)
        self.x=(self.x[0:-1]+self.x[1:])/2#assigning the middle
        self.y=self.x
    def save_frach5(self,fname=None):
        '''
        Save fractal cloud properties as an hdf5 file.
        '''
        if not(fname==None):
            self.fname=fname
        
        if os.path.isfile(self.fname):
            rd=input('File already exist. Replace?: ')
        else:
            rd='y'
        if rd=='y':
            f=h5py.File(self.fname,'w')
            PCentry=f.create_dataset('tau',data=self.tau)
            PCentry.dims[0].label='x'
            PCentry.attrs["long_name"]='Cloud_optical_thickness'
            
            PCentry=f.create_dataset('orderX',data=self.orderX)
            PCentry.attrs["long_name"]='x_dimension_order'
            
            PCentry=f.create_dataset('lwp',data=self.lwp)
            PCentry.dims[0].label='x'
            PCentry.attrs["long_name"]='liquid_water_path'
            PCentry.attrs["unit"]='g/m^2'
            
            PCentry=f.create_dataset('re',data=self.re)
            PCentry.dims[0].label='x'
            PCentry.attrs["long_name"]='effective_radius'
            PCentry.attrs["unit"]='microns'
            
            PCentry=f.create_dataset('x',data=self.x)
            PCentry.attrs["long_name"]='x_dimension'
            PCentry.attrs["unit"]='km'
            
            f.create_dataset('y',data=self.y)
            f.create_dataset('f',data=self.f)
            
            f.close()
            print(self.fname+' SAVED!')
    def read_frach5(self,fname=None):
        if not(fname==None):
            self.fname=fname
        f=h5py.File(self.fname,'r')
        self.re=f['re'][:]
        self.tau=f['tau'][:]
        self.f=f['f']
        self.orderX=f['orderX']
        self.lwp=f['lwp'][:]
        self.x=f['x'][:]
        self.y=f['y'][:]
        f.close()
        
    def plot_lwp_pdf(self,):
        '''
        Plot lwp PDF for both 1D and 2D case.
        '''
        lwp,lwpax=plt.subplots()
        if np.size(self.lwp.shape)>1:
            lwpval=self.lwp.reshape(self.lwp.shape[0]**2)
        else:
            lwpval=self.lwp
        weights = np.ones_like(lwpval)/len(lwpval)
        lwpax.hist(lwpval,bins = 10 ** np.linspace(1, 3, 50),weights=weights,histtype=u'step')
        #val,bins=np.histogram(self.lwp,bins = 10 ** np.linspace(1, 3, 50),weights=weights)
        #width = np.diff(bins)
        #center = (bins[:-1] + bins[1:]) / 2
        #lwpax.bar(center,val,align='center',width=width)
        lwpax.set_xscale('log')
        #x=np.linspace(ss.lognorm.ppf(0.01,0.954),ss.lognorm.ppf(0.99,0.954),100)
        #lwpax.plot(x,ss.lognorm.pdf(x,0.954))
        lwp.show()
        return lwp,lwpax

    
    def plot_cld(self,fig=None,ax=None):
        
        if ax==None:
            fig,ax=plt.subplots()
        ax.plot(self.x,self.re,'r',label='re')
        ax.set_ylabel(r'Effective radius ($\mu$)')
        ax2=ax.twinx()
        ax2.set_ylabel(r'Optical depth')
        ax2.plot(self.x,self.tau,'g',label=r'$\tau$')
        ax.legend(loc='best',prop={'size':10})
        ax.set_xlabel('km')
        fig.show()
        
        return fig

def find_bows(ang,P11):
    '''
    Find peaks in P11
    #ang: scattering angle in degrees
    #P11: phase function
    '''
    plt.figure()
    axx=plt.gca();figg=plt.gcf()
    axx.plot(P11)
    axx.set_yscale('log')
    figg.show()
    left=input('Left most index:')
    right=input('Right most index:')
    data=P11[left:right]
    xx=ang[left:right]
    plt.figure()
    ax=plt.gca();fig=plt.gcf()
    ax.plot(xx,data)
#    ax.set_yscale('log')
    ixx=find_peaks_cwt(data,np.arange(1,50))
    ax.plot(xx[ixx],data[ixx],'r+')
    fig.show()
    print(xx[ixx])
    pb=input('Select Primarybow angle:')
    sb=input('Select Supernumerarybow angle:')
    return pb,sb

def scat_ang(SZA,VZA,RAA):
    '''
    Return the scattering angle in degrees
    SZA,VZA : Source(ray)_Zenith_Angle, Viewing_Zenith_Angle (in degrees)
    RAA : Source Azimuth Angle - (Apparent)Scattering Azimuth
    '''
#    RAA=0-180-RAA
# scat_a=cos(VZA)*cos(SZA)+sin(VZA)*sin(SZA)*cos(SAA-VAA)
    SZA=np.deg2rad(SZA);
    VZA=np.deg2rad(VZA)
    RAA=np.deg2rad(RAA)
    
    scat_a=np.arccos(np.cos(SZA)*np.cos(VZA)+np.sin(SZA)*np.sin(VZA)*np.cos(RAA))
    scat_a=np.rad2deg(scat_a)
    return scat_a

if __name__=='__main__':
    f1d=frac_physics('Fractal_1D')
    f1d.generate_fractal(12,90.0,xdist=1.024,xorder=10)
    f2d=frac_physics('Fractal_2D')
    f2d.generate_fractal2D(12,90.0,xydist=1.024,xyorder=10)
    
    fig1,ax1=plt.subplots(subplot_kw={'aspect':'equal'})
    
    ctf=ax1.contourf(f2d.x,f2d.y,f2d.lwp,np.linspace(0,300,50))
#    ctf=ax1.imshow(f2d.tau)
    fig1.colorbar(ctf,ticks=np.arange(0,301,100),extend='both')
    fig1.show()
