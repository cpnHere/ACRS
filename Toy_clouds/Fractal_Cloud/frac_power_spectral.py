#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Oct  2 15:30:51 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Power spectral analysis for the spectral cloud.
"""

import matplotlib.pyplot as plt
import numpy as np
import cpnMicrophylib
import cpnLES_MSCARTlib   
from textwrap import wrap
import pickle
from scipy.fftpack import fft

class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]

with open('fractal_cld_b865re12ve05_x40km.pkl', "rb") as f:
    pkl=pickle.load(f)
    MSCART_field=pickle.load(f)
    fcld=pickle.load(f)

fcld.plot_cld()

# Reading Fractal cloud case for Dan's thesis ################################
fcld_djm=cpnMicrophylib.frac_physics('physical_scene_case_new.h5')   
fcld_djm.read() 
fcld_djm.plot_cld()

#Reflectances from MSCART
fracSZAx=cpnLES_MSCARTlib.POLCARTdset('fracSZAx','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/')
fracSZAx.readPOLCARThdf5('fractal_cld_b865re12ve05_x40km_MSCART_SZA120_SAA000_VAA000plus_NPH1e5.hdf5',dpath='results/fracb865/')
#frac=cpnLES_MSCARTlib.LES_field(fracSZAx.fname.split('_MSCART',2)[0]+'.nc')#DYCOM field file for MSCART
#frac.readLES_field()

frac1D=cpnLES_MSCARTlib.POLCARTdset('frac1D','/umbc/xfs1/zzbatmos/users/charaj1/Toy_clouds_MSCART/Fractal_Cloud/1Druns/fracb865_bins/')
frac1D.readPOLCARThdf5('fractal_cld_b865re12ve05_x40km_1D_bins_MSCART_SZA120_SAA000_VAA000plus_NPH1e5.hdf5',\
                              dpath='1Druns/results/fracb865_bins/')

#Fourier Tranform and power spectral analysis
def linear_fit(ax,x_fit,y_fit,lst='k-',linewidth=2.0,points=True):
    pz=np.polyfit(np.log10(x_fit[1:]),np.log10(y_fit[1:]),1)
    p=np.poly1d(pz)
    ax.plot(10**np.log10(x_fit),10**p(np.log10(x_fit)),lst,label=r'$\beta=$%0.2f'%pz[0],linewidth=linewidth)
    print('Beta: %0.2f'%pz[0])
    if points:
        ax.plot(x_fit,y_fit,'k*',markersize=10,label='Ocatave binned average')
def find_octave_binned(xf,YF,ax=None):
    #binning points into octave bins
    Xorder=int(np.log2(xf.size))
    i_lef=0;x_fit=[];y_fit=[]
    print('Octave binning.....')
    for i in np.arange(0,Xorder,1):
        bin_width=2**i
#        print('%d to %d'%(i_lef,i_lef+bin_width))
        x_fit=np.append(x_fit,np.mean(xf[i_lef:i_lef+bin_width]))
        y_fit=np.append(y_fit,np.mean(YF[i_lef:i_lef+bin_width]))
        if not(ax==None):
            ax.plot(xf[i_lef:i_lef+bin_width], YF[i_lef:i_lef+bin_width],'.')
        i_lef=i_lef+bin_width
    return x_fit,y_fit
def plot_psd(ax,fig_ttl,xf,YF):
    x_fit,y_fit=find_octave_binned(xf,YF,ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'1/x $(km^{-1})$')
    ax.set_ylabel('E[1/x]')
    ax.set_title("\n".join(wrap(fig_ttl,60)))

    return x_fit,y_fit

#Optical depth field beta
fig1_ttl=MSCART_field.fname.split('.',1)[0]+'_tau'
N=fcld.x.size
T=40.0/(2**12)#40km/16384
xf=np.linspace(0.0,1.0/(2.0*T),N/2)
yf=fft(np.squeeze(fcld.tau))
fig1,ax1=plt.subplots()
YF=(np.abs(yf[0:N//2]))**2
x_fit,y_fit=plot_psd(ax1,fig1_ttl,xf,YF)
linear_fit(ax1,x_fit,y_fit)
ax1.legend(prop={'size':10})
ax1.grid()

#Getting 1.66
tempfcld=cpnMicrophylib.frac_physics('higer_order_temp')

lwp=90#gm^{-2}
re=12.0#microns

tempfcld.f=0.51#fractal parameter fn=fc^n
tempfcld.generate_fractal(re,lwp,xorder=15,xdist=40.0)
N=tempfcld.x.size
T=40.0/(2**15)#40km/16384
xf=np.linspace(0.0,1.0/(2.0*T),N/2)
yf=fft(np.squeeze(tempfcld.tau))
YF=(np.abs(yf[0:N//2]))**2


fig1_ttl='x24'



fig1,ax1=plt.subplots()

x_fit,y_fit=plot_psd(ax1,'title',xf,YF)
linear_fit(ax1,x_fit,y_fit,linewidth=1.0)
ax1.legend(prop={'size':10})
ax1.grid()
ax1.set_xscale('log')
ax1.set_yscale('log')



x_fit_set=[];y_fit_set=[]
for i in np.arange(0,50):
    tempfcld.f=0.51#fractal parameter fn=fc^n
    tempfcld.generate_fractal(re,lwp,xorder=24,xdist=40.0)
    N=tempfcld.x.size
    T=40.0/(2**24)#40km/16384
    xf=np.linspace(0.0,1.0/(2.0*T),N/2)
    yf=fft(np.squeeze(tempfcld.tau))
    YF=(np.abs(yf[0:N//2]))**2
    x_fit,y_fit=find_octave_binned(xf,YF)
    x_fit_set=np.append(x_fit_set,x_fit[1:])
    y_fit_set=np.append(y_fit_set,y_fit[1:])

fig2,ax2=plt.subplots()

linear_fit(ax2,x_fit_set,y_fit_set,linewidth=0.5)
ax2.legend(prop={'size':10})
ax2.grid()
ax2.set_xscale('log')
ax2.set_yscale('log')



#====================== combined 4 figures ====================================
fig3,ax3=plt.subplots(2,2,figsize=(12,12))
fig3_ttl=MSCART_field.fname.split('.',1)[0]+'_'+fracSZAx.fname.split('.',1)[0]+'_I3D_and_I1D'

fig3.suptitle("\n".join(wrap(fig3_ttl,60)))
ax3[0,0].plot(fcld.x,fcld.tau,'g',label=r'$\tau$')
ax3[0,0].set_ylabel(r'Optical depth')
ax3[0,0].set_xlabel('x (km)')
ax3[0,0].set_title(r'$\tau$')

N=fcld.x.size
T=40.0/(2**15)#40km/16384
xf=np.linspace(0.0,1.0/(2.0*T),N/2)
yf=fft(np.squeeze(fcld.tau))
fig1,ax1=plt.subplots()
YF=(np.abs(yf[0:N//2]))**2
x_fit,y_fit=plot_psd(ax3[0,1],r'$\tau$',xf,YF)
linear_fit(ax3[0,1],x_fit,y_fit)
ax3[0,1].legend(prop={'size':10})
ax3[0,1].grid()

#3D reflectance
I3D=np.mean(fracSZAx.MeanPRad[fracSZAx.VZA==0,:,0],axis=0)
yf=fft(np.squeeze(I3D))
YF=(np.abs(yf[0:N//2]))**2
x_fit,y_fit=plot_psd(ax3[1,0],r'$I_{3D}$',xf,YF)
linear_fit(ax3[1,0],x_fit,y_fit)
#linear_fit(ax3[1,0],x_fit[0:7],y_fit[0:7],'g-')
#linear_fit(ax3[1,0],x_fit[6:],y_fit[6:],'b-')
ax3[1,0].grid()
ax3[1,0].legend(prop={'size':10})

#1D reflectances
I1D=np.mean(frac1D.MeanPRad[frac1D.VZA==0,:,0],axis=0)
yf=fft(np.squeeze(I1D))
YF=(np.abs(yf[0:N//2]))**2
x_fit,y_fit=plot_psd(ax3[1,1],r'$I_{1D}$',xf,YF)
linear_fit(ax3[1,1],x_fit,y_fit)
#linear_fit(ax3,x_fit[0:7],y_fit[0:7],'g-')
#linear_fit(ax3,x_fit[6:],y_fit[6:],'b-')
ax3[1,1].grid()
ax3[1,1].legend(prop={'size':10})


#fcld_djm PSD
#fig2_ttl=fcld_djm.fname.split('.',1)[0]+'_tau'
#N=fcld_djm.x.size
#T=0.01
#xf=np.linspace(0.0,1.0/(2.0*T),N/2)
#yf=fft(np.squeeze(fcld_djm.tau))
#fig2,ax2=plt.subplots()
#YF=(np.abs(yf[0:N//2]))**2
#plot_psd(ax2,fig2_ttl,xf,YF)
#ax2.grid()

def savefig(fig,fig_ttl):
    fig.savefig('figures/'+fig_ttl+'.png',format='png',dpi=200)
    print('figures/'+fig_ttl+'.png SAVED!')
#plt.show()
