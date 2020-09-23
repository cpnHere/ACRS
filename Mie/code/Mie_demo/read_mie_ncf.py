from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import netCDF4

def read_mie_ncdf(fname):
    f = netCDF4.Dataset(fname,'r')
    d = f.variables['Diameter'][:]
    ang = f.variables['PhaseFunctionAngle'][:]
    wvl = f.variables['Wavelength'][:]
    alb = f.variables['SingleScatteringAlbedo'][:,:]
    qe  = f.variables['ExtinctionEfficiency'][:,:]
    P11 = f.variables['P11'][:,:,:]
    P12 = f.variables['P12'][:,:,:]
    P33 = f.variables['P33'][:,:,:]
    P34 = f.variables['P34'][:,:,:]
    f.close()
    return d, ang, wvl, alb, qe, P11, P12, P33, P34

if __name__ =='__main__':

    nf = 'Single_Scattering_Mie.nc'
    d, ang, wvl, alb, qe, P11, P12, P33, P34 = read_mie_ncdf(nf)

    print('number of size',d.size)
    print('number of angle',ang.size)
    print('number of wavelength', wvl.size)
    print('Qe & Alb shape',qe.shape)
    print('P11 shape',P11.shape)

    fig,ax = plt.subplots()
    ax.plot(d,qe[:,0])
    ax.set_xscale('log')
    ax.set_xlim([1,100])
    ax.set_xlabel('size (um)')
    ax.set_ylabel('Qe')

    fig,ax = plt.subplots(2,1)
    ax[0].plot(ang,P11[0,0,:])
    ax[1].plot(ang,-P12[0,0,:]/P11[0,0,:])
    ax[1].set_xlabel('Scattering Angle [degree]')
    ax[0].set_ylabel('normalized P11')
    ax[1].set_ylabel('DOLP -P12/P11')

    fig,ax = plt.subplots(2,1)
    ax[0].plot(ang,P11[-1,0,:]); ax[0].set_yscale('log')
    ax[1].plot(ang,-P12[-1,0,:]/P11[-1,0,:])
    ax[1].set_xlabel('Scattering Angle [degree]')
    ax[0].set_ylabel('normalized P11')
    ax[1].set_ylabel('DOLP -P12/P11')

    plt.show()
