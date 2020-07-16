#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Apr 13 08:39:22 2017
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To run Mie code and do the averaging. (Intermediate version)
PSDs object:
    To handle size distributions (conversions and stuff...)
11/06/2017:
    n_lnV2n_lnN added.
04/03/2018:
    save_bulk_Mie() and read_bulk_Mie() added to bulk_Mie object 
"""
import matplotlib.pyplot as plt
import netCDF4, h5py
import numpy as np
from subprocess import call

class MieSet(object):
    '''
    For Mie runs
    out_name: Mie ouput name without extension
    path(None): file path
    '''
    def __init__(self, out_name,path=None):
        self.d=[]; self.ang=[]; self.wvl=[]; self.alb=[]; self.qe=[]
        self.P11=[]; self.P12=[]; self.P33=[]; self.P34=[]; self.nr=[];
        self.ni=[]; self.g=[]; 
        self.out_name=out_name
        self.path=path
        self.avqe=[]
        self.avalb=[]
        self.avCex=[]
        self.CDFavCex=[]
        
    def runMie(self, lam, nr, nc, D ):
        print('Mie code .....')
        if np.size(lam)==1:
            self.out_name=sgl_refidx(self.out_name,lam,nr,nc,D)
        else:
            mul_refidx(self.out_name,lam, nr, nc, D)
    
    def readMie(self,):
        if self.path==None:
            file=self.out_name
        else:
            file=self.path+self.out_name
        self.d, self.ang, self.wvl, self.alb, self.qe, self.P11, self.P12, self.P33, self.P34, \
        self.nr, self.ni, self.g= read_mie_cmplt(file+'.nc')

def read_mie_ncdf(fname):
    # read output from the Mie code
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
    
def read_mie_cmplt(fname):
    '''
    read all parameters from the Mie code
    uses read_mie_ncdf()
    '''
    d, ang, wvl, alb, qe, P11, P12, P33, P34=read_mie_ncdf(fname)
    f = netCDF4.Dataset(fname,'r')
    nr=f.variables['Refr_real'][:]
    ni=f.variables['Refr_img'][:]
    g=f.variables['AsymmetryFactor'][:]
    f.close()
    return d, ang, wvl, alb, qe, P11, P12, P33, P34, nr, ni, g

def sgl_refidx(out_name,lam,nr,nc,D):
    #Single ref index version. 
    #This version save the out.nc file based on ref index.  
    #size parameter x=piD/lam
    #----------------------------------
#    lam=str(0.532)#wavelength
#    nr=str(1.337115)#real n
#    nc=str(1.818175e-9)#complex n
#    D=np.linspace(0,80,10000)#diameter
    D=np.hstack((np.array([10000]),D))
    #--------------------------------------
    #D=np.hstack((np.array([1]),D))

    out_name=out_name+'nr'+nr+'_nc'+nc+'.nc'
    
    f1=open("filename.dat","w")
    f1.write(out_name+'\n')
    f1.close()
    
    f2=open("wl_ref.dat","w")
    f2.write('           1\n  '+lam+'       '+nr+'     '+nc+'\n')
    f2.close()
    
    np.savetxt('size.dat',D,delimiter='\n',fmt='%.7f')
    
    call(["./mie_single_size.exe"])
    
    print('File name: '+out_name)
    print('nr: '+str(nr))
    print('ni: '+str(nc))
    print('D : from '+str(D[1])+' to '+str(D[-1]))
    f1.close()
    f2.close()
    return out_name.split('.n',1)[0]
    
def mul_refidx(out_name,lam, nr, nc, D):
    #Multiple reffractive indices
    #out_name is the output file name

#    lam=[3.7,11]
#    nr=[1.356769,1.127233 ]#Real part of n
#    nc=[0.3589187E-02,0.9690998E-01]#imag. part of n
#    D=np.linspace(0,60,10000)#diameter
    #D=[3.1415]
    D=np.hstack((np.array([np.size(D)]),D))
    #D=np.hstack((np.array([1]),D))
    out_name=out_name+'.nc'
    
    f1=open("filename.dat","w")
    f1.write(out_name+'\n')
    f1.close()
    
    f2=open("wl_ref.dat","w")
    f2.write('           '+str(np.size(nr))+'\n')
    for i in range(0,np.size(nr)):
        f2.write('  '+str(lam[i])+'       '+str(nr[i])+'     '+str(nc[i])+'\n')
        print(str(i))
    f2.close()

    np.savetxt('size.dat',D,delimiter='\n',fmt='%.7f')
    call(["mie_single_size.exe"])
     
    print('File name: '+out_name)
    print('nr: '+str(nr))
    print('ni: '+str(nc))    
    print('D : from '+str(D[1])+' to '+str(D[-1]))
    f1.close()
    f2.close()
    
class PSDs(object):
    '''
    Make sure D is a numpy float when you directly input a d array that read from 
    Mie output *.nc file.
    Convert input D by using np.asarray(arr,dtype=float)
    D: diameter in microns
    '''
    def __init__(self,tag,D=None):
        self.tag=tag
        if D is None:
            self.D=[]
            self.r=[]
            self.n_lnV=[]
        else:
            self.D=D
            self.r=D/2
    def mod_gamma_norm(self,re,ve):
        '''
        re in microns
        Out put normalized modified gamma distributon for re,ve (Hansen and Travis, 1974)
        '''
        f=self.r**((1.0-3.0*ve)/ve)*np.exp(-self.r/(re*ve))
        C=np.trapz(f,self.r)
        self.re=re;self.ve=ve
        return f/C
    def set_mod_gamma_norm(self,re,ve):
        '''
        මෙහෙමයි....
        re in microns
        Setup modified gamma distribution
        dN/dr
        '''
        self.n_N=self.mod_gamma_norm(re,ve)
        
    def freshagedsal(self,):
        #A data set of PSDs2 
        print('Reading Probability Size Distribution Functions......')
        self.fpsd='ryder_output.txt'
        self.fl=open(self.fpsd,'r')
        self.str1=self.fl.readline();
        self.data=self.fl.readlines()[:]
        self.psd=np.loadtxt(self.data)
        self.fresh_n_lnV=self.psd[:,1]
        self.aged_n_lnV=self.psd[:,2]
        self.sal_n_lnV=self.psd[:,3]
        self.r=self.psd[:,0]/2
        self.D=self.psd[:,0]  
        
    def convertPSD(self):
        #converting factors from n_lnV
        self.fn_lnD=6/np.pi/self.D**3#Particle number
        self.fn_SlnD=1.5/self.D#Crossection
        self.fn_D=6/np.pi/self.D**4
        self.fn_r=self.fn_D*2
        
    def setlogNorm(self, Cf,Cc,sigf,sigc,rf,rc):
        self.n_lnV=Cf/np.sqrt(2*np.pi)/np.log(sigf)*\
            np.exp(-(np.log(self.r)-np.log(rf))**2/2/np.log(sigf)**2)+\
            Cc/np.sqrt(2*np.pi)/np.log(sigc)*\
            np.exp(-(np.log(self.r)-np.log(rc))**2/2/np.log(sigc)**2)
        
    def normalizefreshagedsal(self,):
        print('Normalizing to unit volume.')
        self.fresh_n_lnV=self.fresh_n_lnV/np.trapz(self.fresh_n_lnV,np.log(self.D))
        self.aged_n_lnV=self.aged_n_lnV/np.trapz(self.aged_n_lnV,np.log(self.D))
        self.sal_n_lnV=self.sal_n_lnV/np.trapz(self.sal_n_lnV,np.log(self.D))

    def normalize(self,):
        self.n_lnV=self.n_lnV/np.trapz(self.n_lnV,np.log(self.r))
    
    def n_lnV2n_lnN(self,):
        #dN/dlnr
        #V_l(r)dlnr=4/3*pi*r^3*n_l(r)dlnr
        if hasattr(self,'n_lnV'):
            self.n_lnN=3.0/(4.0*np.pi*self.r**3)*self.n_lnV
        else:
            print('Palamuwa n_lnV artha dakwanna :p')
    def findn_N(self,):
        if hasattr(self,'n_lnN'):
            self.n_N=self.n_lnN/self.r
        elif hasattr(self,'n_lnV'):
            self.n_lnV2n_lnN()
            self.n_N=self.n_lnN/self.r
        else:
            print('No n_lnN (or n_lnV) found!')
            
    def findln_lnS(self,):
        #dS/dlnr S is cross section area
        if hasattr(self,'n_lnN'):
            self.n_lnS=np.pi*self.r**2*self.n_lnN
        elif hasattr(self,'n_lnV'):
            self.n_lnV2n_lnN()
            self.n_lnS=np.pi*self.r**2*self.n_lnN
        else:
            print('First define n_lnN or n_lnV!')
    def plot_nN(self,fig=None,ax=None,ls='k-'):
        if fig is None:
            fig,ax=plt.subplots()
            ax.set_xlabel(r'Radius ($\mu m$)')
            ax.set_ylabel(r'Normalized PDF [$cm^{-3}\mu m^{-1}$]')
        ax.plot(self.r,self.n_N,ls)
        fig.show()
        return fig,ax
        
class bulk_Mie(object):
    def __init__(self, fname,fdpath=None,psd=None,Mie=None):
        '''
        To read or creat bulk_Mie object
        fname: filename
        fdpath: data path
        psd: cpnMielib.PSDs object
        Mie: cpnMielig.MieSet object
        '''
        if psd==None and Mie==None:
            self.read_bulk_Mie(fname,fpath=fdpath)
        else:
            self.psd=psd; self.Mie=Mie;
    def cal_bulk_Pmat(self,):
        if not(hasattr(self.psd,'n_N')):
            if hasattr(self.psd,'n_lnV'):
                self.psd.n_lnV2n_lnN()
                self.psd.n_N=self.psd.n_lnN/self.psd.r
                print('n_N not found. Calculated from n_lnV.')
            else:
                print('No n_N (or n_lnV) found!')
        avP11=np.zeros((self.Mie.ang.size,self.Mie.wvl.size))
        avP12=np.zeros((self.Mie.ang.size,self.Mie.wvl.size))
        avP33=np.zeros((self.Mie.ang.size,self.Mie.wvl.size))
        avP34=np.zeros((self.Mie.ang.size,self.Mie.wvl.size))
        print('Computing bulk Pmat........')
        for k in range(0,self.Mie.wvl.size):
        #fb=r**2*qe.T*alb.T*n_r
            fb=self.psd.r**2*self.Mie.qe[:,k].T*self.Mie.alb[:,k].T*self.psd.n_N
            ft=np.einsum('i,ij->ij',fb,self.Mie.P11[:,k,:])
            avP11[:,k]=np.trapz(ft,self.psd.r,axis=0)/np.trapz(fb,self.psd.r)
            ft=np.einsum('i,ij->ij',fb,self.Mie.P12[:,k,:])
            avP12[:,k]=np.trapz(ft,self.psd.r,axis=0)/np.trapz(fb,self.psd.r)
            
            ft=np.einsum('i,ij->ij',fb,self.Mie.P33[:,k,:])
            avP33[:,k]=np.trapz(ft,self.psd.r,axis=0)/np.trapz(fb,self.psd.r)
            
            ft=np.einsum('i,ij->ij',fb,self.Mie.P34[:,k,:])
            avP34[:,k]=np.trapz(ft,self.psd.r,axis=0)/np.trapz(fb,self.psd.r)
        self.bulk_P11=avP11
        self.bulk_P12=avP12
        self.bulk_P33=avP33
        self.bulk_P34=avP34
    def plot_bulk_Pmat(self,):
        '''
        Plot bulk_averaged Phase matrix
        return fig,ax
        '''
        fig,ax=plt.subplots(2,2)
        ax[0,0].set_title('P11')
        ax[0,0].plot(self.Mie.ang,self.bulk_P11)
        ax[0,0].set_yscale('log')
        ax[0,0].legend(self.Mie.wvl.astype('|S5'),loc='best')
        ax[0,1].set_title('P33/P11')
        ax[0,1].plot(self.Mie.ang,self.bulk_P33/self.bulk_P11)
        ax[1,0].set_title('-P12/P11')
        ax[1,0].plot(self.Mie.ang,-self.bulk_P12/self.bulk_P11)
        ax[1,1].set_title('P34/P11')
        ax[1,1].plot(self.Mie.ang,self.bulk_P34/self.bulk_P11)
        fig.show()
        return fig,ax

    def cal_bulk_albnQe(self,):
        print('Averaging for bulk SSA and Qe........')
        avalb=np.zeros(self.Mie.wvl.size)
        avQe=np.zeros(self.Mie.wvl.size)
        for i in range(0,self.Mie.wvl.size):
            #fb=qe*n_V/r   
            fb=self.Mie.qe[:,i]*4.0/3.0*np.pi*self.psd.r**2*self.psd.n_N
            #ft=fb*alb
            ft=fb*self.Mie.alb[:,i]
            avalb[i]=np.trapz(ft,self.psd.r)/np.trapz(fb,self.psd.r)
            fb=np.pi*self.psd.r**2*self.psd.n_N
            ft=self.Mie.qe[:,i]*fb
            avQe[i]=np.trapz(ft,self.psd.r)/np.trapz(fb,self.psd.r)
        self.bulk_alb=avalb
        self.bulk_Qe=avQe
        
    def save_bulk_Mie(self,out_name,obj=None):
        '''
        save bulk_Mie object (including bulk avg SS and psd) as *.hdf5
        out_name: Name without extension
        ***Note that this will not save all the attributes of the Mie object. Add it***
        '''
        if obj==None:
            obj=self
        f=h5py.File(out_name+'.hdf5','w')
        name={1:'bulk_P11',2:'bulk_P12',3:'bulk_P33',4:'bulk_P34'}
        data={1:obj.bulk_P11,2:obj.bulk_P12,3:obj.bulk_P33,4:obj.bulk_P34}
        for i in np.arange(1,5):
            PCentry=f.create_dataset(name[i],data=data[i])
            PCentry.dims[0].label='ang'
            PCentry.dims[1].label='wvl'
            PCentry.attrs['long_name']='Phase_matrix_element'
        
        PCentry=f.create_dataset('bulk_Qe',data=obj.bulk_Qe)
        PCentry.dims[0].label='wvl'
        PCentry.attrs['long_name']='bulk_average_extinction_efficiency'
        
        PCentry=f.create_dataset('bulk_alb',data=obj.bulk_alb)
        PCentry.dims[0].label='wvl'
        PCentry.attrs['long_name']='bulk_average_extinction_albedo'
        
        PCentry=f.create_dataset('ang',data=obj.Mie.ang)
        PCentry.attrs['long_name']='Scattering_angle'
        
        PCentry=f.create_dataset('wvl',data=obj.Mie.wvl)
        PCentry.attrs['long_name']='wavelength'
        
        PCentry=f.create_dataset('n_N',data=obj.psd.n_N)
        PCentry.dims[0].label='r'
        PCentry.attrs['long_name']='dN/dr'
        
        if hasattr(obj.psd,'n_lnN'):
            PCentry=f.create_dataset('n_lnN',data=obj.psd.n_lnN)
            PCentry.dims[0].label='r'
            PCentry.attrs['long_name']='dN/dlnr'
        
        PCentry=f.create_dataset('n_lnV',data=obj.psd.n_lnV)
        PCentry.dims[0].label='r'
        PCentry.attrs['long_name']='dV/dlnV'
        
        PCentry=f.create_dataset('r',data=obj.psd.r)
        PCentry.attrs['long_name']='radii'
        
        print(out_name+'.hdf5 SAVED !!')
        f.close()
    
    def read_bulk_Mie(self,fname,fpath=None):
        '''
        Reading bulk_Mie object *.hdf5 file
        '''
        if fpath==None:
            filename=fname
        else:
            filename=fpath+fname
        

        f=h5py.File(filename,'r')
        self.Mie=MieSet('from_'+fname)
        self.psd=PSDs('from_'+fname)
        self.bulk_P11=f['bulk_P11'][:]
        self.bulk_P12=f['bulk_P12'][:]
        self.bulk_P33=f['bulk_P33'][:]
        self.bulk_P34=f['bulk_P34'][:]
        self.bulk_Qe=f['bulk_Qe'][:]
        self.bulk_alb=f['bulk_alb'][:]
        self.Mie.ang=f['ang'][:]
        self.Mie.wvl=f['wvl'][:]
        self.psd.n_lnV=f['n_lnV'][:]
        self.psd.r=f['r'][:]
        self.psd.n_N=f['n_N'][:]    
        for i in np.arange(1,13):
            try:
                self.psd.n_lnN=f['n_lnN'][:]
            except KeyError:
                print('n_lnN key does not exisit in '+filename)
