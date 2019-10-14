# -*- coding: utf-8 -*-
"""
Created on Mon Oct. 09, 2017

@author: Zhibo Zhang

Almost the same original PDA_Util.py. May be a few additional minor changes. 
Added some comments.   
-cpn-
    
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

def read_rsp_ref(input_rsp_file):
    f=FortranFile(input_rsp_file,'r')
    maxview,maxlayer,maxkern,nview,nlayer,nkern,iint,isrf=f.read_ints(dtype=np.int32)
#    print(maxview,nview)
    phi0,xmu0 = f.read_reals(dtype=np.float)
    tmp = (f.read_reals(dtype=np.float)).reshape(3,maxview)
    thetav = tmp[0,:] 
    tmp = (f.read_reals(dtype=np.float)).reshape(3,maxview)
    rv11,rv21,rv31 = tmp[0,:],tmp[1,:],tmp[2,:] 
    tmp = (f.read_reals(dtype=np.float)).reshape(3,maxview)
    rsrf11,rsrf21,rsrf31 = tmp[0,:],tmp[1,:],tmp[2,:]
    return maxview,phi0,xmu0,nview,thetav,rv11,rv21,rv31,rsrf11,rsrf21,rsrf31

class PDA_Wavelength(object):
	def __init__(self, NLAM, LAM):
		'''
		Input:
		NLAM: number of incident wavelength
		LAM:  array with NLAM elements corresponding to wavelength [um]
		'''
		self.NLAM = NLAM
		self.LAM  = np.array(LAM)

class PDA_Atmos(object):
	def __init__(self, NLAYER, DELP, PDA_Wavelength, NGAS=0, TAUABS=0.0):
		'''
		Input:
		NLAYER: number of layers
		DELP:  array with NLAYER elements corresponding to layer pressure thickness [mb]
		NGAS: If this parameter is set to be greater than 0 then absorption optical depths for gases need to be provided for each wavelength and vertical layer
		'''
		self.NLAYER = NLAYER
		self.DELP   = np.array(DELP)
		self.NGAS   = NGAS 
		self.NLAM   = PDA_Wavelength.NLAM
		self.LAM    = PDA_Wavelength.LAM
		self.TAUABS = np.array(TAUABS) 
		
		

class PDA_Surface(object):
	def __init__(self, PDA_Wavelength, ALBEDO, SURFFILE, OUTFILE):
		'''
		Input:
		NLAM: number of incident wavelength
		LAM:  array with NLAM elements corresponding to wavelength [um]
		ALBEDO: surface albedo at each wavelength
		SURF_FILE: files that contain the surface reflectance when Albedo <0.0
		Albedo>=0.0 - Use a Lambertian surface albedo given by the value of Albedo
		-1.0 < Albedo < 0.0 - Lower boundary is read in from a file but does not have an analytic form and so cannot be removed and then added back in analytically
		Albedo < -1.0 - Lower boundary is read in from a file and does have an analytic form so its direct beam contribution can be removed and then added back in analytically
		'''
		self.NLAM = PDA_Wavelength.NLAM
		self.LAM  = PDA_Wavelength.LAM
		self.ALBEDO = np.array(ALBEDO)
		self.SURFFILE = SURFFILE
		self.OUTFILE  = OUTFILE




class PDA_SS(object):
	def __init__(self, LAM, PIZERC, ANG, P11, F22, F33, F44, F12, F34,\
				 A=1.0, B=0.1, G=np.pi, NR=1.33, NI=0.0,\
				 KEXT = 1.0 ,KSCA =1.0 ,COSBAR=0.8):
		'''
		Inputs:
		KEXT: extinction coefficient
		KSCA: scattering coefficient
		PIZERC: single scattering albedo
		COSBAR: asymetry factor
		Ang: scattering angle [degree]
		P11: phase function
		F22, F33, F44, F12, F34: ratio of scattering phase fuctions to P11

		'''
		self.LAM=LAM
		self.A=A
		self.B=B
		self.G=G
		self.NR=NR
		self.NI=NI
		self.KEXT = KEXT
		self.KSCA = KSCA
		self.PIZERC= PIZERC
		self.COSBAR = COSBAR
		self.ANG = np.array(ANG)
		self.P11 = np.array(P11)
		self.F22 = np.array(F22)
		self.F33 = np.array(F33)
		self.F44 = np.array(F44)
		self.F12 = np.array(F12)
		self.F34 = np.array(F34)
  
	def plot_PM(self,fig1_ttl='PM'):
		fig1,ax1=plt.subplots(3,2,sharex=True)
		ax1[0,0].plot(self.ANG,self.P11);ax1[0,0].set_title('P11')
		ax1[0,0].set_yscale('log')
		ax1[0,1].plot(self.ANG,self.F22);ax1[0,1].set_title('P22/P11')
		ax1[1,0].plot(self.ANG,self.F33);ax1[1,0].set_title('P33/P11')
		ax1[1,1].plot(self.ANG,self.F44);ax1[1,1].set_title('P44/P11')
		ax1[2,0].plot(self.ANG,self.F12);ax1[2,0].set_title('P12/P11')
		ax1[2,1].plot(self.ANG,self.F34);ax1[2,1].set_title('P34/P11')
 
          
		ax1[2,0].set_xlabel('Scattering angle')
		ax1[2,1].set_xlabel('Scattering angle')
		fig1.suptitle(fig1_ttl)
		fig1.show()
 
class PDA_Aerosol(object):
	def __init__(self, NTYPE, NSD, A, B, R1, R2):
		'''
		NTYPE: number of aerosol/cloud type
		NSD:	
		1 - Two parameter Gamma distribution H&T 2.56
		2 - Bimodal gamma distribution H&T 2.59
		3 - Log normal distribution H&T 2.60; A=RG,B=SIGMA
		4 - Power law distribution N(R) = CONST*R**(-A)  FROM RMIN=B TO RMAX=C  AND ZERO OTHERWISE
		5 - Deirmendjian's modified Gamma after de Rooij and van der Stap - Astron and Astrophys. 131 pp237-248 (1984) A=alpha,B=gamma,C=rc
		6 - Read in phase matrices from individual files provided for each aerosol type and wavelength
		7 - Read in scattering database  that also provides the name of the precalculated Fourier decomposition of the same phase matrices
		R1:		Lower limit for size distribution integration
		R2:		Upper limit for size distribution integration
		'''
		self.NTYPE= NTYPE
		self.NSD  = np.array(NSD)
		self.A    = np.array(A)
		self.B    = np.array(B)
		self.R1   = np.array(R1)
		self.R2   = np.array(R2)

	def get_Ref_Index(self,PDA_Wavelength, NR, NI, SS_FILE):
		'''
		NR: real part of refractive index (Ncol: number of wavelength, Nrow:Number of aerosol type)
		NI: imag part of refractive index (Ncol: number of wavelength, Nrow:Number of aerosol type)


		'''
		self.NLAM = PDA_Wavelength.NLAM
		self.LAM  = PDA_Wavelength.LAM
		
		self.NR   = np.array(NR)
		self.NI   = np.array(NI)
		self.SS_FILE = np.array(SS_FILE)


	def get_Loading(self,PDA_Atmos, NDZ1):
		'''
          NDZA: [[species1],[species2],[species3],..]
              [species]=[bottom_layer,bottom_layer+1,bottom_layer+2,...]
        
              -cpn-
		'''
		self.NLAYER = PDA_Atmos.NLAYER
		self.NDZ1 = np.array(NDZ1)


	

class PDA_Parameters(object):
	def __init__(self, DPHI,XMU0,\
		         MCAP=90, NCAP=18, NCAP2=18, NTAU=24, NTAU2=24, MTOT=15, N3BY3=7,\
		         QSTOP= 1.000E-20, QSTOP2= 1.000E-20, NTHETA =  56,\
		         NEXTRA=  1, EXTRA_MU = 1.0, NPHI = 256,\
		         ERRBNDRP=2.000E-4 ,ERRBNDRI=1.000E-4,NPERS =  2,\
		         IPRT=1,IREMV=1):
		 self.DPHI     = DPHI
		 self.XMU0     = XMU0
		 self.MCAP     = MCAP 
		 self.NCAP     = NCAP 
		 self.NCAP2    = NCAP2 
		 self.NTAU     = NTAU 
		 self.NTAU2    = NTAU2 
		 self.MTOT     = MTOT 
		 self.N3BY3    = N3BY3
		 self.QSTOP    = QSTOP 
		 self.QSTOP2   = QSTOP2
		 self.NTHETA   = NTHETA
		 self.NEXTRA   = NEXTRA 
		 self.EXTRA_MU = EXTRA_MU
		 self.NPHI     = NPHI
		 self.ERRBNDRP = ERRBNDRP
		 self.ERRBNDRI = ERRBNDRI
		 self.NPERS    = NPERS
		 self.IPRT     = IPRT
		 self.IREMV    = IREMV
		  
def Write_Drive_File(filename, PDA_Parameters,PDA_Wavelength,PDA_Atmos, PDA_Surface, PDA_Aerosol):
	
	f = open(filename, 'w')
	f.write("Set up parameters for doubling/adding code: quadrature etc.\n") #7(I3,7X)
	fmt = "{:3d}       "*7 + "\n"
	f.write(fmt.format(PDA_Parameters.MCAP,PDA_Parameters.NCAP,PDA_Parameters.NCAP2,\
		               PDA_Parameters.NTAU,PDA_Parameters.NTAU2,PDA_Parameters.MTOT,PDA_Parameters.N3BY3))
	fmt = "        {:10.2e}"*2 + "\n" #2(8X,E10.2E2)
	f.write(fmt.format(PDA_Parameters.QSTOP,PDA_Parameters.QSTOP2))
	fmt="GAUSS DIVISIONS AND WEIGHTS FOR NTHETA ={:4d}\n"
	f.write(fmt.format(PDA_Parameters.NTHETA))
	fmt = "NEXTRA={:3d}\n"
	f.write(fmt.format(PDA_Parameters.NEXTRA))
	fmt = "{:10.5f}\n"
	f.write(fmt.format(PDA_Parameters.EXTRA_MU))
	fmt = "GAUSS DIVISIONS AND WEIGHTS FOR PHI    =  {:4d}\n"
	f.write(fmt.format(PDA_Parameters.NPHI))
	f.write("Relative solar azimuth and cosine_solar_zenith\n")
	fmt = "{:10.5f}" * 2 + "\n"
	f.write(fmt.format(PDA_Parameters.DPHI,PDA_Parameters.XMU0))
	fmt = "ERRBNDR={:10.2e}ERRBNDP={:10.2e}\n"
	f.write(fmt.format(PDA_Parameters.ERRBNDRP,PDA_Parameters.ERRBNDRI))
	fmt = "NPERS ={:3d}\n"
	f.write(fmt.format(PDA_Parameters.NPERS))
	f.write("Number of aerosol types, wavelengths, vertical layers and flag for absorption\n")
	fmt = "       {:4d}"*6 + "\n"
	f.write(fmt.format(PDA_Aerosol.NTYPE, PDA_Wavelength.NLAM, PDA_Atmos.NLAYER,PDA_Atmos.NGAS,PDA_Parameters.IPRT,PDA_Parameters.IREMV ))
	f.write("Definitions of size distributions\n")
	fmt = "       {:8.5f}" * 4 +"      {:4d}\n"
	for it in range(PDA_Aerosol.NTYPE):
		f.write(fmt.format(PDA_Aerosol.A[it],PDA_Aerosol.B[it],\
			    PDA_Aerosol.R1[it],PDA_Aerosol.R2[it],\
			    PDA_Aerosol.NSD[it]))
	f.write("Wavelengths lower bcs and output files\n")
	fmt = "       {:8.5f}" * 2 +"{:40s} " *2 + "\n"
	print(PDA_Surface.LAM.shape)
	for il in range(PDA_Surface.NLAM):
		f.write(fmt.format(PDA_Surface.LAM[il],    PDA_Surface.ALBEDO[il],\
			               PDA_Surface.SURFFILE[il],PDA_Surface.OUTFILE[il]))
	f.write("Refractive indices for each aerosol type and each wavelength, or an input file for each\n")
	fmt = "       {:8.5f}       {:8.3f}" + "  {:80s}\n"
	print(PDA_Aerosol.NR.shape)
	for il in range(PDA_Aerosol.NLAM):
		for it in range(PDA_Aerosol.NTYPE):
			f.write(fmt.format(PDA_Aerosol.NR[il,it],PDA_Aerosol.NI[il,it],PDA_Aerosol.SS_FILE[il,it]))
	f.write("List of pressures and aerosol loads by type for each layer (across the row)\n")
	
	f.write("       ")
	for j in range(PDA_Atmos.NLAYER): f.write(" {:6.1f}".format(PDA_Atmos.DELP[j]))
	f.write("\n")
	for it in range(PDA_Aerosol.NTYPE):
		f.write("       ")
		for j in range(PDA_Atmos.NLAYER):f.write(" {:10.2e}".format(PDA_Aerosol.NDZ1[it,j]))
		f.write("\n")
	if PDA_Atmos.NGAS > 0:
		f.write("Absorption optical depths for each wavelength and layer if NGAS>0\n")
		for il in range(PDA_Atmos.NLAM):
			f.write("       ")
			for j in range(PDA_Atmos.NLAYER): f.write(" {:10.2e}".format(PDA_Aerosol.TAUABS[il,j]))
		f.write("\n")
	 
	f.close()
	print(filename+' SAVED.')
		
def Write_SS_File(filename, PDA_SS, Info1 = ' ', Info2 =' ',Info3 = ' ', Info4=' '):
	
	f = open(filename, 'w')
	f.write(Info1+'\n')
	f.write(Info2+'\n')
	f.write(Info3+'\n')
	print(PDA_SS.A, PDA_SS.B, PDA_SS.G, PDA_SS.LAM, PDA_SS.NR, PDA_SS.NI)
	fmt="Reff={:12.5e}    Veff={:6.4f}    Area={:12.5e}    Lambda={:12.5e}    NR={:12.5e}    NI={:12.5e}\n"
	f.write(fmt.format(PDA_SS.A, PDA_SS.B, PDA_SS.G, PDA_SS.LAM, PDA_SS.NR, PDA_SS.NI))
	f.write('KEXT={:12.5e}    KSCA={:12.5e}    PIZERC={:8.6f}    COSBAR={:8.6f}\n'.\
		    format(PDA_SS.KEXT,PDA_SS.KSCA,PDA_SS.PIZERC,PDA_SS.COSBAR))
	f.write('Angle (degs)          P11          P22/P11         P33/P11         P44/P11         P12/P11         P34/P11\n')

	for i in range(PDA_SS.ANG.size):
		f.write('{:12.5e}    {:12.5e}    {:12.5e}    {:12.5e}    {:12.5e}    {:12.5e}    {:12.5e}    \n'.\
		         format(PDA_SS.ANG[i],PDA_SS.P11[i],PDA_SS.F22[i],PDA_SS.F33[i],PDA_SS.F44[i],PDA_SS.F12[i],PDA_SS.F34[i]))

	f.close()
	print(filename+' SAVED!')

def Write_Surf_File(filename, PDA_Parameters, LAM, ISURF, Fpar, Surf_filename, IGAUSS=2):
	f = open(filename, 'w')
	f.write(Surf_filename+'\n')
	fmt="{:3d}       " * 7+"\n"
	f.write(fmt.format(PDA_Parameters.MCAP, PDA_Parameters.N3BY3,PDA_Parameters.NTHETA,\
					   PDA_Parameters.NEXTRA,\
	                   PDA_Parameters.NPHI, ISURF, IGAUSS))
	fmt = "{:10.5f}\n"
	f.write(fmt.format(PDA_Parameters.EXTRA_MU))
	fmt="       {:9.6f}" *4 + "\n"
	f.write(fmt.format(LAM,Fpar[0],Fpar[1],Fpar[2]))
	f.write(fmt.format(Fpar[3],Fpar[4],Fpar[5],Fpar[6]))
	fmt="       {:9.6f}" *2 + "\n"
	f.write(fmt.format(Fpar[7],Fpar[8]))
	f.close()

def test_pda_driver():
	

	par   = PDA_Parameters(90.0, 0.5) 
	wvl   = PDA_Wavelength(1,[0.86])
	atmos = PDA_Atmos(1, [1.0], wvl)
	surf  = PDA_Surface(wvl, [0.03], ['surf_wvl1'], ['out_wvl1'])
	aer   = PDA_Aerosol(1, [6], [10.0], [0.03], [1], [60.0])
	aer.get_Ref_Index(wvl,[[1.33]],[[0.0]],[['benchmark_aerosol.PDA']])
	aer.get_Loading(atmos,[[0.3262]])
	d=np.loadtxt('Kokhanovsky_benchmark_aerosol.PDA',skiprows=6)
	benchmark_ss = PDA_SS(wvl.LAM[0], 1.0, d[:,0], d[:,1], d[:,2], d[:,3], d[:,4], d[:,5], d[:,6])
	Write_SS_File('benchmark_aerosol.PDA',benchmark_ss)

	Write_Drive_File('test_driver.info', par, wvl,atmos, surf, aer)

	import os
	os.system("PDA_new ./ test_driver.info ./")
	maxview,phi0,xmu0,nview,thetav,\
	rv111,rv211,rv311, rsrf111,rsrf211,rsrf311=read_rsp_ref('out_wvl1.rsp')
	plt.figure()
	plt.plot(thetav[0:nview/2],rv111[0:nview/2],c='r')
	plt.plot(-thetav[nview/2:],rv111[nview/2:],c='b')

	plt.figure()
	plt.plot(thetav[0:nview/2],-rv211[0:nview/2],c='r')
	plt.plot(-thetav[nview/2:],-rv211[nview/2:],c='b')
	#plt.plot(thetav,-rv212,c='b')

	plt.show()

def test_surf_driver():
	par   = PDA_Parameters(90.0, 0.5) 
	wvl   = 0.67
	ISURF = 1
	Fpar  = np.array([0.0, -0.1, 0.0, 0.0, 0.0,1.33, 0.0, 0.0207, 1.0])

	Write_Surf_File('test.info', par, wvl, ISURF, Fpar, 'CoxMunk_wind7.5m_wl0.67_nr1.33_ni0')


 

if __name__ == '__main__':
	test_surf_driver()
