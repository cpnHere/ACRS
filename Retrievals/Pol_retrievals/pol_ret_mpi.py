#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Thu Apr 11 15:39:52 2019
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
To do polarimetric retrievals with Python MPI
Assign appropriate nmpiJ depending on the available resources
"""
from mpi4py import MPI
import numpy as np
import time,os,sys
from cpnCommonlib import load_obj, save_obj


start=time.time()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nmpiJ=16 # Number of ranks
print('My rank is ',rank)
if rank==0:
    df_name=str(sys.argv[1])   
    data=load_obj(df_name)
else:
    data=None
data = comm.bcast(data, root=0)

x=data['x']
y=data['y']
P=data['P']
ygabc=data['ygabc']
nxy=y.shape[2]#Number of pixels along one side
Nperjob=int(nxy**2/nmpiJ)#Number of jobs per rank
ret_Re=np.zeros((nxy,nxy),dtype=float)
ret_Ve=np.zeros((nxy,nxy),dtype=float)
Qls   =np.zeros((nxy,nxy),dtype=float)
Rsq   =np.zeros((nxy,nxy),dtype=float)
flags =np.zeros((nxy,nxy),dtype=int)
abc   =np.zeros((nxy,nxy,3),dtype=float)
yAll  =np.zeros((nxy,nxy,x.size),dtype=float)
jobs=int(rank)*Nperjob+np.arange(0,Nperjob,1)
if P.method=='Breon':
    from cpnRetrievalslib import fitBreon_noRsq as do_fitting
#    print("\t First job:%d"%(jobs[0]))
#    print("\t Last job:%d"%(jobs[-1]))
for n in jobs:
    i=int(n/nxy)
    j=n%nxy
    ret_Re[i,j],ret_Ve[i,j],abc[i,j,:],Qls[i,j],Rsq[i,j],flags[i,j]=do_fitting(x,np.squeeze(y[P.Q_a1:P.Q_a2,i,j]),P,ygabc[i,j,:])
    yAll[i,j,:]=y[P.Q_a1:P.Q_a2,i,j]
if rank>0:
    out={'ret_Re':ret_Re,'ret_Ve':ret_Ve,'abc':abc,'Qls':Qls,'Rsq':Rsq,'yAll':yAll,'flags':flags}
    comm.send(out, dest=0)
    
if rank==0:
    for id in np.arange(nmpiJ-1):
        out=comm.recv(source=id+1)
        abc=abc+out['abc']
        ret_Re=ret_Re+out['ret_Re']
        ret_Ve=ret_Ve+out['ret_Ve']
        Qls   =Qls+out['Qls']
        Rsq   =Rsq+out['Rsq']
        yAll  =yAll+out['yAll']
        flags =flags+out['flags'] 
    saveO={'ret_Re':ret_Re,'ret_Ve':ret_Ve,'abc':abc,'Qls':Qls,'Rsq':Rsq,'yAll':yAll,'x':x,'flags':flags}
    save_obj(saveO,data['savename']+'_MPI',rp=True)
    os.system('rm '+'mpi_data_'+data['savename']+'.pkl')
    end=time.time()
    print('%0.2f mins elapsed'%(end/60-start/60))
    