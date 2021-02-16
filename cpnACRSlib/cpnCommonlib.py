#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
********************************************
Created on Mon Apr  2 10:14:27 2018
by
Chamara Rajapakshe
(cpn.here@umbc.edu)
********************************************
Frequently used common python library
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import os,string
import scipy.signal as signal
import scipy.optimize as opt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import pickle, sys
import mpl_toolkits.basemap as bm
import matplotlib as mp

class pkl_classes(object):
    def __init__(self,):
        self.class_names=[]
# Disable
def blockPrint(_sys):
    '''
    To block printing text 
    '''
    _sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint(_sys):
    '''
    To enable printing text
    '''
    _sys.stdout = _sys.__stdout__
def fit_a_curve(xdata,ydata,label1='Data',label2='Fit',fig=None,ax=None,pltData=True):
    '''
    pltData: False to only plot the fit
    Do non-linear fit for x and y data.
    '''
    if fig is None:
        fig,ax = plt.subplots()
    def func(x, a, b, c):
        '''
        This is the function we are trying to fit to the data.
        '''
        return a * np.exp(-b * x) + c
    # Plot the actual data
    if pltData:
        ax.plot(xdata, ydata, "k.", label=label1);

    # The actual curve fitting happens here
    optimizedParameters, pcov = opt.curve_fit(func, xdata, ydata);

    # Use the optimized parameters to plot the best fit
    ax.plot(xdata, func(xdata, *optimizedParameters),'k-',label=label2);

    # Show the graph
    ax.legend()
    if fig is not None:
        return fig,ax
    else:
        return None
def movingaverage (values, window):
    '''
    Return a new array with the moving averaged values of the given window.
    values: The array of the values
    window: moving average window
    '''
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma
def movingaverage2D(values,window):
    '''
    Return a new array with the moving averaged values of the given window.
    values: The 2D array of the values
    window: (intiger) Dimension of a side of the moving average window. (ex. 3 for 3by3 moving average)
    '''    
    weights = np.ones((window,window),dtype=float)/window/window
    sma = signal.convolve2d(values, weights, mode='same',boundary='wrap')
    return sma
    
def find_CDF(data,bins=None):
    '''
    To find cumulative distribution function (CDF)
    data: 1D data array
    bins: bins for the histogram
    '''
    data=data[~np.isnan(data)]#removing nans
    if np.sum(np.isnan(data))>0:
        print('Nan s ingnored!!')
    weights=np.ones_like(data)/len(data)
    if bins==None:
        val, base = np.histogram(data, weights=weights)
    else:
        val, base = np.histogram(data, bins=bins,weights=weights)
    x=base[1:]
    cdf=np.cumsum(val)
    return x, cdf/cdf.max()

def rmvxtrms(dM):
    '''
    remove extreim value from a data set
    '''
    if dM.size!=0:
        q3=np.percentile(dM,75);q1=np.percentile(dM,25)
        dMmin=q1-1.5*(q3-q1);dMmax=q3+1.5*(q3-q1)
        dM=dM[dM>dMmin];dM=dM[dM<dMmax]
#    else:
#        print('Zero size dM!!!')
    return dM
    
def idx_closest(arr,val):
    '''
    Returns the clossest index of arr to value val.
    '''
    return abs(arr-val).argmin()
def write_latex_table(X,table_name,colnames=None,rownames=None,stl='%0.2f',row_separation=4):
    '''
    To write arrays as latex tables
    colnames should have ALL column names, including the LABEL for ROWs.
    '''
    if colnames is None:
        colnames=['col%d'%(i) for i in range(X.shape[1])]
    if rownames is None:
        rownames=['row%d'%(i) for i in range(X.shape[0])]
    colnames=[r.replace('_','\_') for r in colnames]
    rownames=[r.replace('_','\_') for r in rownames]
    cols=colnames[:-1]
    f1=open(table_name,'w')
    [f1.write(i+' & ') for i in cols]
    f1.write(colnames[-1]+' \\\\'+'\n\hline\n')
    row_count=0
    for i in range(np.size(rownames)):
        f1.write(rownames[i])
        row_count+=1
        for j in range(np.size(colnames[1:])):
            f1.write('&'+stl%(X[i,j]))
        f1.write(' \\\\'+'\n')
        if row_count%row_separation==0:
            f1.write('\hline\n')
    f1.write('\hline\n')
    f1.close()
'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
def setup_figures(plt):
    '''
        xtick.labelsize: 16
        ytick.labelsize: 16
        font.size: 15
        figure.autolayout: True
        figure.figsize: 7.2,4.45
        axes.titlesize : 16
        axes.labelsize : 17
        lines.linewidth : 2
        lines.markersize : 6
        legend.fontsize: 13
        mathtext.fontset: stix
        font.family: STIXGeneral
        Ref:-http://aeturrell.com/2018/01/31/publication-quality-plots-in-python/
    '''
#    plt.rc('font', family='serif',weight='bold',size=16)
    plt.rc('font', family='STIXGeneral',weight='bold',size=15)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.rc('lines',linewidth=2.0,markersize=6)
    plt.rc('legend',fontsize=14,frameon=False)
    plt.rc('axes',titlesize=16,labelsize=17)
    plt.rc('mathtext',fontset='stix')
    
def sub_labels(ax,clr='k',x=0.01,y=0.9,size=15):
    '''
    Add subfigure labels
    '''
    axs=ax.flat
    for n, ax in enumerate(axs):
        ax.text(x, y, '('+string.ascii_lowercase[n]+')', transform=ax.transAxes,
                size=size,color=clr)
        
def add_common_cb(fig,ctf,ts=None,label=None):
    '''
    To add a common color bar for a subplot
    fig: Figure object matplotlib.figure.Figure
    ctf:contourf output instance
    '''
    fig.subplots_adjust(right=0.81)
    cbar_ax = fig.add_axes([0.84, 0.15, 0.01, 0.7])
    if ts is not None and label is not None:
        fig.colorbar(ctf, cax=cbar_ax, ticks=ts,label=label)
    elif label is not None:
        fig.colorbar(ctf, cax=cbar_ax,label=label)
    else:
        fig.colorbar(ctf,cax=cbar_ax)
def add_cb(fig,ctf,ax,ticks=None,orientation='horizontal',label='label',pad=0.2):
    divider = make_axes_locatable(ax)
    if orientation=='horizontal':
        cax = divider.append_axes("bottom", size="5%", pad=pad)
    elif orientation=='vertical':
        cax = divider.append_axes("right", size="5%", pad=pad)
    if ticks is not None:
        fig.colorbar(ctf, cax=cax,ticks=ticks,orientation=orientation,label=label)
    else:
        fig.colorbar(ctf, cax=cax,orientation=orientation,label=label)
    
def savefig(fig,fig_ttl,path=None,rp=False):
    '''
    fig: Figure object matplotlib.figure.Figure
    fig_ttl: figure title string (some specila characterls will be removed from the file name)
    '''
    for ch in [' ','[',']','=']:
        if ch in fig_ttl:
            fig_ttl=fig_ttl.replace(ch,'_')
    fig_ttl=fig_ttl.replace('.','p')
    if path==None:
        filename=fig_ttl
    else:
        filename=path+fig_ttl
    if not(rp):
        if os.path.isfile(filename+'.png'):
            usr=input('Replace existing file?: ')
            if usr=='y':
                rp=True
        else:
            rp=True
    if rp:
        fig.savefig(filename+'.png', format='png', dpi=200)
        print(filename+'.png SAVED.')

def add_subplot_axes(ax,rect,axisbg='w'):
    '''
    example1():
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    rect = [0.2,0.2,0.7,0.7]
    ax1 = add_subplot_axes(ax,rect)
    ax2 = add_subplot_axes(ax1,rect)
    ax3 = add_subplot_axes(ax2,rect)
    plt.show()
    '''
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

def Corr_plot_frames(ax,lineS='b-',linewidth=1.0):
    '''
    Draws frames for a correlation plot. 
    '''
    xb=ax.get_xbound()
    yb=ax.get_ybound()
    xmin,xmax=float(xb[0]),float(xb[1])
    ymin,ymax=float(yb[0]),float(yb[1])
    bmx=np.max([xmax,ymax])
    bmn=np.min([xmin,ymin])    
    ax.plot([bmn,bmx],[bmn,bmx],lineS,linewidth=linewidth)
    ax.plot([0.0,0.0],[bmn,bmx],lineS,linewidth=linewidth)
    ax.plot([bmn,bmx],[0.0,0.0],lineS,linewidth=linewidth)
    ax.set_xlim(bmn,bmx)
    ax.set_ylim(bmn,bmx)

def mapProject(ax,lat={'mn':-90,'mx':90},lon={'mn':-180,'mx':180},line_step={'lon':45.0,'lat':30.0}):
    lat_int=line_step['lat']
    lon_int=line_step['lon']

    mapproj = bm.Basemap(ax=ax,projection='cyl',llcrnrlat= lat['mn'], llcrnrlon= lon['mn'],urcrnrlat= lat['mx'], urcrnrlon= lon['mx'])
    
    latlines = np.arange(lat['mn'],lat['mx'],lat_int)
    lonlines = np.arange(lon['mn'],lon['mx'],lon_int)
    
    mapproj.drawcoastlines()
    mapproj.drawparallels(latlines, labels=[1,0,0,0])
    mapproj.drawmeridians(lonlines, labels=[0,0,0,1])
    
    return mapproj

def array3D_slice(axis,array,x,y,z,slice_i,fig,ax,mnmx=None,cb=True,cb_label='label'):  
    '''
    To plot slices of 3D volumetric data array
    axis: 'x' / 'y' / 'z'
    array:arr[x,y,z]  # parameter
    x,y,z:arr,arr,arr # spacial coordinates
    slice_i: slice position index
    mnmx={vmin:{0},vmax:{1}}
    '''
    def x_slice(array,xslice_i,x,y,z):
        x_cut = array[xslice_i,:,:]
        Y, Z = np.meshgrid(y,z)
        Y, Z = Y.T, Z.T
        X = x[xslice_i] * np.ones((n_y, n_z))
        return x_cut,X,Y,Z    
    def y_slice(array,yslice_i,x,y,z):
        y_cut = array[:,yslice_i,:]
        X, Z = np.meshgrid(x,z)
        X, Z = X.T, Z.T
        Y = y[yslice_i] * np.ones((n_x, n_z))
        return y_cut,X,Y,Z
    def z_slice(array,zslice_i,x,y,z):
        z_cut = array[:,:,zslice_i]
        X, Y = np.meshgrid(x,y)
        X, Y = X.T, Y.T
        Z = z[zslice_i] * np.ones((n_x, n_y))
        return z_cut,X,Y,Z        
    if mnmx is not None:
        min_val=mnmx['vmin']
        max_val=mnmx['vmax']
    else:        
        min_val = array.min()
        max_val = array.max()
    n_x, n_y, n_z = array.shape
    colormap = plt.cm.jet
    if axis=='x':
        cut,X,Y,Z=x_slice(array,slice_i,x,y,z)
    elif axis=='y':
        cut,X,Y,Z=y_slice(array,slice_i,x,y,z)
    elif axis=='z':
        cut,X,Y,Z=z_slice(array,slice_i,x,y,z)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=colormap((cut-min_val)/(max_val-min_val)), shade=False)    
    norm = mp.colors.Normalize(vmin=min_val, vmax=max_val)
    m = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=norm)
    m.set_array([])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    if cb:
        fig.colorbar(m,ax=ax,label=cb_label)
    fig.show() 
def axBorderOff(ax,sides=['top']):
    '''
    To remove box without removing axis labels.
    '''
    for s in sides:
        ax.spines[s].set_visible(False)
def multiple_ticks(ax_XY,major_step,minor_step=None,major_format='%0.1f'):
    '''
    ax_XY: matplotlib.axis.YAxis or matplotlib.axis.XAxis
            ax[0,0].yaxis or ax[0,0].xaxis
    '''
    ax_XY.set_major_locator(MultipleLocator(major_step))
    ax_XY.set_major_formatter(FormatStrFormatter(major_format))
    if minor_step is not None:
        ax_XY.set_minor_locator(MultipleLocator(minor_step))
    else:
        ax_XY.set_minor_locator(AutoMinorLocator())

def shade_quadrant(ax,q=[1,0,0,0]):
    '''
    Use after plotting all the data
    q=[I,II,III,IV] quadrant. 1 to mask 0 to not
    '''
    arr=np.array([[q[1],q[0]],[q[2],q[3]]])
    x0,x1=ax.get_xlim()
    y0,y1=ax.get_ylim()
    vmax = np.abs(np.array([x0,x1,y0,y1])).max() 
    ax.autoscale(False)
    ax.imshow(arr,extent=[vmax*-1,vmax, vmax*-1,vmax],cmap=plt.cm.Greys, interpolation='none', alpha=.2)
