# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 10:05:19 2016

@author: yungkuo
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import tifffile as tff
import lmfit
#import math
#import scipy.ndimage as ndi
#from sub import polygon
#from PIL import Image
#import matplotlib.animation as animation

filePath = '/Users/yungkuo/Google Drive/011716 flow anepps/'
fileName = 'QD_VA_Period_20_cycle_7_1'
maskName = 'QD_VA_Period_20_cycle_7_1_mask3.npy'
savefig = 1
abc = 'ptmask'
dt = 0.03028 # 30ms = 32.352Hz, 100ms =
bgscan = 30 #pixels
bgpoint = [196, 349] #coordinate of the left bottom corner of the bg box
period = 20  #seconds
click_transi = 1 #define transition time points as the max of the background curve's derivative
transiperiod = [100,200] # discard (transiperiod) frames

#%%
mask = np.load(filePath+'results/'+maskName)
datapath = filePath+'raw data/'
def listdir_nohidden(path):
    for f in os.listdir(path):
        if f.startswith(fileName):
            yield f

frame = 0
I = []
bg = []
for count, file in enumerate(listdir_nohidden(datapath)):
    current_file = os.path.join(datapath, file)
    if count == 0:
        tiffimg = tff.TiffFile(current_file)
        data = tiffimg.asarray().shape
        bgmask = np.zeros((data[1],data[2]))
        bgmask[(bgpoint[1]-bgscan):bgpoint[1], bgpoint[0]:(bgpoint[0]+bgscan)] = 1
        image = np.zeros((data[1],data[2]), dtype = int)
        for i in range(data[0]/100):
            image = image + tiffimg[i].asarray()
        fig1, ax = plt.subplots(1,2)
        cax = ax[0].imshow(image, cmap='afmhot')
        #divider = make_axes_locatable(ax[0])
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig1.colorbar(cax, ax=ax[0], fraction=0.05, pad=0.05)
        ax[0].imshow(bgmask, interpolation='none', cmap='gray', alpha=0.2)
        cax1 = ax[1].imshow(mask, interpolation='none', cmap='gray')
        cbar1 = fig1.colorbar(cax1, ax=ax[1], fraction=0.05, pad=0.05)
        ax[0].set_xlim(0, data[1])
        ax[0].set_ylim(data[2], 0)
        fig1.canvas.draw()
    tiffimg = tff.TiffFile(current_file)
    data = tiffimg.asarray().shape

    for i in range(data[0]):
        img = tiffimg[i].asarray()
        I1 = np.sum(img*mask)/np.sum(mask)
        bg1 = np.sum(img*bgmask,dtype='d')/np.sum(bgmask) #Anepps: [(512-bgscan):, :bgscan] #QD [63:121, 180:280]
        I = np.append(I,I1)
        bg = np.append(bg, bg1)
        if i%20000 == 0:
            fig0, ax = plt.subplots()
            cax = ax.imshow(img, interpolation='none', cmap='afmhot')
            cbar = fig0.colorbar(cax, ax=ax)
            ax.imshow(bgmask, interpolation='none', cmap='gray', alpha=0.2)
            cax1 = ax.imshow(mask, interpolation='none', cmap='gray', alpha=0.2)
            cbar1 = fig0.colorbar(cax1, ax=ax)
    frame = frame + data[0]
t = np.arange(0,frame*dt, dt)
fig2, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(t, I, 'b', label='I')
ax[0].plot(t, bg, 'g', label='bg')
ax[0].plot(t, I-bg, 'r', label= 'I-bg')
ax[1].plot(t, I-bg, 'r', label= 'I-bg')
ax[0].set_title('Masked ensemble average')
ax[1].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)

I = I-bg
#%%
def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
def localmax1D(trace, pts, scan):
    x = []
    for i in range(len(pts)):
        localmax = np.max(trace[(pts[i]-scan):(pts[i]+scan)])
        x1 = np.where(trace == localmax)
        x = np.append(x, x1)
    return x
window_size=100
diff_bg = movingaverage(np.gradient(movingaverage(bg, window_size)), window_size)
x = np.arange(0,frame,1)
fig3, ax = plt.subplots(2, sharex=True)
ax[0].plot(x, I, 'c')
ax[0].plot(x, movingaverage(I, window_size), 'r', alpha=0.5)
ax[1].plot(x, np.gradient(movingaverage(I, window_size)), 'b')
ax[1].plot(x, diff_bg, 'y')
ax[1].set_xlim(0, x.max())
if click_transi == 1:
    print 'click to select transition period'
    pts =  plt.ginput(0, timeout=0)
    pts = np.array(pts)
    print pts[:,0]
    pts = localmax1D(diff_bg,pts[:,0],period/dt/4)
    print pts
elif click_transi == 0:
    pts = np.arange(0.5,7,1)*(period/dt*2)
#%%
for i in range(len(pts)):
    ax[0].axvline(x = pts[i], c='0.3', alpha=0.2)
    ax[1].axvline(x = pts[i], c='0.3', alpha=0.2)
fig3.canvas.draw()
#%%
transiperiod = [100,250]
Lweight = np.zeros((frame))
Hweight = np.zeros((frame))
for i in range(len(pts)):
    Hweight[pts[i]+transiperiod[0]:pts[i]+period/dt-transiperiod[0]] = 1
    Lweight[pts[i]-period/dt+transiperiod[1]:pts[i]-transiperiod[0]] = 1
#%%
# fit 2 exponential
def EXP2(x, A1, tau1, A2, tau2):
    return A1*np.exp(-tau1*x)+A2*np.exp(-tau2*x)
mod = lmfit.Model(EXP2)
params = mod.make_params()
params['A1'].set(value = I.max()-I.min(), min=0)
params['tau1'].set(value = -np.log(I[frame-1]/I[0])/t[frame-1])
params['A2'].set(value = I.max()-I.min(), min=0)
params['tau2'].set(value = -np.log(I[frame-1]/I[0])/t[frame-1])
#params['b'].set(value = I.min(), min=0)
result_L = mod.fit(I, x=t, weights=Lweight, **params)
residual_L = I - result_L.best_fit
result_H = mod.fit(I, x=t, weights=Hweight, **params)
residual_H = I - result_H.best_fit

fig4, ax = plt.subplots(2, sharex=True)
ax[0].plot(t, I, 'r-', label='Raw Data')
ax[0].plot(t, result_L.best_fit, '-', label='Lfit', color='k')
ax[0].plot(t, result_H.best_fit, '-', label='Hfit', color='k')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].set_xlim(0, t.max())
ax[1].plot(t, residual_L*Lweight, 'c-', alpha=0.5, label='L')
ax[1].plot(t, residual_H*Hweight, 'm-', alpha=0.5, label='H')
ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
for i in range(len(pts)):
    ax[0].axvline(x = t[pts[i]], c='0.3', alpha=0.2)
    ax[1].axvline(x = t[pts[i]], c='0.3', alpha=0.2)
ax[1].plot(t,Hweight*5,'r')
ax[1].plot(t,Lweight*-5,'b')
ax[1].set_ylim(-10, 10)
fig4.canvas.draw()
#%%
dFF_HLfit = (result_H.best_fit - result_L.best_fit)/result_L.best_fit
dFF_Lfit = residual_L*Hweight/result_L.best_fit
dFF_Hfit = residual_H*Lweight/result_L.best_fit
fig5, ax = plt.subplots()
ax.plot(t, dFF_Lfit, 'r', label='L')
ax.plot(t, -dFF_Hfit, 'g', label='H')
ax.plot(t, dFF_HLfit, 'b', label='HL')
ax.set_ylim(-0.5,1)
ax.legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
#%%
if savefig == 1:
    fig1.savefig(filePath+'results/'+fileName+abc+'_fig1_img.png', format='png', bbox_inches = 'tight')
    fig2.savefig(filePath+'results/'+fileName+abc+'_fig2_fluor.png', format='png', bbox_inches = 'tight')
    fig3.savefig(filePath+'results/'+fileName+abc+'_fig3_transiperiod.png', format='png', bbox_inches = 'tight')
    fig4.savefig(filePath+'results/'+fileName+abc+'_fig4_fit and residuals.png', format='png', bbox_inches = 'tight')
    fig5.savefig(filePath+'results/'+fileName+abc+'_fig5_dFF.png', format='png', bbox_inches = 'tight')
