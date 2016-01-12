# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 16:12:45 2015

@author: yungkuo
"""

import numpy as np
import matplotlib.pyplot as plt
import libtiff
from sub import polygon
import os
#import matplotlib.animation as animation

filePath = '/Users/yungkuo/Documents/Data/K+ cell/flow chamber/101815 flow cell/ZNR/'
savefig = 1
abc = 'a'
dt = 1/32.352
bgscan = 30
assignmask = 0
thresholdmask = 0
if assignmask == 1:
    if thresholdmask == 1:
        mask = np.load(filePath+'_threshold_mask.npy')
    elif thresholdmask == 0:
        mask = np.load(filePath+'_mask.npy')

def listdir_nohidden(path):
    for f in os.listdir(path):
        if f.startswith('ZNR_low K_60s_high K_60s_X'):
            yield f
frame = 0
I = []
bg = []
for count, file in enumerate(listdir_nohidden(filePath)):
    current_file = os.path.join(filePath, file)
    mov = libtiff.TiffFile(current_file)
    mov = mov.get_tiff_array()
    mov = np.array(mov[:,:,:],dtype='d')
    frame1 = len(mov[:,0,0])
    nrow = len(mov[0,:,0])
    ncol = len(mov[0,0,:])
    if count == 0:
        if assignmask == 0:
            fig1, ax = plt.subplots()
            cax = ax.imshow(mov[0,:,:], interpolation='none', cmap='afmhot')
            cbar = fig1.colorbar(cax, ax=ax)
            print("Choose polygon ROI by clicking multiple points")
            print("If finished, press enter")
            pts = plt.ginput(0, timeout=0)
            pts = np.array(pts)
            mask = polygon.mask_polygon(mov,pts)
            print("Choose background")
            bgpts = plt.ginput(1, timeout=0)
            bgpts = np.array(bgpts)
            ax.plot(pts[:,0],pts[:,1], 'b-+')
            ax.plot([bgpts[:,0]+bgscan,bgpts[:,0]+bgscan,bgpts[:,0]-bgscan,bgpts[:,0]-bgscan,bgpts[:,0]+bgscan],[bgpts[:,1]+bgscan,bgpts[:,1]-bgscan,bgpts[:,1]-bgscan,bgpts[:,1]+bgscan,bgpts[:,1]+bgscan], '-o', color='w')
            ax.set_xlim([0,ncol])
            ax.set_ylim([nrow,0])
            fig1.canvas.draw()
        elif assignmask == 1:
            fig1, ax = plt.subplots()
            cax = ax.imshow(mov[0,:,:], interpolation='none', cmap='afmhot')
            cbar = fig1.colorbar(cax, ax=ax)
            cax1 = ax.imshow(mask, interpolation='none', cmap='gray', alpha=0.5)
            cbar1 = fig1.colorbar(cax1, ax=ax)
            if thresholdmask == 0:
                print("Choose background")
                bgpts = plt.ginput(1, timeout=0)
                bgpts = np.array(bgpts)
                ax.plot([bgpts[:,0]+bgscan,bgpts[:,0]+bgscan,bgpts[:,0]-bgscan,bgpts[:,0]-bgscan,bgpts[:,0]+bgscan],[bgpts[:,1]+bgscan,bgpts[:,1]-bgscan,bgpts[:,1]-bgscan,bgpts[:,1]+bgscan,bgpts[:,1]+bgscan], '-o', color='w')
                ax.set_xlim([0,ncol])
                ax.set_ylim([nrow,0])
                fig1.canvas.draw()
    maskedimg = np.multiply(mov,np.tile(mask, (frame1,1,1)))
    I1 = np.sum(np.sum(maskedimg, axis=1), axis=1)/mask.sum()
    I = np.append(I,I1)
    if thresholdmask == 0:
        bg1 = np.mean(np.mean(mov[:,bgpts[:,1]-bgscan:bgpts[:,1]+bgscan,bgpts[:,0]-bgscan:bgpts[:,0]+bgscan], axis=1), axis=1)
    elif thresholdmask == 1:
        tsum = np.sum(np.sum(maskedimg, axis=1), axis=1)
        tmean = tsum/mask.sum()
        bg1 = (np.sum(np.sum(mov, axis=1), axis=1)-tsum)/(nrow*ncol-mask.sum())
    bg = np.append(bg, bg1)
    frame = frame + frame1
t = np.arange(0,frame*dt, dt)
fig2, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(t, I, 'b', label='I')
ax[0].plot(t, bg, 'g', label='bg')
ax[0].plot(t, I-bg, 'r', label= 'I-bg')
ax[1].plot(t, I-bg, 'r', label= 'I-bg')
ax[0].set_title('Polygon manual select membrane')
ax[1].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
'''
fig3, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(t, tmean, label='I')
ax[0].plot(t, tbg, label='bg')
ax[0].plot(t, tmean-tbg, label= 'I-bg')
ax[1].plot(t, tmean-tbg, label= 'I-bg')
ax[0].set_title('Threshold mask')
ax[1].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
'''
'''
fig3, ax = plt.subplots()
ims = []
for i in range(frame):
    if i%100 == 0:
        im=ax.imshow(mov[i,np.min(pts[:,1]):np.max(pts[:,1]),np.min(pts[:,0]):np.max(pts[:,0])],vmin=mov[:,:,:].min(),vmax=mov[:,:,:].max(),cmap='hot')
        ims.append([im])
ani = animation.ArtistAnimation(fig3, ims, interval=dt*1000, blit=True,repeat_delay=1000)
'''
'''

#fig0, ax = plt.subplots(2,1)
#ax[0].hist(np.reshape(mov[0,:,:],[nrow*ncol]), bins=1000, histtype='step')
#threshold = plt.ginput(1)[0][0]
#ax[0].vlines(x= threshold, ymin= 0, ymax= 1000, color='r')
#tmask = np.zeros((nrow,ncol), dtype=int)
#for i in range(nrow):
#    for j in range(ncol):
#        if mov[0,i,j] > threshold:
#            tmask[i,j] = 1
#        else:
#            tmask[i,j] = 0
#tmask3d = np.tile(tmask, (frame,1,1))
#tmaskedimg = np.multiply(mov,tmask3d)
#tsum = np.sum(np.sum(tmaskedimg, axis=1), axis=1)
#tmean = tsum/tmask.sum()
#tbg = (np.sum(np.sum(mov, axis=1), axis=1)-tsum)/(nrow*ncol-tmask.sum())
#cax = ax[1].imshow(tmaskedimg[0,:,:], interpolation='none', cmap='afmhot')
#cbar = fig0.colorbar(cax, ax=ax[1])
#fig0.canvas.draw()
'''

if savefig == 1:
    np.save(filePath+'I(thresholdmask={}).npy'.format(thresholdmask), I)
    np.save(filePath+'bg(thresholdmask={}).npy'.format(thresholdmask), bg)
    #fig0.savefig(filePath+fileName+'_'+abc+'_fig0_hist.pdf', format='pdf', bbox_inches = 'tight')
    fig1.savefig(filePath+'fig1_img(thresholdmask={}).pdf'.format(thresholdmask), format='pdf', bbox_inches = 'tight')
    fig2.savefig(filePath+'fig2_fluor(thresholdmask={}).pdf'.format(thresholdmask), format='pdf', bbox_inches = 'tight')
