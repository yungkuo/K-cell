# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 16:12:45 2015

@author: yungkuo
"""

import numpy as np
import matplotlib.pyplot as plt
import libtiff
from sub import polygon
#import matplotlib.animation as animation

filePath = '/Users/yungkuo/Documents/Data/K+ cell/091715 DiO DPA K potential/'
fileName = '002_add high K (28.1, 158)'
savefig = 0
abc = 'a'
dt = 1/32.352
scan = 10
assignpts = 1
if assignpts == 1:
    pts = np.load(filePath+fileName+'_'+abc+'_'+'pts.npy')

mov = libtiff.TiffFile(filePath+fileName+'.tif')
mov = np.array(mov.get_tiff_array()[:,:,:], dtype='d')
#mov = mov[:,230:360, 150:320]
frame=len(mov[:,0,0])
nrow=len(mov[0,:,0])
ncol=len(mov[0,0,:])
t = np.arange(0,frame*dt, dt)

fig0, ax = plt.subplots(2,1)
ax[0].hist(np.reshape(mov[0,:,:],[nrow*ncol]), bins=1000, histtype='step')
threshold = plt.ginput(1)[0][0]
ax[0].vlines(x= threshold, ymin= 0, ymax= 1000, color='r')
tmask = np.zeros((nrow,ncol), dtype=int)
for i in range(nrow):
    for j in range(ncol):
        if mov[0,i,j] > threshold:
            tmask[i,j] = 1
        else:
            tmask[i,j] = 0
tmask3d = np.tile(tmask, (frame,1,1))
tmaskedimg = np.multiply(mov,tmask3d)
tsum = np.sum(np.sum(tmaskedimg, axis=1), axis=1)
tmean = tsum/tmask.sum()
tbg = (np.sum(np.sum(mov, axis=1), axis=1)-tsum)/(nrow*ncol-tmask.sum())
cax = ax[1].imshow(tmaskedimg[0,:,:], interpolation='none', cmap='afmhot')
cbar = fig0.colorbar(cax, ax=ax[1])
fig0.canvas.draw()

fig1, ax = plt.subplots()
cax = ax.imshow(mov[0,:,:], interpolation='none', cmap='afmhot')
cbar = fig1.colorbar(cax, ax=ax)
if assignpts == 0:
    print("Choose polygon ROI by clicking multiple points")
    print("If finished, press enter")
    pts = plt.ginput(0,0)
    pts = np.array(pts)
print("Choose background")
bgpts = plt.ginput(1)
bgpts = np.array(bgpts)
ax.plot(pts[:,0],pts[:,1], 'b-+')
ax.plot([bgpts[:,0]+scan,bgpts[:,0]+scan,bgpts[:,0]-scan,bgpts[:,0]-scan,bgpts[:,0]+scan],[bgpts[:,1]+scan,bgpts[:,1]-scan,bgpts[:,1]-scan,bgpts[:,1]+scan,bgpts[:,1]+scan], '-o', color='w')
ax.set_xlim([0,ncol])
ax.set_ylim([nrow,0])
fig1.canvas.draw()
I = polygon.mean3d_polygon(mov, pts)
bg = np.mean(np.mean(mov[:,bgpts[:,1]-scan:bgpts[:,1]+scan,bgpts[:,0]-scan:bgpts[:,0]+scan], axis=1), axis=1)

fig2, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(t, I, label='I')
ax[0].plot(t, bg, label='bg')
ax[0].plot(t, I-bg, label= 'I-bg')
ax[1].plot(t, I-bg, label= 'I-bg')
ax[0].set_title('Polygon manual select membrane')
ax[1].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)

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
fig3, ax = plt.subplots()
ims = []
for i in range(frame):
    if i%5 ==0:
        im=ax.imshow(mov[i,np.min(pts[:,1]):np.max(pts[:,1]),np.min(pts[:,0]):np.max(pts[:,0])],vmin=mov[:,:,:].min(),vmax=mov[:,:,:].max(),cmap='hot')
        ims.append([im])
ani = animation.ArtistAnimation(fig3, ims, interval=dt*1000, blit=True,repeat_delay=1000)
'''

if savefig == 1:
    np.save(filePath+fileName+'_'+abc+'_'+'pts.npy', pts)
    fig0.savefig(filePath+fileName+'_'+abc+'_fig0_hist.pdf', format='pdf', bbox_inches = 'tight')
    fig1.savefig(filePath+fileName+'_'+abc+'_fig1_img.pdf', format='pdf', bbox_inches = 'tight')
    fig2.savefig(filePath+fileName+'_'+abc+'_fig2_fluor.pdf', format='pdf', bbox_inches = 'tight')
