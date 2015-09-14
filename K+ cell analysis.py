# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 16:12:45 2015

@author: yungkuo
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import libtiff
from sub import polygon

savefig = 1
abc = 'c'
squareROI = 1 #1 = squareROI with area: (scan*2)*(scan*2); 0 = polygon ROI
scan = 30
filePath = '/Users/yungkuo/Documents/Data/K+ cell/091215 TMRM HEK K potential/TMRM in DMEM_250nM_20min__4/Pos0/'
param = open(filePath+'metadata.txt', 'r')
for count, line in enumerate(param.readlines()):
    if count == 4:
        print line
        interval = float(line.split(':')[1].split(',')[0])
    if count == 39:
        print line
        frame = float(line.split(':')[1].split(',')[0])

def listdir_nohidden(path):
    for f in os.listdir(path):
        if f.endswith('_Green_000.tif'):
            yield f

I = []
t = []
bg = []
a = int(np.sqrt(frame))
fig2, axarr = plt.subplots(a+1, a, figsize=(10,10))
for count, file in enumerate(listdir_nohidden(filePath)):
    current_file = os.path.join(filePath, file)
    mov = libtiff.TiffFile(current_file)
    mov = np.array(mov.get_tiff_array()[:,:,:], dtype='d')
    mov = np.squeeze(mov)
    S = file.split('_')[1]
    t = np.append(t, S)
    if squareROI == 1:
        if count == 0:
            #frame=len(mov[:,0,0])
            nrow=len(mov[:,0])
            ncol=len(mov[0,:])
            fig1, ax = plt.subplots()
            cax = ax.imshow(mov, interpolation='none', cmap='afmhot')
            cbar = fig1.colorbar(cax)
            print("SQUARE ROI!!!! Click the center of square ROI")
            pts = plt.ginput(1)
            pts = np.array(pts)
            print("click background")
            bgpts = plt.ginput(1)
            bgpts = np.array(bgpts)
            ax.plot([pts[:,0]+scan,pts[:,0]+scan,pts[:,0]-scan,pts[:,0]-scan,pts[:,0]+scan],[pts[:,1]+scan,pts[:,1]-scan,pts[:,1]-scan,pts[:,1]+scan,pts[:,1]+scan], '-o')
            ax.plot([bgpts[:,0]+scan,bgpts[:,0]+scan,bgpts[:,0]-scan,bgpts[:,0]-scan,bgpts[:,0]+scan],[bgpts[:,1]+scan,bgpts[:,1]-scan,bgpts[:,1]-scan,bgpts[:,1]+scan,bgpts[:,1]+scan], '-o', color='w')
            ax.set_xlim([0,ncol])
            ax.set_ylim([nrow,0])
        I1 = np.mean(mov[pts[:,1]-scan:pts[:,1]+scan,pts[:,0]-scan:pts[:,0]+scan])
        bg1 = np.mean(mov[bgpts[:,1]-scan:bgpts[:,1]+scan,bgpts[:,0]-scan:bgpts[:,0]+scan])
        I = np.append(I, I1)
        bg = np.append(bg, bg1)
        cax = axarr[count/a, count%a].imshow(
              mov[pts[:,1]-scan:pts[:,1]+scan,pts[:,0]-scan:pts[:,0]+scan],
              interpolation='none', cmap='afmhot', vmin=0, vmax=10000)
    if squareROI == 0:
        if count == 0:
            #frame=len(mov[:,0,0])
            nrow=len(mov[:,0])
            ncol=len(mov[0,:])
            fig1, ax = plt.subplots()
            cax = ax.imshow(mov, interpolation='none', cmap='afmhot')
            cbar = fig1.colorbar(cax)
            print("Choose polygon ROI by clicking multiple points")
            print("If finished, press enter")
            pts = plt.ginput(0,0)
            pts = np.array(pts)
            ax.plot(pts[:,0],pts[:,1], '-o')
            ax.set_xlim([0,ncol])
            ax.set_ylim([nrow,0])
        I1, img = polygon.mean_polygon(mov, pts)
        I = np.append(I, I1)
        cax = axarr[count/a, count%a].imshow(
              img[np.min(pts[:,1]):np.max(pts[:,1]),np.min(pts[:,0]):np.max(pts[:,0])],
              interpolation='none', cmap='afmhot', vmin=0, vmax=10000)
    axarr[count/a, count%a].set_xticks([])
    axarr[count/a, count%a].set_yticks([])
    fig2.canvas.draw()
t = np.array(t, dtype=float)*interval
cbar = fig2.colorbar(cax)
fig3, ax = plt.subplots()
ax.plot(t, I,'-o')
if squareROI == 1:
    ax.plot(t, I-bg, '-rs')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Fluorescence intensity')
if savefig == 1:
    fig1.savefig(filePath+abc+'.fig1_img.pdf', format='pdf', bbox_inches = 'tight')
    fig2.savefig(filePath+abc+'.fig2_ROI.pdf', format='pdf', bbox_inches = 'tight')
    fig3.savefig(filePath+abc+'.fig3_IT curve.pdf', format='pdf', bbox_inches = 'tight')
