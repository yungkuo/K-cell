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

filePath = '/Users/yungkuo/Documents/Data/K+ cell/081415 TMRM HEK K potential/'

def listdir_nohidden(path):
    for f in os.listdir(path):
        if f.startswith('004'):
            yield f

I = []
t = []
fig2, axarr = plt.subplots(1, 14, figsize=(14,1))
for count, file in enumerate(listdir_nohidden(filePath)):
    current_file = os.path.join(filePath, file)
    mov = libtiff.TiffFile(current_file)
    mov = np.array(mov.get_tiff_array()[:,:,:], dtype='d')
    mov = np.squeeze(mov)
    S = file.split('_')[2]
    t = np.append(t, S)

    if count == 0:
        #frame=len(mov[:,0,0])
        nrow=len(mov[:,0])
        ncol=len(mov[0,:])
        fig1, ax = plt.subplots()
        cax = ax.imshow(mov, interpolation='none', cmap='afmhot')
        cbar = fig1.colorbar(cax)
        print("Choose ROI by clicking multiple points")
        print("If finished, press enter")
        pts = plt.ginput(0,0)
        pts = np.array(pts)
        ax.plot(pts[:,0],pts[:,1], '-o')
        ax.set_xlim([0,ncol])
        ax.set_ylim([nrow,0])
    I1, img = polygon.mean_polygon(mov, pts)
    I = np.append(I, I1)
    cax = axarr[count].imshow(
          img[0,np.min(pts[:,1]):np.max(pts[:,1]),np.min(pts[:,0]):np.max(pts[:,0])],
          interpolation='none', cmap='afmhot', vmin=0, vmax=10000)
    axarr[count].set_xticks([])
    axarr[count].set_yticks([])
    fig2.canvas.draw()
cbar = fig2.colorbar(cax)
fig, ax = plt.subplots()
ax.plot(t, I,'o')




