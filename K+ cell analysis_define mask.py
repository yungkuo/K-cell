# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 16:12:45 2015

@author: yungkuo
"""

import numpy as np
import matplotlib.pyplot as plt
import tifffile as tff
from sub import polygon
import os
#import matplotlib.animation as animation

filePath = '/Users/yungkuo/Documents/Data/K+ cell/011116 peristaltic pump/'
savemask = 1
mask_type = 2 # 1 = polygon mask, 2 = threshold mask, 3 = multi spot

def listdir_nohidden(path):
    for f in os.listdir(path):
        if f.startswith('Anepps_VA_Period_40_cycle_10_X'):
            yield f
frame = 0
I = []
bg = []
for count, file in enumerate(listdir_nohidden(filePath)):
    current_file = os.path.join(filePath, file)
    tiffimg = tff.TiffFile(current_file)
    data = tiffimg.asarray().shape
    if count == 0:
        image = np.zeros((data[1],data[2]), dtype = int)
        for i in range(data[0]/500):
            image = image + tiffimg[i].asarray()
        fig, ax = plt.subplots(figsize=(10,10))
        cax = ax.imshow(image, cmap='afmhot')
        cbar = fig.colorbar(cax, ax=ax)
        if mask_type == 1:
            pts = plt.ginput(0, timeout=0)
            pts = np.array(pts)
            mask = polygon.mask2d_polygon(image,pts)
        elif mask_type == 2:
            fig, ax = plt.subplots(2,1)
            ax[0].hist(np.reshape(image[:,:],[data[1]*data[2]]), bins=1000, histtype='step')
            threshold = plt.ginput(1)[0][0]
            ax[0].vlines(x= threshold, ymin= 0, ymax= 1000, color='r')
            mask = np.zeros((data[1],data[2]), dtype=int)
            for i in range(data[1]):
                for j in range(data[2]):
                    if image[i,j] > threshold:
                        mask[i,j] = 1
                    else:
                        mask[i,j] = 0
            maskedimg = np.multiply(image,mask)
            cax = ax[1].imshow(maskedimg, cmap = 'afmhot')
            cbar = fig.colorbar(cax, ax=ax[1])
        elif mask_type == 3:
            mask = np.zeros((data[1],data[2]))
            scan = 5
            pts = plt.ginput(0, timeout=0)
            pts = np.array(pts)
            for i in range(len(pts)):
                mask[pts[i,1]-scan:pts[i,1]+scan,pts[i,0]-scan:pts[i,0]+scan] = 1
            cax = ax.imshow(mask, alpha=0.5,cmap='gray')
            cbar = fig.colorbar(cax, ax=ax)
            ax.set_xlim(0, data[1])
            ax.set_ylim(data[2], 0)
            fig.canvas.draw()

    for i in range(data[0]):
        if i%2000 == 0:
            if mask_type == 2:
                image = tiffimg[i].asarray()
                fig, ax = plt.subplots(figsize=(10,10))
                cax = ax.imshow(image, cmap='afmhot')
                cbar = fig.colorbar(cax, ax=ax)
                cax1 = ax.imshow(mask, interpolation='none', cmap='gray', alpha=0.5)
                cbar1 = fig.colorbar(cax1, ax=ax)
            elif mask_type == 2:
                image = tiffimg[i].asarray()
                fig, ax = plt.subplots()
                cax = ax.imshow(image[pts[:,1].min():pts[:,1].max(),pts[:,0].min():pts[:,0].max()], interpolation='none', cmap='afmhot')
                cbar = fig.colorbar(cax, ax=ax)
                cax1 = ax.imshow(mask[pts[:,1].min():pts[:,1].max(),pts[:,0].min():pts[:,0].max()], interpolation='none', cmap='gray', alpha=0.5)
                cbar1 = fig.colorbar(cax1, ax=ax)
            elif mask_type == 3:
                image = tiffimg[i].asarray()
                fig, ax = plt.subplots()
                cax = ax.imshow(image[pts[:,1].min()-scan:pts[:,1].max()+scan,pts[:,0].min()-scan:pts[:,0].max()+scan], interpolation='none', cmap='afmhot')
                cbar = fig.colorbar(cax, ax=ax)
                cax1 = ax.imshow(mask[pts[:,1].min()-scan:pts[:,1].max()+scan,pts[:,0].min()-scan:pts[:,0].max()+scan], interpolation='none', cmap='gray', alpha=0.5)
                cbar1 = fig.colorbar(cax1, ax=ax)
#%%
if savemask == 1:
    np.save(filePath+'_mask{}.npy'.format(mask_type), mask)
