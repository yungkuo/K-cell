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
abc = 'a'
squareROI = 0 #1 = squareROI with area: (scan*2)*(scan*2); 0 = polygon ROI
scan = 20
DiOboundary = 0 #Use DiO channel to help click the coundary of cells
assignpts = 0 #1 = provide pts, 0 = ginput
filePath = '/Users/yungkuo/Documents/Data/K+ cell/091715 Anepps HEK K potential/Anepps_0.4uM_addVA_0.2uM_addK__1/Pos0/'

if DiOboundary == 1:
    squareROI = 0
if assignpts == 1:
    DiOboundary = 1
pts = np.array([[ 128.23610224,   73.51916933],
 [ 132.32555911,   89.87699681],
 [ 138.45974441,  109.30191693],
 [ 152.77284345,  122.59265176],
 [ 169.13067093,  136.9057508 ],
 [ 186.51086262,  145.08466454],
 [ 203.89105431,  147.12939297],
 [ 222.29361022,  136.9057508 ],
 [ 236.60670927,  123.61501597],
 [ 251.94217252,  108.27955272],
 [ 258.07635783,   89.87699681],
 [ 257.05399361,   77.6086262 ],
 [ 246.83035144,   64.31789137],
 [ 232.5172524 ,   54.0942492 ],
 [ 219.22651757,   46.93769968],
 [ 203.89105431,   44.89297125],
 [ 189.57795527,   47.9600639 ],
 [ 174.24249201,   50.00479233],
 [ 157.88466454,   50.00479233],
 [ 139.48210863,   54.0942492 ],
 [ 129.25846645,   54.0942492 ],
 [ 119.03482428,   50.00479233],
 [ 124.14664537,   63.29552716]])


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
if assignpts == 0:
    if DiOboundary == 1:
        mov1 = libtiff.TiffFile(filePath+'img_000000000_Blue_000.tif')
        mov1 = np.array(mov1.get_tiff_array()[:,:,:], dtype='d')
        mov1 = np.squeeze(mov1)
        mov2 = libtiff.TiffFile(filePath+'img_000000000_Green_000.tif')
        mov2 = np.array(mov2.get_tiff_array()[:,:,:], dtype='d')
        mov2 = np.squeeze(mov2)

        fig4, ax = plt.subplots()
        cax1 = ax.imshow(mov1, interpolation='none', cmap='Greens', alpha=0.5)
        cax2 = ax.imshow(mov2, interpolation='none', cmap='Reds', alpha=0.5)
        cbar = fig4.colorbar(cax1)
        cbar = fig4.colorbar(cax2)
        print("choose DiO boundary")
        pts = plt.ginput(0,0)
        pts = np.array(pts)
        print pts

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
              interpolation='none', cmap='afmhot', vmin=0, vmax=20000)
    if DiOboundary == 0:
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
                mask = polygon.mask_polygon(mov, pts)
            I1, img = polygon.mean_polygon(mov, mask)
            I = np.append(I, I1)
            cax = axarr[count/a, count%a].imshow(
                  img[np.min(pts[:,1]):np.max(pts[:,1]),np.min(pts[:,0]):np.max(pts[:,0])],
                  interpolation='none', cmap='afmhot', vmin=0, vmax=20000)
        axarr[count/a, count%a].set_xticks([])
        axarr[count/a, count%a].set_yticks([])
        fig2.canvas.draw()
    if DiOboundary == 1:
        if count == 0:
            #frame=len(mov[:,0,0])
            nrow=len(mov[:,0])
            ncol=len(mov[0,:])
            fig1, ax = plt.subplots()
            cax = ax.imshow(mov, interpolation='none', cmap='afmhot')
            cbar = fig1.colorbar(cax)
            ax.plot(pts[:,0],pts[:,1], '-o')
            ax.set_xlim([0,ncol])
            ax.set_ylim([nrow,0])
            mask = polygon.mask_polygon(mov, pts)
        I1, img = polygon.mean_polygon(mov, mask)
        I = np.append(I, I1)
        cax = axarr[count/a, count%a].imshow(
            img[np.min(pts[:,1]):np.max(pts[:,1]),np.min(pts[:,0]):np.max(pts[:,0])],
            interpolation='none', cmap='afmhot', vmin=0, vmax=20000)
        axarr[count/a, count%a].set_xticks([])
        axarr[count/a, count%a].set_yticks([])
        fig2.canvas.draw()


t = np.array(t, dtype=float)*interval
cbar = fig2.colorbar(cax)
fig3, ax = plt.subplots()
ax.plot(t, I,'-o')
if squareROI == 1:
    ax.plot(t, I-bg, '-rs')
    ax.plot(t, bg, '-g^')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Fluorescence intensity')

np.save(filePath+'IT_result.npy', I)
if savefig == 1:
    fig1.savefig(filePath+abc+'.fig1_img.pdf', format='pdf', bbox_inches = 'tight')
    fig2.savefig(filePath+abc+'.fig2_ROI.pdf', format='pdf', bbox_inches = 'tight')
    fig3.savefig(filePath+abc+'.fig3_IT curve.pdf', format='pdf', bbox_inches = 'tight')
