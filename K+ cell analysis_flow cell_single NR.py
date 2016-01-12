# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 16:12:45 2015

@author: yungkuo
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import tifffile as tff
import lmfit
import scipy.ndimage as ndi
#from PIL import Image
#import matplotlib.animation as animation

filePath = '/Users/yungkuo/Documents/Data/K+ cell/flow chamber/101815 flow cell/ZNR/'
savefig = 0
abc = 'QD0'
dt = 1/32.352
scan = 5
bgscan = 30
assignpts = 1
period = 60
transiperiod = 300 # discard (transiperiod) frames
                  # befroe and after the transition of low/high K buffer
                  # e.g. discard 1941-5 to 1941+5 timetrace points when calculating the mean
if assignpts == 1:
    pts = np.array([[372, 297],
                    [244, 218],
                    [305, 330],
                    [394, 370],
                    [287, 453],
                    [175, 126],
                    [254, 332],
                    [384, 402],
                    [219, 400],
                    [116, 318]])
    pts = pts[[int(s) for s in abc.split('D') if s.isdigit()]]

def listdir_nohidden(path):
    for f in os.listdir(path):
        if f.startswith('ZNR_low K_60s_high K_60s_X'):
            yield f
frame = 0
I = []
bg = []
for count, file in enumerate(listdir_nohidden(filePath)):
    current_file = os.path.join(filePath, file)
    if count == 0:
        tiffimg = tff.TiffFile(current_file)
        data = tiffimg.asarray().shape
        image = np.zeros((data[1],data[2]), dtype = int)
        for i in range(data[0]/100):
            image = image + tiffimg[i].asarray()
        fig1, ax = plt.subplots(figsize=(10,10))
        cax = ax.imshow(image, cmap='afmhot')
        cbar = fig1.colorbar(cax, ax=ax)
        if assignpts == 0:
            pts = plt.ginput(0,0)
            pts = np.array(pts)
        ax.plot((pts[:,0]-scan, pts[:,0]-scan, pts[:,0]+scan, pts[:,0]+scan, pts[:,0]-scan),
                (pts[:,1]+scan, pts[:,1]-scan, pts[:,1]-scan, pts[:,1]+scan, pts[:,1]+scan), '-+', color='b')
        ax.set_xlim(0, 512)
        ax.set_ylim(512, 0)
        fig1.canvas.draw()
    tiffimg = tff.TiffFile(current_file)
    data = tiffimg.asarray().shape
    for i in range(data[0]):
        img = tiffimg[i].asarray()
        I1 = np.mean(img[pts[:,1]-scan:pts[:,1]+scan,pts[:,0]-scan:pts[:,0]+scan])
        bg1 = np.mean(img[:bgscan, 512-bgscan:],dtype='d')
        I = np.append(I,I1)
        bg = np.append(bg, bg1)
        if i%2000 == 0:
            fig0, ax = plt.subplots()
            ax.imshow(img[pts[:,1]-scan:pts[:,1]+scan,pts[:,0]-scan:pts[:,0]+scan], interpolation='none', cmap='gray')
    frame = frame + data[0]
t = np.arange(0,frame*dt, dt)
fig2, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(t, I, 'b', label='I')
ax[0].plot(t, bg, 'g', label='bg')
ax[0].plot(t, I-bg, 'r', label= 'I-bg')
ax[1].plot(t, I-bg, 'r', label= 'I-bg')
ax[0].set_title('Manul select QD')
ax[1].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)

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

if savefig == 1:
    np.save(filePath+abc+'_I.npy', I)
    np.save(filePath+abc+'_bg.npy', bg)
    #fig0.savefig(filePath+'fig0_QD.pdf', format='pdf', bbox_inches = 'tight')
    fig1.savefig(filePath+abc+'_fig1_img.png', format='png', bbox_inches = 'tight')
    fig2.savefig(filePath+abc+'_fig2_fluor.png', format='png', bbox_inches = 'tight')
#%%
# fit 2 exponential
def EXP2(x, A1, tau1, A2, tau2, b):
    return A1*np.exp(-tau1*x)+A2*np.exp(-tau2*x)+b
gmod = lmfit.Model(EXP2)
params = gmod.make_params()
params['A1'].set(value = I.max()-I.min(), min=0)
params['tau1'].set(value = -np.log(I[frame-1]/I[0])/t[frame-1])
params['A2'].set(value = I.max()-I.min(), min=0)
params['tau2'].set(value = -np.log(I[frame-1]/I[0])/t[frame-1])
params['b'].set(value = I.min(), min=0)
result = gmod.fit(I, x=t, **params)
fI = EXP2(t, result.best_values['A1'], result.best_values['tau1'],result.best_values['A2'], result.best_values['tau2'], result.best_values['b'])
fig1, ax = plt.subplots(2)
#result.plot_fit(ax=ax[0], datafmt='-', fitfmt='o-', initfmt='--', numpoints=40)
result.plot_residuals(ax=ax[1], datafmt='c-')
ax[0].plot(t, I, 'r-', label='Raw Data', alpha=0.5)
ax[0].plot(t, fI, '-', label='fit', color='k')
ax[0].plot(t, EXP2(t, I.max()-I.min(), -np.log(I[frame-1]/I[0])/t[frame-1],I.max()-I.min(), -np.log(I[frame-1]/I[0])/t[frame-1], I.min()), '--', label='initial guess', color='b')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].set_xlim(0, t.max())
ax[1].set_xlim(0, t.max())

lk_trace = [(I-fI)[s] for s in range(frame) if s/int(period/dt)%2 == 0 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)]
hk_trace = [(I-fI)[s] for s in range(frame) if s/int(period/dt)%2 == 1 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)]
lk_mean = []
hk_mean = []
dff = []
for i in range(int(t.max()/period/2)+1):
    lk_mean1 = np.mean(lk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    hk_mean1 = np.mean(hk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    dff1 = np.abs(lk_mean1-hk_mean1)/np.mean(fI[i*2*int(period/dt-transiperiod*2):(i+1)*2*int(period/dt-transiperiod*2)])
    lk_mean = np.append(lk_mean, lk_mean1)
    hk_mean = np.append(hk_mean, hk_mean1)
    dff = np.append(dff, dff1)
print np.abs(lk_mean-hk_mean)
print dff

fig2, ax = plt.subplots()
ax.plot(t, (I-fI), 'c.', markersize=1, label='residual', alpha=0.5)
for i in range(int(t.max()/period/2)+1):
    l_tfrag = t[i*2*int(period/dt)+transiperiod:(i*2+1)*int(period/dt)-transiperiod]
    h_tfrag = t[(i*2+1)*int(period/dt)+transiperiod:(i*2+2)*int(period/dt)-transiperiod]
    ax.plot(l_tfrag, np.repeat(lk_mean[i],len(l_tfrag)), 'r-', label='low K mean {}'.format(i)+ '= {}'.format(round(lk_mean[i],4)))
    ax.plot(h_tfrag, np.repeat(hk_mean[i],len(h_tfrag)), 'b-', label='high K mean {}'.format(i)+ '= {}'.format(round(hk_mean[i],4)))
for i in range(int(t.max()/period)+1):
    ax.axvline(x = period*i, c='0.3', alpha=0.2)
ax.legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax.set_xlim(0, t.max())
ax.set_title('fit 2 exponential')
plt.xticks(np.arange(0, t.max(), period))
for i in range(int(t.max()/period/2)):
    ax.annotate(r'$\Delta$F/F={}'.format(round(dff[i],3)), xy=(0,0), xytext=(period*2*i,10), xycoords='axes pixels', fontsize=10)
fig2.canvas.draw()
#%%
# fit polynomial
polyfit = np.polyfit(t, I, 9)
fI = np.polyval(polyfit, t)
fig3, ax = plt.subplots(2)
ax[0].plot(t, I, 'r-', label='Raw Data', alpha=0.5)
ax[0].plot(t, fI, '-', label='fit(polynomial)', color='k')
ax[1].plot(t, I-fI, '-', label='residual', color='c')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].set_xlim(0, t.max())
ax[1].set_xlim(0, t.max())
lk_trace = [(I-fI)[s] for s in range(frame) if s/int(period/dt)%2 == 0 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)]
hk_trace = [(I-fI)[s] for s in range(frame) if s/int(period/dt)%2 == 1 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)]
lk_mean = []
hk_mean = []
dff = []
for i in range(int(t.max()/period/2)+1):
    lk_mean1 = np.mean(lk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    hk_mean1 = np.mean(hk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    dff1 = np.abs(lk_mean1-hk_mean1)/np.mean(fI[int(period/dt)*2*i:int(period/dt)*(2*i+2)])
    lk_mean = np.append(lk_mean, lk_mean1)
    hk_mean = np.append(hk_mean, hk_mean1)
    dff = np.append(dff, dff1)
print np.abs(lk_mean-hk_mean)
print dff

fig4, ax = plt.subplots()
ax.plot(t, (I-fI), 'c.', markersize=1, label='residual', alpha=0.5)
for i in range(int(t.max()/period/2)+1):
    l_tfrag = t[i*2*int(period/dt)+transiperiod:(i*2+1)*int(period/dt)-transiperiod]
    h_tfrag = t[(i*2+1)*int(period/dt)+transiperiod:(i*2+2)*int(period/dt)-transiperiod]
    ax.plot(l_tfrag, np.repeat(lk_mean[i],len(l_tfrag)), 'r-', label='low K mean {}'.format(i)+ '= {}'.format(round(lk_mean[i],4)))
    ax.plot(h_tfrag, np.repeat(hk_mean[i],len(h_tfrag)), 'b-', label='high K mean {}'.format(i)+ '= {}'.format(round(hk_mean[i],4)))
for i in range(int(t.max()/period)+1):
    ax.axvline(x = period*i, c='0.3', alpha=0.2)
ax.legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax.set_xlim(0, t.max())
ax.set_title('fit polynomial')
plt.xticks(np.arange(0, t.max(), period))
for i in range(int(t.max()/period/2)):
    ax.annotate(r'$\Delta$F/F={}'.format(round(dff[i],3)), xy=(0,0), xytext=(period*2*i,10), xycoords='axes pixels', fontsize=10)
fig4.canvas.draw()
#%%
# Gaussian filter
'''
#artificial trace for test
I = np.zeros((16506))
for s in range(16506):
    if s/int(period/dt)%2 == 0:
        I[s] = 1
I = I+10+np.random.rand(16506)*10
'''
fI = ndi.filters.gaussian_filter1d(I, 1300)
fig5, ax = plt.subplots(2)
ax[0].plot(t, I, 'r-', label='Raw Data', alpha=0.5)
ax[0].plot(t, fI, '-', label='gaussian filtered', color='g')
ax[1].plot(t, I-fI, '-', label='residual', color='c')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].set_xlim(0, t.max())
ax[1].set_xlabel('Time (s)')
ax[1].set_ylabel('residuals')
ax[1].set_xlim(0, t.max())

lk_trace = [(I-fI)[s] for s in range(frame) if s/int(period/dt)%2 == 0 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)]
hk_trace = [(I-fI)[s] for s in range(frame) if s/int(period/dt)%2 == 1 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)]
lk_mean = []
hk_mean = []
dff = []
for i in range(int(t.max()/period/2)+1):
    lk_mean1 = np.mean(lk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    hk_mean1 = np.mean(hk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    dff1 = np.abs(lk_mean1-hk_mean1)/np.mean(fI[int(period/dt)*2*i:int(period/dt)*(2*i+2)])
    lk_mean = np.append(lk_mean, lk_mean1)
    hk_mean = np.append(hk_mean, hk_mean1)
    dff = np.append(dff, dff1)
print np.abs(lk_mean-hk_mean)
print dff

fig6, ax = plt.subplots()
ax.plot(t, (I-fI), 'c.', markersize=1, label='residual', alpha=0.5)
for i in range(int(t.max()/period/2)+1):
    l_tfrag = t[i*2*int(period/dt)+transiperiod:(i*2+1)*int(period/dt)-transiperiod]
    h_tfrag = t[(i*2+1)*int(period/dt)+transiperiod:(i*2+2)*int(period/dt)-transiperiod]
    ax.plot(l_tfrag, np.repeat(lk_mean[i],len(l_tfrag)), 'r-', label='low K mean {}'.format(i)+ '= {}'.format(round(lk_mean[i],4)))
    ax.plot(h_tfrag, np.repeat(hk_mean[i],len(h_tfrag)), 'b-', label='high K mean {}'.format(i)+ '= {}'.format(round(hk_mean[i],4)))
for i in range(int(t.max()/period)+1):
    ax.axvline(x = period*i, c='0.3', alpha=0.2)
ax.legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax.set_xlim(0, t.max())
ax.set_title('gaussian filter')
plt.xticks(np.arange(0, t.max(), period))
for i in range(int(t.max()/period/2)):
    ax.annotate(r'$\Delta$F/F={}'.format(round(dff[i],3)), xy=(0,0), xytext=(period*2*i,10), xycoords='axes pixels', fontsize=10)
fig5.canvas.draw()
#%%
if savefig == 1:
    fig1.savefig(filePath+abc+'_fig1_2expfit.png', format='png', bbox_inches = 'tight')
    fig2.savefig(filePath+abc+'_fig2_Fluor.png', format='png', bbox_inches = 'tight')
    fig3.savefig(filePath+abc+'_fig3_polyfit.png', format='png', bbox_inches = 'tight')
    fig4.savefig(filePath+abc+'_fig4_Fluor.png', format='png', bbox_inches = 'tight')
    fig5.savefig(filePath+abc+'_fig5_Gfilter.png', format='png', bbox_inches = 'tight')
    fig6.savefig(filePath+abc+'_fig6_Fluor.png', format='png', bbox_inches = 'tight')
#%%
# fit low K range and high K range separately to 2 exponential
def EXP2(x, A1, tau1, A2, tau2, b):
    return A1*np.exp(-tau1*x)+A2*np.exp(-tau2*x)+b
lk_I = np.array([(I)[s] for s in range(frame) if s/int(period/dt)%2 == 0 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])
lk_t = np.array([(t)[s] for s in range(frame) if s/int(period/dt)%2 == 0 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])
hk_I = np.array([(I)[s] for s in range(frame) if s/int(period/dt)%2 == 1 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])
hk_t = np.array([(t)[s] for s in range(frame) if s/int(period/dt)%2 == 1 and
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])

gmod = lmfit.Model(EXP2)
params = gmod.make_params()
params['A1'].set(value = lk_I.max()-lk_I.min(), min=0)
params['tau1'].set(value = -np.log(lk_I[len(lk_I)-1]/lk_I[0])/lk_t[len(lk_I)-1])
params['A2'].set(value = lk_I.max()-lk_I.min(), min=0)
params['tau2'].set(value = -np.log(lk_I[len(lk_I)-1]/lk_I[0])/lk_t[len(lk_I)-1])
params['b'].set(value = lk_I.min(), min=0)
lk_result = gmod.fit(lk_I, x=lk_t, **params)
lk_fI = EXP2(lk_t, lk_result.best_values['A1'], lk_result.best_values['tau1'],lk_result.best_values['A2'], lk_result.best_values['tau2'], lk_result.best_values['b'])
lk_fI_tot = EXP2(t, lk_result.best_values['A1'], lk_result.best_values['tau1'],lk_result.best_values['A2'], lk_result.best_values['tau2'], lk_result.best_values['b'])
params['A1'].set(value = hk_I.max()-hk_I.min(), min=0)
params['tau1'].set(value = -np.log(hk_I[len(hk_I)-1]/hk_I[0])/hk_t[len(hk_I)-1])
params['A2'].set(value = hk_I.max()-hk_I.min(), min=0)
params['tau2'].set(value = -np.log(hk_I[len(hk_I)-1]/hk_I[0])/hk_t[len(hk_I)-1])
params['b'].set(value = hk_I.min(), min=0)
hk_result = gmod.fit(hk_I, x=hk_t, **params)
hk_fI = EXP2(hk_t, hk_result.best_values['A1'], hk_result.best_values['tau1'],hk_result.best_values['A2'], hk_result.best_values['tau2'], hk_result.best_values['b'])
hk_fI_tot = EXP2(t, hk_result.best_values['A1'], hk_result.best_values['tau1'],hk_result.best_values['A2'], hk_result.best_values['tau2'], hk_result.best_values['b'])

fig7, ax = plt.subplots(2)
#result.plot_fit(ax=ax[0], datafmt='-', fitfmt='o-', initfmt='--', numpoints=40)
lk_result.plot_residuals(ax=ax[1], datafmt='c.')
hk_result.plot_residuals(ax=ax[1], datafmt='g.')
ax[0].plot(lk_t, lk_I, 'r.', markersize=2, label='lk_Raw Data', alpha=0.5)
ax[0].plot(t, lk_fI_tot, '-', label='lk_fit', color='k')
#ax[0].plot(lk_t, EXP2(lk_t, lk_I.max()-lk_I.min(), -np.log(lk_I[len(lk_I)-1]/lk_I[0])/lk_t[len(lk_I)-1], lk_I.max()-lk_I.min(), -np.log(lk_I[len(lk_I)-1]/lk_I[0])/lk_t[len(lk_I)-1], lk_I.min()), '--', label='lk_initial guess', color='b')
ax[0].plot(hk_t, hk_I, 'm.', markersize=2, label='hk_Raw Data', alpha=0.5)
ax[0].plot(t, hk_fI_tot, '-', label='hk_fit', color='0.7')
#ax[0].plot(hk_t, EXP2(hk_t, hk_I.max()-hk_I.min(), -np.log(hk_I[len(hk_I)-1]/hk_I[0])/hk_t[len(hk_I)-1], hk_I.max()-hk_I.min(), -np.log(hk_I[len(hk_I)-1]/hk_I[0])/hk_t[len(hk_I)-1], hk_I.min()), '--', label='hk_initial guess', color='c')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].set_xlim(0, t.max())
ax[1].set_xlim(0, t.max())

lk_trace = lk_I-lk_fI
hk_trace = hk_I-hk_fI
lk_mean = []
hk_mean = []
dff = []
for i in range(int(t.max()/period/2)+1):
    lk_mean1 = np.mean(lk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    hk_mean1 = np.mean(hk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    lk_f1 = np.mean(lk_fI[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    hk_f1 = np.mean(hk_fI[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    dff1 = np.abs(lk_mean1-hk_mean1)/(0.5*(lk_f1+hk_f1))
    lk_mean = np.append(lk_mean, lk_mean1)
    hk_mean = np.append(hk_mean, hk_mean1)
    dff = np.append(dff, dff1)
    print 0.5*(lk_f1+hk_f1)
print np.abs(lk_mean-hk_mean)
print dff

fig8, ax = plt.subplots()
ax.plot(lk_t, lk_I-lk_fI, 'c.', markersize=1, label='lk_residual', alpha=0.5)
ax.plot(hk_t, hk_I-hk_fI, 'g.', markersize=1, label='hk_residual', alpha=0.5)
for i in range(int(t.max()/period/2)+1):
    l_tfrag = lk_t[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)]
    h_tfrag = hk_t[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)]
    ax.plot(l_tfrag, np.repeat(lk_mean[i],len(l_tfrag)), 'r-', label='low K mean {}'.format(i)+ '= {}'.format(round(lk_mean[i],4)))
    ax.plot(h_tfrag, np.repeat(hk_mean[i],len(h_tfrag)), 'b-', label='high K mean {}'.format(i)+ '= {}'.format(round(hk_mean[i],4)))
for i in range(int(t.max()/period)+1):
    ax.axvline(x = period*i, c='0.3', alpha=0.2)
ax.legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax.set_xlim(0, t.max())
ax.set_title('fit 2 exponential separately')
plt.xticks(np.arange(0, t.max(), period))
for i in range(int(t.max()/period/2)):
    ax.annotate(r'$\Delta$F/F={}'.format(round(dff[i],3)), xy=(0,0), xytext=(period*2*i,10), xycoords='axes pixels', fontsize=10)
fig8.canvas.draw()
#%% calculate dF/F
lk_fI = EXP2(t, lk_result.best_values['A1'], lk_result.best_values['tau1'],lk_result.best_values['A2'], lk_result.best_values['tau2'], lk_result.best_values['b'])
hk_fI = EXP2(t, hk_result.best_values['A1'], hk_result.best_values['tau1'],hk_result.best_values['A2'], hk_result.best_values['tau2'], hk_result.best_values['b'])
dff = (lk_fI-hk_fI)/((lk_fI+hk_fI)/2)
fig9, ax = plt.subplots()
ax.plot(t,dff)
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'$\Delta$F/F')
ax.set_title(r'$\Delta$F/F v.s. time')
#%%
if savefig == 1:
    fig7.savefig(filePath+abc+'_fig7_2exp_sp.png', format='png', bbox_inches = 'tight')
    fig8.savefig(filePath+abc+'_fig8_Fluor.png', format='png', bbox_inches = 'tight')
    fig9.savefig(filePath+abc+'_fig9_dFF.png', format='png', bbox_inches = 'tight')
#%%
# apply moving average
def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
window_size=10

lk_I = np.array([(I)[s] for s in range(frame) if s/int(period/dt)%2 == 0])
lk_t = np.array([(t)[s] for s in range(frame) if s/int(period/dt)%2 == 0])
hk_I = np.array([(I)[s] for s in range(frame) if s/int(period/dt)%2 == 1])
hk_t = np.array([(t)[s] for s in range(frame) if s/int(period/dt)%2 == 1])
lk_Is = []
hk_Is = []
for i in range(int(t.max()/period/2)+1):
    lk_Is1 = movingaverage(lk_I[i*int(period/dt):(i+1)*int(period/dt)], window_size)
    lk_Is = np.append(lk_Is, lk_Is1)
for i in range(int(t.max()/period/2)):
    hk_Is1 = movingaverage(hk_I[i*int(period/dt):(i+1)*int(period/dt)], window_size)
    hk_Is = np.append(hk_Is, hk_Is1)

fig10, ax = plt.subplots()
ax.plot(lk_t, lk_I, 'r.', markersize=2, label='lk_Raw Data', alpha=0.5)
ax.plot(lk_t, lk_Is, 'k.', markersize=2, label='lk_smooth')
ax.plot(hk_t, hk_I, 'm.', markersize=2, label='hk_Raw Data', alpha=0.5)
ax.plot(hk_t, hk_Is, '.', markersize=2, label='hk_smooth', color='0.7')
ax.legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Fluorescence intensity')
ax.set_xlim(0, t.max())
#%%
# discard data points near the boundary and the buffer transiperiod
def EXP2(x, A1, tau1, A2, tau2, b):
    return A1*np.exp(-tau1*x)+A2*np.exp(-tau2*x)+b
lk_Is = np.array([(lk_Is)[s] for s in range(len(lk_Is)) if
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])
lk_I = np.array([(lk_I)[s] for s in range(len(lk_I)) if
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])
lk_t = np.array([(lk_t)[s] for s in range(len(lk_t)) if
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])
hk_Is = np.array([(hk_Is)[s] for s in range(len(hk_Is)) if
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])
hk_I = np.array([(hk_I)[s] for s in range(len(hk_I)) if
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])
hk_t = np.array([(hk_t)[s] for s in range(len(hk_t)) if
            s%(period/dt) > transiperiod and s%(period/dt) < ((period/dt)-transiperiod)])

if len(lk_I) > len(hk_I):
    lk_I = lk_I[:(len(lk_I)-window_size)]
    lk_Is = lk_Is[:(len(lk_Is)-window_size)]
    lk_t = lk_t[:(len(lk_t)-window_size)]
if len(lk_I) < len(hk_I):
    hk_I = hk_I[:(len(hk_I)-window_size)]
    hk_Is = hk_Is[:(len(hk_Is)-window_size)]
    hk_t = hk_t[:(len(hk_t)-window_size)]
#%%
# fit low K range to 2 exponential and calculate dF/F
gmod = lmfit.Model(EXP2)
params = gmod.make_params()
params['A1'].set(value = lk_Is.max()-lk_Is.min(), min=0)
params['tau1'].set(value = -np.log(lk_Is[len(lk_Is)-1]/lk_Is[0])/lk_t[len(lk_Is)-1])
params['A2'].set(value = lk_Is.max()-lk_Is.min(), min=0)
params['tau2'].set(value = -np.log(lk_Is[len(lk_Is)-1]/lk_Is[0])/lk_t[len(lk_Is)-1])
params['b'].set(value = lk_Is.min(), min=0)
lk_result = gmod.fit(lk_Is, x=lk_t, **params)
lk_fIs = EXP2(lk_t, lk_result.best_values['A1'], lk_result.best_values['tau1'],lk_result.best_values['A2'], lk_result.best_values['tau2'], lk_result.best_values['b'])
hk_fIs = EXP2(hk_t, lk_result.best_values['A1'], lk_result.best_values['tau1'],lk_result.best_values['A2'], lk_result.best_values['tau2'], lk_result.best_values['b'])
tot_fIs = EXP2(t, lk_result.best_values['A1'], lk_result.best_values['tau1'],lk_result.best_values['A2'], lk_result.best_values['tau2'], lk_result.best_values['b'])

lk_trace = lk_Is-lk_fIs
hk_trace = hk_Is-hk_fIs
lk_mean = []
hk_mean = []
for i in range(int(t.max()/period/2)+1):
    lk_mean1 = np.mean(lk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    hk_mean1 = np.mean(hk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    lk_f1 = np.mean(lk_fIs[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    hk_f1 = np.mean(hk_fIs[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
    lk_mean = np.append(lk_mean, lk_mean1)
    hk_mean = np.append(hk_mean, hk_mean1)

fig11, ax = plt.subplots(2, figsize = (10,10))
ax[0].plot(lk_t, lk_I, 'r.', markersize=2, label='lk_raw Data', alpha=0.5)
ax[0].plot(lk_t, lk_Is, 'k.', markersize=2, label='lk_smooth', alpha=0.5)
ax[0].plot(hk_t, hk_I, 'm.', markersize=2, label='hk_raw Data', alpha=0.5)
ax[0].plot(hk_t, hk_Is, '.', markersize=2, label='hk_smooth', color='0.7')
ax[0].plot(t, tot_fIs, 'b-', label='lk_smooth fit')
ax[1].plot(lk_t, lk_I  - lk_fIs, 'r.', markersize=2, label='lk_residual', alpha=0.5)
ax[1].plot(lk_t, lk_Is - lk_fIs, 'k.', markersize=2, label='lk_smooth residual', alpha=0.5)
ax[1].plot(hk_t, hk_I  - hk_fIs, 'm.', markersize=2, label='hk_residual', alpha=0.5)
ax[1].plot(hk_t, hk_Is - hk_fIs, '.', markersize=2, label='hk_smooth residual', color='0.7')
ax[1].axhline(y=0, c='b')
ax[0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
ax[0].set_title('fit low K part to 2 exponential')
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Fluorescence intensity')
ax[0].set_xlim(0, t.max())
ax[1].set_xlabel('Time (s)')
ax[1].set_ylabel('Residuals')
ax[1].set_xlim(0, t.max())

for i in range(int(t.max()/period/2)+1):
    l_tfrag = lk_t[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)]
    h_tfrag = hk_t[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)]
    ax[1].plot(l_tfrag, np.repeat(lk_mean[i],len(l_tfrag)), 'c-', linewidth=4, label='low K mean {}'.format(i)+ '= {}'.format(round(lk_mean[i],4)))
    ax[1].plot(h_tfrag, np.repeat(hk_mean[i],len(h_tfrag)), 'g-', linewidth=4, label='high K mean {}'.format(i)+ '= {}'.format(round(hk_mean[i],4)))
for i in range(int(t.max()/period)+1):
    ax[0].axvline(x = period*i, c='0.3', alpha=0.2)
    ax[1].axvline(x = period*i, c='0.3', alpha=0.2)
#ax[1].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
plt.xticks(np.arange(0, t.max(), period))
#for i in range(int(t.max()/period/2)):
#    ax[1].annotate(r'$\Delta$F/F={}'.format(round(dff[i],3)), xy=(0,0), xytext=(period*2*i,10), xycoords='axes pixels', fontsize=10)
fig11.tight_layout()
fig11.canvas.draw()