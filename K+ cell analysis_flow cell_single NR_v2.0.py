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
#from mpl_toolkits.mplot3d import Axes3D

filePath = '/Users/yungkuo/Documents/Data/K+ cell/flow chamber/101815 flow cell/ZNR/'
savefig = 1
dt = 1/32.352
scan = 5
bgscan = 30
assignpts = 1
period = 60
transiperiod = 200  # discard (transiperiod) frames
                    # befroe and after the transition of low/high K buffer
                    # e.g. discard 1941-(transiperiod) to 1941+(transiperiod)
                    #      timetrace points when calculating the mean
def listdir_nohidden(path):
    for f in os.listdir(path):
        if f.startswith('ZNR_low K_60s_high K_60s_X'):
            yield f

# assign QD coordinates
if assignpts == 0:
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

if assignpts == 1:
    ptss = np.array([[372, 297],
                    [244, 218],
                    [305, 330],
                    [394, 370],
                    [287, 453],
                    [175, 126],
                    [254, 332],
                    [384, 402],
                    [219, 400],
                    [116, 318]])

# Extract QD signals and background
for i in range(len(ptss)):
    abc = 'QD{}'.format(i)
    pts = ptss[[i]]
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
            #if i%2000 == 0:
            #    fig0, ax = plt.subplots()
            #    ax.imshow(img[pts[:,1]-scan:pts[:,1]+scan,pts[:,0]-scan:pts[:,0]+scan], interpolation='none', cmap='gray')
        frame = frame + data[0]
    t = np.arange(0,frame*dt, dt)

    # plot fluorescence and background v.s. time
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
    I = I-bg
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

    # plot data and smoothed data
    fig10, ax = plt.subplots()
    ax.plot(lk_t, lk_I, 'r.', markersize=2, label='lk_Raw Data', alpha=0.5)
    ax.plot(lk_t, lk_Is, 'k.', markersize=2, label='lk_smooth')
    ax.plot(hk_t, hk_I, 'm.', markersize=2, label='hk_Raw Data', alpha=0.5)
    ax.plot(hk_t, hk_Is, '.', markersize=2, label='hk_smooth', color='0.7')
    ax.legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Fluorescence intensity')
    ax.set_xlim(0, t.max())

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

    if len(lk_I)/int(period/dt) > len(hk_I)/int(period/dt):
        lk_I = lk_I[:(len(lk_I)-window_size)]
        lk_Is = lk_Is[:(len(lk_Is)-window_size)]
        lk_t = lk_t[:(len(lk_t)-window_size)]
    if len(lk_I)/int(period/dt) < len(hk_I)/int(period/dt):
        hk_I = hk_I[:(len(hk_I)-window_size)]
        hk_Is = hk_Is[:(len(hk_Is)-window_size)]
        hk_t = hk_t[:(len(hk_t)-window_size)]
    if len(lk_I)/int(period/dt) == len(hk_I)/int(period/dt):
        lk_I = lk_I[:(len(lk_I)-window_size)]
        lk_Is = lk_Is[:(len(lk_Is)-window_size)]
        lk_t = lk_t[:(len(lk_t)-window_size)]
        hk_I = hk_I[:(len(hk_I)-window_size)]
        hk_Is = hk_Is[:(len(hk_Is)-window_size)]
        hk_t = hk_t[:(len(hk_t)-window_size)]

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

    lk_trace = lk_I-lk_fIs
    hk_trace = hk_I-hk_fIs
    lks_trace = lk_Is-lk_fIs
    hks_trace = hk_Is-hk_fIs
    lk_mean = []
    hk_mean = []
    lks_mean = []
    hks_mean = []
    for i in range(int(t.max()/period/2)+1):
        lk_mean1 = np.mean(lk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
        hk_mean1 = np.mean(hk_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
        lks_mean1 = np.mean(lks_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
        hks_mean1 = np.mean(hks_trace[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
        lk_mean = np.append(lk_mean, lk_mean1)
        hk_mean = np.append(hk_mean, hk_mean1)
        lks_mean = np.append(lks_mean, lks_mean1)
        hks_mean = np.append(hks_mean, hks_mean1)
    # Calcualte dF/F
    # EXTRA ATTENTION !! CASE SPECIFIC !! Depending on which period chosen to compare.
    # in this case, we skipped first low K period and compared highK1 to lowK2, highK2 to lowK3...
    dF = np.array(hks_mean[~np.isnan(hks_mean)]-lks_mean[1:])
    F = []
    for i in range(int(t.max()/period/2)+1):
        F1 = np.mean(hk_fIs[(i)*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)])
        F = np.append(F, F1)
    dFF = dF/F[~np.isnan(F)]
    x = np.arange(period*1.5, t.max(), period*2)
    # plot results
    fig11, ax = plt.subplots(2,5, figsize = (20,10))
    ax[0,0] = plt.subplot2grid((2,5), (0,0), colspan=3)
    ax[1,0] = plt.subplot2grid((2,5), (1,0), colspan=3)
    ax[0,1] = plt.subplot2grid((2,5), (0,3), colspan=2)
    ax[1,1] = plt.subplot2grid((2,5), (1,3))
    ax[1,2] = plt.subplot2grid((2,5), (1,4))
    # plot fit and residuals
    ax[0,0].plot(lk_t, lk_I, 'r.', markersize=2, label='lk_raw Data', alpha=0.5)
    ax[0,0].plot(lk_t, lk_Is, 'k.', markersize=2, label='lk_smooth', alpha=0.5)
    ax[0,0].plot(hk_t, hk_I, 'm.', markersize=2, label='hk_raw Data', alpha=0.5)
    ax[0,0].plot(hk_t, hk_Is, '.', markersize=2, label='hk_smooth', color='0.7')
    ax[0,0].plot(t, tot_fIs, 'b-', label='lk_smooth fit')
    ax[1,0].plot(lk_t, lk_I  - lk_fIs, 'r.', markersize=2, label='lk_residual', alpha=0.5)
    ax[1,0].plot(lk_t, lk_Is - lk_fIs, 'k.', markersize=2, label='lk_smooth residual', alpha=0.5)
    ax[1,0].plot(hk_t, hk_I  - hk_fIs, 'm.', markersize=2, label='hk_residual', alpha=0.5)
    ax[1,0].plot(hk_t, hk_Is - hk_fIs, '.', markersize=2, label='hk_smooth residual', color='0.7')
    ax[1,0].axhline(y=0, c='b')
    ax[0,0].legend(bbox_to_anchor=(1, 1), frameon=False, fontsize=10)
    ax[0,0].set_title('fit low K part to 2 exponential')
    ax[0,0].set_xlabel('Time (s)')
    ax[0,0].set_ylabel('Fluorescence intensity')
    ax[0,0].set_xlim(0, t.max())
    ax[1,0].set_xlabel('Time (s)')
    ax[1,0].set_ylabel('Residuals')
    ax[1,0].set_xlim(0, t.max())
    for i in range(int(t.max()/period/2)+1):
        l_tfrag = lk_t[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)]
        h_tfrag = hk_t[i*int(period/dt-transiperiod*2):(i+1)*int(period/dt-transiperiod*2)]
        ax[1,0].plot(l_tfrag, np.repeat(lk_mean[i],len(l_tfrag)), 'c-', linewidth=2, label='low K mean {}'.format(i)+ '= {}'.format(round(lk_mean[i],4)))
        ax[1,0].plot(h_tfrag, np.repeat(hk_mean[i],len(h_tfrag)), 'g-', linewidth=2, label='high K mean {}'.format(i)+ '= {}'.format(round(hk_mean[i],4)))
    # plot dF/F
    ax[0,1].plot(np.arange(period*0.5, t.max(), period*2), lk_mean, 'ro-')
    ax[0,1].plot(np.arange(period*0.5, t.max(), period*2), lks_mean, 'c^-')
    ax[0,1].plot(np.arange(period*1.5, t.max(), period*2), hk_mean[~np.isnan(hk_mean)], 'mo-')
    ax[0,1].plot(np.arange(period*1.5, t.max(), period*2), hks_mean[~np.isnan(hks_mean)], 'g^-')
    ax[0,1].set_xlabel('Time(s)')
    ax[0,1].set_ylabel('Residuals mean')
    ax01_2 = ax[0,1].twinx()
    ax01_2.plot(x, dFF, 'bs-')
    ax01_2.set_ylabel(r'$\Delta$F/F', color='b')
    for tl in ax01_2.get_yticklabels():
        tl.set_color('b')
    # plot residual distribution histogram
    n, bins, patches = ax[1,1].hist((lk_I  - lk_fIs), bins=100, histtype='stepfilled', color = 'r', alpha=0.5, orientation='horizontal')
    n, bins, patches = ax[1,1].hist((hk_I  - hk_fIs), bins=bins, histtype='stepfilled', color = 'm', alpha=0.5, orientation='horizontal')
    n, bins, patches = ax[1,2].hist((lk_Is - lk_fIs), bins=100, histtype='stepfilled', color = 'k', alpha=0.5, orientation='horizontal')
    n, bins, patches = ax[1,2].hist((hk_Is - hk_fIs), bins=bins, histtype='stepfilled', color = '0.7', alpha=0.5, orientation='horizontal')
    ax[1,1].axhline(y=np.mean(lk_trace), color='r')
    ax[1,1].axhline(y=np.mean(hk_trace), color='m')
    ax[1,2].axhline(y=np.mean(lks_trace), color='c')
    ax[1,2].axhline(y=np.mean(hks_trace), color='g')
    ax[1,1].set_xlabel('Counts')
    ax[1,2].set_xlabel('Counts')
    # grids and ticks matching the period
    ax[0,0].set_xticks(np.arange(0, t.max(), period))
    ax[1,0].set_xticks(np.arange(0, t.max(), period))
    ax[0,1].set_xticks(np.arange(0, t.max()+period, period))
    ax[0,0].grid(color='0.1', linestyle='-', axis='x', linewidth=0.5, alpha=0.2)
    ax[1,0].grid(color='0.1', linestyle='-', axis='x', linewidth=0.5, alpha=0.2)
    ax[0,1].grid(color='0.1', linestyle='-', axis='x', linewidth=0.5, alpha=0.2)
    #for i in range(int(t.max()/period/2)):
    #    ax[1].annotate(r'$\Delta$F/F={}'.format(round(dff[i],3)), xy=(0,0), xytext=(period*2*i,10), xycoords='axes pixels', fontsize=10)
    fig11.tight_layout()
    fig11.canvas.draw()

    # save figures
    if savefig == 1:
        fig1.savefig(filePath+'fig1 and 2/'+abc+'_fig1_img.png', format='png', bbox_inches = 'tight')
        fig2.savefig(filePath+'fig1 and 2/'+abc+'_fig2_fluor.png', format='png', bbox_inches = 'tight')
        fig10.savefig(filePath+'fig10/'+abc+'_fig10_movingavg.png', format='png', bbox_inches = 'tight')
        fig11.savefig(filePath+'fig11/'+abc+'_fig11_fit and residuals.png', format='png', bbox_inches = 'tight')
