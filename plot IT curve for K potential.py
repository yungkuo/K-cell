# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 17:07:49 2015

@author: yungkuo
"""

import numpy as np
import matplotlib.pyplot as plt
filePath1 = '/Users/yungkuo/Documents/Data/K+ cell/091515 TMRM DiO HEK K potential/TMRM low K_500nM_30min_DiO_12uM_VA_0.2uM_add K__1/Pos0/'
filePath2 = '/Users/yungkuo/Documents/Data/K+ cell/091515 TMRM DiO HEK K potential/TMRM low K_500nM_30min_DiO_12uM_VA_0.2uM_add K__2/Pos0/'
filePath3 = '/Users/yungkuo/Documents/Data/K+ cell/091515 TMRM DiO HEK K potential/TMRM low K_500nM_30min_DiO_12uM_VA_0.2uM_add K__3/Pos0/'
I1 = np.load(filePath1+'IT_result.npy')
I2 = np.load(filePath2+'IT_result.npy')
I3 = np.load(filePath3+'IT_result.npy')
I = np.append(np.append(I1,I2),I3)
t = np.arange(0,len(I),1)*30
fig, ax = plt.subplots()
ax.plot(t, I, '-o')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Fluorescence Intensity (a.u.)')
