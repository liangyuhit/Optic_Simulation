# -*- coding: utf-8 -*-
'''
Created on Feb 22, 2021

@author: yl
'''
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from plt_style import *
colors = color_set_medium

'''
    readout
'''
Lamda = 633e-9
npzfile = np.load('1D_shift.npy')
print(npzfile.files)
M_m_set = npzfile['M_m_set']
w_0_set = npzfile['w_0_set']
contrast_set = npzfile['contrast_set']

'''
    plot
''' 
fig = plt.figure('光斑平移')
# plt.gcf().set_size_inches(15/2.54, 7/2.54)
colorset = [colors['black'],colors['blue'],colors['red'],colors['orange'],colors['green']]
markerset = ['o','s','D','^','v']

for i in range(len(w_0_set)):
    plt.plot(M_m_set*1e3, contrast_set[i], color=colorset[i],label='光斑直径%i$\mathrm{(mm)}$'%(w_0_set[i]*1e3*2))    
plt.ylabel('对比度')
plt.xlabel('平移量$\mathrm{(mm)}$')
plt.legend(loc='upper left')

plt.show()
# plt.savefig(r'/Users/yl/Documents/Eclipse Workspace/Optic_Simulation/Heterodyne_Contrast/1D_shift.jpg', dpi=600) #指定分辨率保存
# plt.savefig(r'C:\Users\yuxiaoyang\Desktop\Fig_4_10.jpg', dpi=600) #指定分辨率保存
