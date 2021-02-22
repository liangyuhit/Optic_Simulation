# -*- coding: utf-8 -*-
'''
Created on Feb 22, 2021

@author: yl
'''
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
from plt_style import *
colors = color_set_medium

'''
    readout
'''
npzfile = np.load('2D_shift_tilt.npy')
print(npzfile.files)
M_m_set = npzfile['M_m_set']
V_m_x_set = npzfile['V_m_x_set']
w_0_set = npzfile['w_0_set']
contrast_set = npzfile['contrast_set']
print(np.shape(contrast_set))

'''
    plot
''' 
fig = plt.figure('平行度 & 平移')
plt.gcf().set_size_inches(15/2.54, 14/2.54)
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax2 = fig.add_subplot(2, 2, 2, projection='3d',sharex=ax1,sharey=ax1,sharez=ax1)
ax3 = fig.add_subplot(2, 2, 3, projection='3d',sharex=ax1,sharey=ax1,sharez=ax1)
ax4 = fig.add_subplot(2, 2, 4, projection='3d',sharex=ax1,sharey=ax1,sharez=ax3)

V_m_x_set, M_m_set = np.meshgrid(V_m_x_set, M_m_set)
print(np.shape(M_m_set))
print(np.shape(V_m_x_set))
print(np.shape(contrast_set[0]))
print(1)
axes = [ax1,ax2,ax3,ax4]
for i in range(3):
    ax = axes[i]
    ax.plot_surface(M_m_set*1e3, V_m_x_set*1e6, contrast_set[i], cmap='jet', rstride=1, cstride=1,
                       linewidth=0, antialiased=False)
    ax.set_ylabel('角度偏差$\mathrm{(\u03BCrad)}$')
    ax.set_xlabel('平移量$\mathrm{(mm)}$')
#     ax.set_zlabel('对比度')
ax1.view_init(20, -50)
ax2.view_init(20, -50)
ax3.view_init(20, -50)
ax4.view_init(20, -50)
plt.show()
# plt.savefig(r'/Users/yl/Documents/Eclipse Workspace/Optic_Simulation/Heterodyne_Contrast/2D_tilt_shift.jpg', dpi=600) #指定分辨率保存
# plt.savefig(r'C:\Users\yuxiaoyang\Desktop\Fig_4_10.jpg', dpi=600) #指定分辨率保存
