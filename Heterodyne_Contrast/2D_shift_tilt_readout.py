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

'''
    plot
''' 
fig1 = plt.figure('平行度 & 平移（3D）')
plt.gcf().set_size_inches(25/2.54, 8/2.54)
ax1 = fig1.add_subplot(1, 3, 1, projection='3d')
ax2 = fig1.add_subplot(1, 3, 2, projection='3d',sharex=ax1,sharey=ax1,sharez=ax1)
ax3 = fig1.add_subplot(1, 3, 3, projection='3d',sharex=ax1,sharey=ax1,sharez=ax1)
 
V_m_x_set, M_m_set = np.meshgrid(V_m_x_set, M_m_set)
axes = [ax1,ax2,ax3]
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
# plt.tight_layout()
# plt.show()
plt.savefig(r'/Users/yl/Documents/Eclipse Workspace/Optic_Simulation/Heterodyne_Contrast/2D_tilt_shift_surf.jpg', dpi=600) #指定分辨率保存
# plt.savefig(r'C:\Users\yuxiaoyang\Desktop\Fig_4_10.jpg', dpi=600) #指定分辨率保存



fig2 = plt.figure('平行度 & 平移（2D）')
plt.gcf().set_size_inches(25/2.54, 8/2.54)
ax4 = fig2.add_subplot(1, 3, 1)
ax5 = fig2.add_subplot(1, 3, 2, sharex=ax4,sharey=ax4)
ax6 = fig2.add_subplot(1, 3, 3, sharex=ax4,sharey=ax4)
V_m_x_set, M_m_set = np.meshgrid(V_m_x_set, M_m_set)
axes = [ax4,ax5,ax6]
for i in range(3):
    ax = axes[i]
    extent = []
    ax.imshow(contrast_set[i],origin='lower', cmap='jet')
    ax.set_ylabel('角度偏差$\mathrm{(\u03BCrad)}$')
    ax.set_xlabel('平移量$\mathrm{(mm)}$')
plt.tight_layout()
# plt.show()   
plt.savefig(r'/Users/yl/Documents/Eclipse Workspace/Optic_Simulation/Heterodyne_Contrast/2D_tilt_shift_imshow.jpg', dpi=600) #指定分辨率保存

    
    
    
    
