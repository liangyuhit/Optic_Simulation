# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

screen_diameter = 2e-3
I_0 = 1
angle_beta = np.arcsin(4*Lamda/screen_diameter)
angle_alpha = angle_beta / 2
print('Angle_Alpha = %e [rad]'%angle_alpha)

L = 0.2
D = np.linspace(0, 1, num=301) * Lamda
d = np.linspace(screen_diameter/2*(-1), screen_diameter/2, num=301)


X, Y = np.meshgrid(d, D)

D_d = Y - X * np.tan(angle_alpha)
L_d = L + Y * np.tan(angle_alpha)

diff_L_d = 2 * D_d + 2 * L_d * np.tan(angle_alpha)**2 / (1 + np.tan(angle_alpha)**2)
Z = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_d / Lamda))

if 1:
    plt.figure(1)
    c = plt.gca().pcolor(X*1000, Y/Lamda, Z, cmap='RdBu')
    plt.xlabel('d [mm]')
    plt.ylabel('D/Lamda')
    plt.colorbar(c, ax=plt.gca())
#     fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.title('Intensity Distribution with moving mirror (One tilted mirror)')
    plt.show()


