
# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np


c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]



angle_beta = np.linspace(0.001, 0.01, 10001)
# angle_beta = 3.165e-3
angle_alpha = angle_beta / 2
# print('Angle_Beta = %e [rad]'%angle_beta)
d_repeat = Lamda/np.sin(2*angle_alpha)


if 1: # Center Point
    plt.figure(1)
    plt.title('fringe distance')
    plt.plot(angle_beta, d_repeat*1000, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)

    plt.xlabel('beta [rad]')
    plt.ylabel('d_repeat [mm]')
#     [i.set_color("blue") for i in plt.gca().get_yticklabels()]
#     plt.ylim(-0.2,4.2)
    plt.grid(which = 'both')
    
    plt.show()
