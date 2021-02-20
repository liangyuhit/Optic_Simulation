
# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np


c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]



angle_beta_2 = np.linspace(0, np.pi/200, 10001)
angle_alpha_2 = angle_beta_2 / 2
# print('Angle_Beta_2 = %e [rad]'%angle_beta_2)
k = np.cos(angle_alpha_2)**2
# print(k)


if 1: # Center Point
    plt.figure(1)
    plt.title('(Non)Linearty')
    plt.plot(angle_beta_2, k, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('beta_2 [rad]')
#     plt.ylabel('Linearity')
#     [i.set_color("blue") for i in plt.gca().get_yticklabels()]
#     plt.ylim(0.99995,1)
    plt.grid(which = 'both')
    
    plt.show()
