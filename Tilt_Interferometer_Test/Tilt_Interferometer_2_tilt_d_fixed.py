# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np



j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

screen_diameter = 2e-3
I_0 = 1
angle_beta = np.arcsin(10*Lamda/screen_diameter)
angle_alpha = angle_beta / 2
print('Angle_alpha = %e [rad]'%angle_alpha)

angle_beta_2 = angle_beta * (0.5)
angle_alpha_2 = angle_beta_2 / 2
print('Angle_Alpha_2 = %e [rad]'%angle_alpha_2)


L = 0.2
D = 0.1
d = np.linspace(screen_diameter/2*(-1), screen_diameter/2, num=301)

D_d = D - d * np.tan(angle_alpha) + d * np.tan(angle_alpha_2)
L_d = L + d * np.tan(angle_alpha)

# diff_L_d_2 = D + d * (np.tan(angle_alpha_2) - np.tan(angle_alpha)) + (L + D + d*np.tan(angle_alpha_2)) * (1 - np.tan(angle_alpha_2)**2) / (1 + np.tan(angle_alpha_2)**2) - (L + d*np.tan(angle_alpha)) * (1 - np.tan(angle_alpha)**2) / (1 + np.tan(angle_alpha)**2)
# diff_L_d_2 = (D + d*np.tan(angle_alpha_2)) * (1 + np.cos(2*angle_alpha_2)) + L * (np.cos(2*angle_alpha_2) - np.cos(2*angle_alpha)) - d * np.tan(angle_alpha) * (1 + np.cos(2*angle_alpha))
diff_L_d_2 = D_d + (L_d + D_d) * np.cos(2 * angle_alpha_2) - L_d * np.cos(2 * angle_alpha)
I_beat_d_2 = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda))



if 1: # Center Point
    plt.figure(1)
    plt.title('Intensity Distribution(Two tilted mirror) beta2 = beta')
    plt.plot(d*1000, I_beat_d_2, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Screen [mm]')
    plt.ylabel('Intensity [A.U.]')
#     [i.set_color("blue") for i in plt.gca().get_yticklabels()]
    plt.ylim(-0.2,4.2)
    plt.grid(which = 'both')    
    
    plt.show()
