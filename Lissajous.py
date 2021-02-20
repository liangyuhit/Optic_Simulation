# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np
from Ellipse_Fit_Func import Ellipse_fit


timeline = np.linspace(0, 1, 1001)

fig_num = 9
fig_color = ['red', 'blue', 'green', 'magenta', 'yellow', 'cyan', 'grey', 'orange', 'brown']
fig_set = np.linspace(0, fig_num-1, fig_num)

Amp_1, Phase_1, Freq_1 = 1, 1e-2, 1
Amp_2, Phase_2, Freq_2 = (Amp_1 + 0*fig_set/(fig_num-1) * Amp_1/10), (0+ 1*fig_set/(fig_num-1) * np.pi*2/2), (Freq_1 + 0*fig_set/(fig_num-1) * Freq_1/10)
I_0 = 1

I_1 = I_0 + Amp_1 * np.cos(2 * np.pi * Freq_1 * timeline + Phase_1)

I_2 = []
for k in range(fig_num):
    I = I_0 + Amp_2[k] * np.cos(2 * np.pi * Freq_2[k] * timeline + Phase_2[k])
    I_2.append(I)
    
     
fit_xaxis = fit_yaxis = np.linspace(0, 2, 1001)
fit_xaxis, fit_yaxis = np.meshgrid(fit_xaxis, fit_yaxis)
fit_ellipse = []
x0_set = []
y0_set = []
a_set = []
b_set = []
Delta_set = []



text = []

    
# for i in range(9):
#     fit = Ellipse_fit(I_1, I_2[i])
#     parameters = np.array(fit.fitEllipse())
#     print(parameters)
#     A, B, C, D, F, G = parameters
#     print(parameters/A)
#     fit_result = A*(fit_xaxis**2) + B*fit_xaxis*fit_yaxis + C*(fit_yaxis**2) +D*fit_xaxis + F*fit_yaxis + G
#     fit_ellipse.append(fit_result)



for i in range(9):
    fit = Ellipse_fit(I_1, I_2[i])
    parameters = np.array(fit.fitEllipse())
    A, B, C, D, F, G = parameters
    A, B, C, D, E, F = parameters/A
    fit_result = (fit_xaxis**2) + B*fit_xaxis*fit_yaxis + C*(fit_yaxis**2) +D*fit_xaxis + E*fit_yaxis + F
    x0 = (2*C*D - B*E)/(B**2-4*C)
    y0 = (2*E - B*D)/(B**2-4*C)
    a = ((x0**2 + (y0**2)*C + x0*y0*B -F)/(1-(B**2)/(4*C)))**0.5
    b = ((a**2)/C)**0.5
    Delta = np.arccos(-B/2/(C**0.5))
    
    fit_ellipse.append(fit_result)
    x0_set.append(x0)
    y0_set.append(y0)
    a_set.append(a)
    b_set.append(b)
    Delta_set.append(Delta)
    
    print('(%f,'%x0 + '%f);'%y0 + '(%f,'%a + '%f)'%b + 'Delta=%f'%Delta)
    
plt.figure('Intensity')
plt.plot(timeline, I_1, color='black', linestyle='solid', linewidth=5,  marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
for k in range(fig_num):
    plt.plot(timeline, I_2[k], color=fig_color[k], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Time')
    plt.xlabel('Intensity')
plt.get_current_fig_manager().window.setGeometry(20,100,640, 545)


plt.figure('Lissajous')
for k in range(fig_num):
    plt.subplot(3,3,k+1)
    plt.plot(I_1, I_2[k], color=fig_color[k], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('black')
#     plt.title('Phase Delay = %f [deg]'%(k * np.pi/4 /(2*np.pi) * 360))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.5)
plt.get_current_fig_manager().window.setGeometry(680,30, 1000, 1000)

plt.figure('Linear Fit')
for k in range(fig_num):
    plt.subplot(3,3,k+1)
    plt.contour(fit_xaxis, fit_yaxis, fit_ellipse[k], [0], colors=fig_color[k])     #x**2 + y**2 = 9 ��Բ��
    plt.xlabel('black')
#     plt.title('Phase Delay = %f [deg]'%(k * np.pi/4 /(2*np.pi) * 360))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.5)
plt.get_current_fig_manager().window.setGeometry(680,30, 1000, 1000)


plt.show()