# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np
from Ellipse_Fit_Func import Ellipse_fit

V_y_2_set = []

a_set = []
b_set = []
for i in range(7):
    a_set.append([])
    b_set.append([])

for m in range(16):
    j = complex(0, 1)
    c = 3e8 # 光速 [m/s]
    timeline = np.linspace(0,1,501)
    
    Lamda = 633e-9 # 光波长 [m]
    Fc = c / Lamda # 光频率 [Hz]
    
    I_0 = 1
    
    V_x, V_y, Vz = 0.0002*633/600*1.8, 0+1e-8, 1
    V_x_2, V_y_2, Vz_2 = 0.00, 0.00,  1
    
    
    
    L = 0.2 
    D = 0 + timeline * 0.5 * Lamda
    R = 0.5
    T = 1-R
    
    
    fig_color = ['black', 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'brown', 'grey', 'yellow']
    
    
    dx = [-0.25e-3, 0, 0.25e-3, 0.5e-3]
    # dy = [0]
    dy = [0, -0.5e-3]
    I_origin = []


    if 1:    ###modulation
        V_y_2 = m * 1e-6
        V_y_2_set.append(V_y_2)
    
    print('V_y_2 = %f'%V_y_2)
    infomation = 'Ref. (%f,'%V_x + ' %f,'%V_y + '%f)\n'%Vz + 'Mea. (%f,'%V_x_2 + ' %f,'%V_y_2 + '%f)\n'%Vz_2

    
    ### Origin Values
    for i in range(len(dy)):
        for j in range(len(dx)):
            D_d = D - dx[j] * np.tan(V_x_2 - V_x) + dy[i] * np.tan(V_y_2 - V_y)
            L_d = L + dx[j] * V_x + dy[i] * V_y
            diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
            a = np.array(2 * I_0 *2*(R*T)**0.5 * (1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda)))
            I_origin.append(a)
    
    
    ##### Ecllipse Fitting
    
    fit_xaxis = fit_yaxis = np.linspace(0, 4, 501)
    fit_xaxis, fit_yaxis = np.meshgrid(fit_xaxis, fit_yaxis)
    fit_ellipse = []
    text = []
    
    
    for i in range(len(dx)*len(dy)-1):
        fit = Ellipse_fit(I_origin[0], I_origin[i+1])
        parameters = np.array(fit.ellipse_get_parameters())
        x0 = parameters[0]
        y0 = parameters[1]
        a = parameters[2]
        b = parameters[3]
        phi = parameters[4]
    #     print(text)
    #     print('For ellipse_%i:'%(i+1) + 'Center(%f,'%x0 + '%f),'%y0 + 'a=%f,'%a + 'b=%f,'%b + 'phi=%f[deg]'%(phi/np.pi*180))
        fit_result = ((fit_xaxis - x0) * np.cos(phi) + (fit_yaxis - y0) * np.sin(phi))**2 / a**2 + ((fit_yaxis - y0) * np.cos(phi) - (fit_xaxis - x0) * np.sin(phi))**2 / b**2
        fit_ellipse.append(fit_result)
        if a < b:
            phi = phi * (-1)
            a, b = parameters[3], parameters[2]
        
        a_set[i].append(a)
        b_set[i].append(b)
        text.append(('Center(%f,'%x0 + '%f)\n'%y0 + 'Axes(%f,'%a + '%f)\n'%b + '%f[deg]'%(phi/np.pi*180)))
        
        
    if 1:
        
        plt.figure('Calc.Intensity_%i'%m)
        plt.gcf().set_size_inches(18,10)
        
        plt.subplot(4,4,9)
        plt.text(0, 0, infomation, fontsize=10, style='oblique')
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        
        for j in range(len(dx)*len(dy)):
            plt.subplot(4,4,j+1)
            plt.plot(D/Lamda, I_origin[j], color=fig_color[j], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    #         plt.xlabel('D/Lamda')
    #         plt.axis('scaled')
            plt.title('I_%i'%j)
    #     plt.title('Phase Delay = %f [deg]'%(k * np.pi/4 /(2*np.pi) * 360))
    
        for j in range(len(dx)*len(dy)-1):
            plt.subplot(4,4,j+10)
            plt.plot(I_origin[0], I_origin[j+1], color=fig_color[j+1], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
            plt.contour(fit_xaxis, fit_yaxis, fit_ellipse[j], [1])     #x**2 + y**2 = 9 的圆形
            plt.axis('scaled')
    #         plt.title('I_%i/I_0'%(j+1))
        
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
        plt.get_current_fig_manager().window.setGeometry(20, 50, 1600, 1000)
        plt.savefig('C:/Users/yu03/Desktop/Gabor/Interferometer Sim/Simulation_Result/2. V_y_2 Modulation/Calc.Intensity_%i.jpg'%m, dpi=300)
        plt.close()
        
        
        plt.figure('Fit_Result_%i'%m)
        plt.gcf().set_size_inches(18,10)
        
        plt.subplot(2,4,1)
        plt.text(0, 0, infomation, fontsize=10, style='oblique')
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')

        for j in range(len(dx)*len(dy)-1):
            plt.subplot(2,4,j+2)
    #         plt.plot(I_origin[0], I_origin[j+1], color=fig_color[j+1], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
            plt.contour(fit_xaxis, fit_yaxis, fit_ellipse[j], [1], colors=fig_color[j+1])     #x**2 + y**2 = 9 的圆形
            plt.axis('scaled')
            plt.text(0.5, 1, text[j], fontsize=10, style='oblique')
            plt.title('I_%i/I_0'%(j+1))

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
        plt.get_current_fig_manager().window.setGeometry(20, 50, 1800, 800)
        plt.savefig('C:/Users/yu03/Desktop/Gabor/Interferometer Sim/Simulation_Result/2. V_y_2 Modulation/Fit_Result_%i.jpg'%m, dpi=300)
        plt.close()
        
    #     plt.show()


plt.figure('Semi-axis length')
plt.gcf().set_size_inches(18,10)

plt.subplot(2,4,1)
plt.text(0, 0, 'V_y_2 from 0 to 0.000015', fontsize=10, style='oblique')
plt.xticks([])
plt.yticks([])
plt.axis('off')

for i in range(len(dx)*len(dy)-1):
    plt.subplot(2,4,i+2)
    slope_a = (a_set[i][-1] - a_set[i][0]) / (V_y_2_set[-1] - V_y_2_set[0])
    slope_b = (b_set[i][-1] - b_set[i][0]) / (V_y_2_set[-1] - V_y_2_set[0])
    plt.plot(V_y_2_set, a_set[i], label='semi-major %f'%slope_a, color='red', linestyle='solid', linewidth=1, marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(V_y_2_set, b_set[i], label='semi-minior %f'%slope_b, color='blue', linestyle='solid', linewidth=1, marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('V_y_2')
    plt.legend(loc='best')
plt.savefig('C:/Users/yu03/Desktop/Gabor/Interferometer Sim/Simulation_Result/2. V_y_2 Modulation/Compare_VmeaY.jpg', dpi=300)
plt.close()


    