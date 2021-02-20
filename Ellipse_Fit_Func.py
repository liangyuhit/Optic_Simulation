# -*- coding: utf-8 -*-

#http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
#https://blog.csdn.net/liucc09/article/details/82763635

import numpy as np
from numpy.linalg import eig, inv

class Ellipse_fit():
    
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def fitEllipse(self):
        self.x = self.x[:,np.newaxis]
        self.y = self.y[:,np.newaxis]
        self.D =  np.hstack((self.x*self.x, self.x*self.y, self.y*self.y, self.x, self.y, np.ones_like(self.x)))
        self.S = np.dot(self.D.T, self.D)
        self.C = np.zeros([6,6])
        self.C[0,2] = self.C[2,0] = 2; self.C[1,1] = -1
        self.E, self.V =  eig(np.dot(inv(self.S), self.C))
        self.n = np.argmax(np.abs(self.E))
        self.a = self.V[:, self.n].real
        return np.array(self.a)
     
    def ellipse_get_parameters(self):
        self.fitEllipse()
        self.b, self.c, self.d, self.f, self.g, self.a = self.a[1]/2, self.a[2], self.a[3]/2, self.a[4]/2, self.a[5], self.a[0]
        self.num = self.b * self.b - self.a * self.c
        self.x0 = (self.c * self.d - self.b * self.f) / self.num
        self.y0 = (self.a * self.f - self.b * self.d) / self.num
        
        self.angle = 0.5 * np.arctan(2 * self.b / (self.a - self.c))
        
#         if self.b == 0:
#             if self.a > self.c:
#                 self.angle = 0
#             else:
#                 self.angle = np.pi/2
#         else:
#             if self.a > self.c:
#                 self.angle = np.arctan(2*self.b/(self.a-self.c))/2
#             else:
#                 self.angle = np.pi/2 + np.arctan(2*self.b/(self.a-self.c))/2
        
        self.up = 2 * (self.a * self.f * self.f + self.c * self.d * self.d + self.g * self.b * self.b - 2 * self.b * self.d * self.f - self.a * self.c * self.g)
        self.down1 = (self.b * self.b - self.a * self.c) * ((self.c - self.a) * np.sqrt(1 + 4 * self.b * self.b / ((self.a - self.c) * (self.a - self.c))) - (self.c + self.a))
        self.down2 = (self.b * self.b - self.a * self.c) * ((self.a - self.c) * np.sqrt(1 + 4 * self.b * self.b / ((self.a - self.c) * (self.a - self.c))) - (self.c + self.a))
        self.res1 = np.sqrt(np.abs(self.up / self.down1))
        self.res2 = np.sqrt(np.abs(self.up / self.down2))
#         self.semimajor = np.max(np.array([self.res1, self.res2]))
#         self.semiminior = np.min(np.array([self.res1, self.res2]))
        self.semimajor = self.res1
        self.semiminior = self.res2
        
        return np.array([self.x0, self.y0, self.semimajor, self.semiminior, self.angle])
     
#     def ellipse_angle_of_rotation( a ):
#         b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
#         return 0.5*np.arctan(2*b/(a-c))
#     
#     def ellipse_axis_length( a ):
#         b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
#         up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
#         down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
#         down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
#         res1=np.sqrt(up/down1)
#         res2=np.sqrt(up/down2)
#         return np.array([res1, res2])
#     
#     def ellipse_angle_of_rotation2( a ):
#         b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
#         if b == 0:
#             if a > c:
#                 return 0
#             else:
#                 return np.pi/2
#         else:
#             if a > c:
#                 return np.arctan(2*b/(a-c))/2
#             else:
#                 return np.pi/2 + np.arctan(2*b/(a-c))/2


    def ellipse_new_parameters(self):
        self.fitEllipse()
        self.b, self.c, self.d, self.f, self.g, self.a = self.a[1]/2, self.a[2], self.a[3]/2, self.a[4]/2, self.a[5], self.a[0]
        return (self.a, self.b, self.c, self.d, self.f, self.g)
        