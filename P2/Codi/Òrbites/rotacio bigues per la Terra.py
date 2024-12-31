#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 13:47:47 2024

@author: davidrosamolina
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 




#angles i constants------------------------------------------------------------
alfa=23.5*2*np.pi/360

theta=41.68*np.pi/180
Rt=1
w= 2*np.pi/(24*3600)  #tot en radians per utilitzar np.funcions trigonometriques.

def R_bigues(t):
    x=Rt*(np.cos(theta)*np.cos(w*t))
    y=Rt*(np.cos(theta)*np.sin(w*t)*np.cos(alfa)+np.sin(theta)*np.sin(alfa))
    z=Rt*(-np.cos(theta)*np.sin(w*t)*np.sin(alfa)+np.sin(theta)*np.cos(alfa))
    
    return x,y,z

t=np.linspace(0,24*3600,1000)

x,y,z=R_bigues(t)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

ax.quiver(0, 0, 0, 0, np.sin(alfa),np.cos(alfa), length=1, color='r', label="Eje Z de la Tierra (inclinación)")

ax.plot(x, y, z, label="Trayectoria de la rotación")


ax.set_title('Rotació de la Terra per la latitud de Bigues i Riells del Fai')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.view_init(elev=20, azim=40) 
ax.legend()
plt.show()

