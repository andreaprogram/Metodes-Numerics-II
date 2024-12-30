#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 23:37:19 2024

@author: davidrosamolina
"""

import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('coordenades_bigues_sol.txt', comments='#', delimiter=' ')
angle_incident=np.loadtxt('angles_incidents.txt', comments='#', delimiter=' ')
rx, ry, rz = data[:, 0], data[:, 1], data[:, 2]


modul = np.sqrt(rx**2 + ry**2 + rz**2)
I_0=1361
R_S=6.957e8
eta=0.2
A=2
N=5  #nombre plaques

P = np.zeros_like(angle_incident)
for dia in range(P.shape[0]):
    for hora in range(P.shape[1]):
        cos_max = np.maximum(0, np.cos(angle_incident[dia, hora] * np.pi / 180))
        P[dia, hora] = A * I_0 * (R_S / modul[dia * 24 + hora])**2 * eta * cos_max * N



dies = np.arange(1, P.shape[0] + 1)  
hores = np.arange(24) 
X, Y = np.meshgrid(dies, hores)  #mallat 
Z = P.T  # trsposa


fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, Y, Z, cmap='plasma', edgecolor='none')
fig.colorbar(surf, ax=ax, label="Potència generada(W)")
ax.set_title("Potència generada durant 1 any")
ax.set_xlabel("Dia")
ax.set_ylabel("Hora")
ax.set_zlabel("Potència generada (W)")
plt.show()
