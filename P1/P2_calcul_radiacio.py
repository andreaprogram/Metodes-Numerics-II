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
P=[]
for n in range(len(modul)):
    P.append(A*I_0*(R_S/(modul[n]))**2*eta*np.maximum(0,np.cos(angle_incident[n]*np.pi/180))*N)

print(P)
t = np.arange(24)  # t es el tiempo en horas (de 0 a 23)

P_24 = P[4320:4344]
plt.figure(figsize=(10, 6))
plt.plot(t, P_24, marker='o', color='b', linestyle='-', label="Potencia (W)")
plt.title("Potencia solar en función de la hora del día")
plt.xlabel("Hora del día (t) [h]")
plt.ylabel("Potencia (W)")
plt.grid(True)
plt.legend()
plt.show()