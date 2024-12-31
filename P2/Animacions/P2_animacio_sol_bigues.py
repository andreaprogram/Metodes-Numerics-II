#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 17:09:29 2024

@author: davidrosamolina
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


data = np.loadtxt('coordenades_bigues_sol.txt', comments='#', delimiter=' ')
rx, ry, rz = data[:, 0], data[:, 1], data[:, 2]


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

#límits eixos
max_val = np.max([np.max(rx), np.max(ry), np.max(rz)])
min_val = np.min([np.min(rx), np.min(ry), np.min(rz)])

ax.set_xlim([min(rx), max(rx)])
ax.set_ylim([min(ry), max(ry)])
ax.set_zlim([min_val, max_val])


ax.set_title('Trajectòria del Sol al voltant de Bigues i Riells del Fai')
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")

#inici animació
line, = ax.plot([], [], [], label="Trayectoria del Sol", color='orange', linewidth=1)
sun_marker = ax.scatter([], [], [], s=200, color='orange', label="Posició del Sol")  
bigues_marker = ax.scatter(0, 0, 0, s=100, color='#1E90FF', label="Bigues i Riells del Fai") 

# Hora
time_text = ax.text2D(0.05, 0.95, '', transform=ax.transAxes, fontsize=14, color='black')


def init():
    line.set_xdata([])  
    line.set_ydata([])  
    line.set_3d_properties([])  
    sun_marker._offsets3d = ([], [], [])  
    time_text.set_text('')  
    return line, sun_marker, time_text, bigues_marker


def update(frame):
    # Actualiza trajectoria
    line.set_xdata(rx[:frame])  
    line.set_ydata(ry[:frame])  
    line.set_3d_properties(rz[:frame])  #

 
    sun_marker._offsets3d = (rx[frame:frame+1], ry[frame:frame+1], rz[frame:frame+1])

    # Calcular hora
    hora = frame 
    horas = int(hora) 

    #text en temps
    time_text.set_text(f"Hora: {horas:02d} h")

    return line, sun_marker, time_text, bigues_marker


ani = FuncAnimation(fig, update, frames=len(rx), init_func=init, interval=20, repeat=False)

# Guardar len gif
ani.save("animacio_sol_bigues.gif", writer="imagemagick", fps=60)


plt.legend()
plt.show()
