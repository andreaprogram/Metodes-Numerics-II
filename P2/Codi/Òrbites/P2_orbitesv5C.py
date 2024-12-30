#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 20:01:29 2024

@author: davidrosamolina
"""

#Càlcul òrbites
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


#constants normalització-------------------------------------------------------
t_0=3600 #temps d'una hora en s
d_0= 1.496e11 #distància mitja T-S en m
M_sol=1.98840e30 #massa sol en kg
M=1
G_0=np.power(d_0,3)/((t_0**2)*M_sol) 
G=6.6743e-11/(G_0) #constant gravitació normalitzada
T_0=365*24*60*60  #periode 1 any en segons

T=T_0/t_0


dt=1 #hem normalitzat a hores, fem pasos d'1h


def acceleracio(r):
   return  -G*M*np.array(r)/(np.linalg.norm(r)**3)


pasos= int(T/dt)
#definim quantitats i condicions inicials--------------------------------------
r_e=np.zeros((pasos,2)) #dos dimensions
v_e=np.zeros((pasos,2))
t=np.zeros(pasos)

r_k=np.zeros((pasos,2)) #dos dimensions
v_k=np.zeros((pasos,2))
t_k=np.zeros(pasos)




# CONDICIONS INICIALS--------------------------------

r_0=np.array([-2.546e10/d_0,1.449e11/d_0])
v_0=np.array([-2.982e4*t_0/d_0,-5.280e3*t_0/d_0])

r_e[0]=r_0
v_e[0]=v_0
r_k[0]=r_0
v_k[0]=v_0

#print("r",r)
#print("v",v)

t[0]=0

#Runge Kutta 4-----------------------------------------------------------------
def runge_kutta(r,v,dt):
    a1=acceleracio(r) #primer pas
    k1_r=v
    k1_v=a1
    
    r2 = r + 0.5 * k1_r * dt  #segon pas
    v2 = v + 0.5 * k1_v * dt
    a2 = acceleracio(r2)
    k2_r = v2
    k2_v = a2
    
    r3 = r + 0.5 * k2_r * dt #tercer pas
    v3 = v + 0.5 * k2_v * dt
    a3 = acceleracio(r3)
    k3_r = v3
    k3_v = a3
    
    r4 = r + k3_r * dt  #quart pas
    v4 = v + k3_v * dt
    a4 = acceleracio(r4)
    k4_r = v4
    k4_v = a4
    
    r_pas = r + (dt / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r)
    v_pas = v + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)

    return r_pas, v_pas

for i in range(1, pasos):
    r_k[i], v_k[i] = runge_kutta(r_k[i-1], v_k[i-1], dt)


#Euler-------------------------------------------------------------------------

for j in range(1,pasos):
    t[j]=t[j-1]+dt
    a=acceleracio(r_e[j-1])
    v_e[j]= v_e[j-1] + a*dt #component x
    r_e[j]=r_e[j-1]+v_e[j]*dt



#plot--------------------------------------------------------------------------
plt.figure(figsize=(8, 8))
plt.plot(r_e[:,0]*d_0, r_e[:,1]*d_0, label="Euler")
plt.plot(r_k[:,0]*d_0, r_k[:,1]*d_0, label="RK-4")
plt.scatter(0, 0, s=280, color='orange', label="Sol")  # el sol en aquest SR està al mig
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.legend()


# Aplicar notació científica als eixos
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2, 2))

ax = plt.gca()  # Obtenir l'objecte Axes actual
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)


plt.show()


#fem la rotació d'eixos i canvi de sistema de referència

#angles i constants------------------------------------------------------------
alfa=23.4*2*np.pi/360

theta=41.68*np.pi/180
Rt=12.756e3/2
w= 7.272e-5*t_0 #rad/h
phi=0

#tot en radians per utilitzar np.funcions trigonometriques.

def R_bigues(t): #vector posició bigues pel centre de la terra amb eixos ja rotats
    x=Rt*(np.cos(theta)*np.cos(phi+w*t))
    y=Rt*(np.cos(theta)*np.sin(phi+w*t)*np.cos(alfa)+np.sin(theta)*np.sin(alfa))
    z=Rt*(-np.cos(theta)*np.sin(phi+w*t)*np.sin(alfa)+np.sin(theta)*np.cos(alfa))
    
    return x,y,z


#ara fem el vector posició de bigues pel centre del sol


rx = np.zeros(pasos)
ry = np.zeros(pasos)
rz = np.zeros(pasos)

#on es troba la terra quan t=0
rx[0]=np.array(-2.546e10+Rt*(np.cos(theta)*np.cos(phi)))
ry[0]=np.array(1.449e11+Rt*(np.cos(theta)*np.sin(phi)*np.cos(alfa)+np.sin(theta)*np.sin(alfa)))
rz[0]=np.array(Rt*(np.sin(theta)*np.cos(alfa)))

for i in range(1,pasos):
    bigues=R_bigues(i)

    rx[i]=r_e[i,0]*d_0+bigues[0]
    ry[i]=r_e[i,1]*d_0+bigues[1]
    rz[i]=bigues[2]
    
 

    
plt.figure(figsize=(8, 8))
plt.plot(rx, ry, label="Trajectòria de Bigues i Riells del Fai")

plt.scatter(0, 0, s=280, color='orange', label="Sol")  # el sol en aquest SR està al mig
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.legend(loc='upper right')



plt.show()





# Guardem coordenades
filename = 'coordenades_bigues_sol.txt'


data = np.column_stack((rx, ry, rz)) 
np.savetxt(filename, data, header="rx (m), ry (m), rz (m)", fmt='%.6e')



#volem veure el sol respecte habitatge

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Traza la trayectoria en 3D
ax.plot(-rx, -ry, -rz, label="Trajectòria del Sol", linewidth=1)

# Dibuja el sol (en el origen)
ax.scatter(0, 0, 0, s=200, color='green', label="Bigues i Riells del Fai")

max_val = np.max([np.max(rx), np.max(ry), np.max(rz)])
min_val = np.min([np.min(rx), np.min(ry), np.min(rz)])

# Ajustar los límites Z para que estén dentro del rango de X e Y
ax.set_xlim([min(rx), max(rx)])
ax.set_ylim([min(ry), max(ry)])
ax.set_zlim([min_val, max_val])  # Ajusta Z a un rango similar a X e Y

# Etiquetas
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")

# Leyenda
ax.legend(loc='upper right')

# Mostrar el gráfico
plt.show()
