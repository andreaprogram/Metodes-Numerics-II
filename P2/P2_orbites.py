#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 20:01:29 2024

@author: davidrosamolina
"""

#Càlcul òrbites
import numpy as np
import matplotlib.pyplot as plt
#constants normalització
t_0=3600
d_0=1.496*np.power(10,11)
M_sol=1.98840*np.power(10,30)
M=1
G_0=np.power(d_0,3)/((t_0**2)*M_sol) 
G=6.6743/(np.power(10,11)*G_0)
T_0=365*24*60*60

T=T_0/t_0



r_0=np.array([147*np.power(10,9)/d_0,0])
v_0=np.array([0,29800*t_0/d_0])

dt= 1


def acceleracio(r):
   return  -G*M*r/np.linalg.norm(r)
    

pasos=int(T/dt)
#definim quantitats
r=np.zeros((pasos,2)) #dos dimensions
v=np.zeros((pasos,2))
t=np.zeros(pasos)
a    
#CI
r[0]=r_0
v[0]=v_0
t[0]=0
#Euler-----
for i in range (1,pasos):
    a=acceleracio(r[i-1])
    v[i]=v[i-1]+a*dt
    r[i]=r[i-1]+v[i-1]*dt
    t[i]=t[i-1]+dt




#graficar

plt.figure(figsize=(8, 8))
plt.plot(r[:, 0], r[:, 1], label="Órbita de la Tierra")
plt.scatter(0, 0, color='orange', label="Sol")  # Posición del Sol
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Órbita de la Tierra alrededor del Sol (Euler)")
plt.legend()
plt.axis('equal')
plt.show() 
