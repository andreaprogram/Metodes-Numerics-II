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
d_0=1.496e11
M_sol=1.98840*np.power(10,30)
M=1
G_0=np.power(d_0,3)/((t_0**2)*M_sol) 
G=6.6743e-11/(G_0)
T_0=365*24*60*60

T=T_0/t_0



r_0=np.array([2.531e7/d_0,1.45e8])
v_0=np.array([-2.982e1*t_0/d_0,4.997*t_0/d_0])

dt=1


def acceleracio(r):
   return  -G*M*np.array(r)/(np.linalg.norm(r)**3)


pasos= int(T/dt)
#definim quantitats
r=np.zeros((pasos,2)) #dos dimensions
v=np.zeros((pasos,2))
t=np.zeros(pasos)


#CI
r[0]=r_0
v[0]=v_0

#print("r",r)
#print("v",v)

t[0]=0
#Euler-----



plt.figure(figsize=(8, 8))
plt.scatter(0, 0, color='orange', label="Sol")  # Posición del Sol
plt.xlabel("x ")
plt.ylabel("y ")
plt.title("Órbita de la Tierra alrededor del Sol (Euler)")
plt.legend()
plt.axis('equal')

for j in range(1,pasos):
    t[j]=t[j-1]+dt
    a=acceleracio(r[j-1])
    v[j][0]= v[j-1][0] + a[0]*dt #component x
    v[j][1]= v[j-1][1] + a[1]*dt #component y
    
    r[j][0]=r[j-1][0]+v[j-1][0]*dt
    r[j][1]=r[j-1][1]+v[j-1][1]*dt
    
   # print(r[j][0])
    plt.scatter(r[j-1][0], r[j-1][1], label="Órbita de la Tierra")
    
    
    

   #for i in range (1,pasos):
       #  a=acceleracio(r[i-1])
       #  v[i]=v[i-1]+ax*dt
        # r[i]=r[i-1]+v[i-1]*dt
        # t[i]=t[i-1]+dt




#graficar

plt.show() 
