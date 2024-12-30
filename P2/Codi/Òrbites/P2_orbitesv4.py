#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 20:01:29 2024

@author: davidrosamolina
"""

#Càlcul òrbites
import numpy as np
import matplotlib.pyplot as plt

#constants normalització-------------------------------------------------------
t_0=3600
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
r=np.zeros((pasos,2)) #dos dimensions
v=np.zeros((pasos,2))
t=np.zeros(pasos)




#CI

r_0=np.array([-2.546e10/d_0,1.449e11/d_0])
v_0=np.array([-2.982e4*t_0/d_0,-5.280e3*t_0/d_0])

r[0]=r_0
v[0]=v_0

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

#for i in range(1, pasos):
    #r[i], v[i] = runge_kutta(r[i-1], v[i-1], dt)


#Euler-------------------------------------------------------------------------

for j in range(1,pasos):
    t[j]=t[j-1]+dt
    a=acceleracio(r[j-1])
    v[j]= v[j-1] + a*dt #component x
    r[j]=r[j-1]+v[j]*dt

   




#plot--------------------------------------------------------------------------
plt.figure(figsize=(8, 8))
plt.plot(r[:,0]*d_0, r[:,1]*d_0, label="Òrbita de la Terra")
plt.scatter(0, 0, color='orange', label="Sol")  # el sol en aquest SR està al mig
plt.xlabel("x ")
plt.ylabel("y ")
plt.title("Òrbita de la Terra alr voltant del Sol (Euler)")
plt.legend()
plt.axis('equal')

plt.show() 



