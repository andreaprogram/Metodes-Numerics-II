#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 09:31:06 2024

@author: davidrosamolina
"""


import numpy as np
import matplotlib.pyplot as plt

# Dades---------
T_c=36.5+273.15 #temperatura inicial i de contorn al cos
t_a = 0.025 #temps al que volem trobar el perfil de temperatures
P_ext=944000
k=0.56
rho=1081
c_v=3686
l_0=0.02
T_0=P_ext*(l_0**2)/k #674 , Temperatura de normalització
t_0=(l_0**2)*rho*c_v/k #constant normalització temporal

T_limit1=(50+273.15)/T_0
T_limit2=(80+273.15)/T_0

z=np.linspace(0,1,101)

def condicions (T,x,t):
    zona_malalta=(x>=(0.75/2))&(x<=(1.25/2))
    zona_sana= ~zona_malalta
    
    if np.all(T[zona_malalta]>T_limit1) and np.max(T)<T_limit2:
        if np.all(T[zona_sana]<T_limit1):
            return True
        else:
            return False
    else:
        return False
    
t_opt=None

Dx=0.01
Dt=0.25*(Dx**2)
Ta=np.loadtxt("matriu_euler_exp_0.25.txt")
T_final=Ta[:,-1]


for s in range(Ta.shape[1]):
    columna=Ta[:,s]/T_0 #agafem la columna per un temps s de la matriu
    t=s*Dt
    #veiem si es compleixen condicions
    if condicions (columna,z,t):
        t_opt= t
        
        break
        
incertesa= 0.25*(Dx**2)*t_0

if t_opt is not None:
    print(f"El temps òptim per les condicions del problema és {t_opt*t_0:.3f} ± {incertesa:.3f} s", "amb gamma=", 0.25)
    col_opt=int(t_opt/Dt)
    T_opt=Ta[:,col_opt]
    plt.title("Perfil de temperatures per $t_{opt}$")
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('x (cm)')
    plt.ylabel('T(x,$t_{opt}$) (K)')
    plt.plot(z*2, T_opt, label="Solució numèrica per $\gamma=0.25$")
    plt.axvline(0.75, linestyle='--', color="red", label="Teixit malalt")
    plt.axvline(1.25, linestyle='--', color="red")
    plt.axhline(323.15, linestyle='--', color="purple", label="50°C")
    plt.axvspan(0.75, 1.25, color='red', alpha=0.2)
    plt.legend(loc="center right")
    plt.show()
    

    
    
    
    
    
    
    