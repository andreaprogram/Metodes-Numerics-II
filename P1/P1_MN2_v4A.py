# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 08:48:19 2024

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Dades---------
T_c=36.5+273.15 #temperatura inicial i de contorn al cos
t_a = 0.025 

#SOLUCIÓ ANALÍTICA----------------------------
N_fourier= int(np.power(10,4))
def T_anal(x,t):
    sumatori= sum( (4*(1-np.exp(-((2*n-1)**2)*(np.pi**2)*t))*np.sin((2*n-1)*np.pi*x)) / ((2*n-1)*np.pi)**3 for n in range(1,N_fourier))
    
    return T_c+((sumatori))

x_dib=np.linspace(0,1,101)

ax = plt.gca()
ax.yaxis.get_offset_text().set_visible(False)  # Oculta el "+3.096e2"
plt.ticklabel_format(axis="y", style="plain")  # Fuerza valores completos
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))  # Números planos
plt.plot(x_dib, T_anal(x_dib,t_a))
plt.show()


# SOLUCIONS NUMÈRIQUES-----------------------------------------------------
#Nº punts maiats espaial
N=101

#Increments en x
Dx = 1/(101-1) #delta_x= 0.01

#EULER EXPLÍCIT --------------------------------

#j va de 1 a 101
# n comença al 1

#Gamma = deltat/(deltax)^2
gamma=[0.51,0.49,0.25]


sol_exp=[]
error_exp_list= []

for g in gamma:
  p_x=101
  Dt=g*(Dx**2)
  p_t=int(t_a/(Dt))+1

  Texp=np.zeros((p_x,p_t))


  Texp[0, :] = T_c   # CC: T(t,0)=T_c
  Texp[-1, :] = T_c  # CC: T(t,1)=T_c
  Texp[:,0] = T_c    # CI: T(0,x)=T_c


  for n in range(0,p_t-1,1): #itera sobre el temps
    for j in range(1,100,1): #itera sobre la posicio
         Texp[j,n+1]= g*(Texp[j+1,n]-2*Texp[j,n]+Texp[j-1,n]) + g*(Dx**2) + Texp[j,n]


  sol_exp.append((Texp[:,-1], g)*309)

    

for i in range(3):
    
    plt.title("Euler explícit")
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('x/L')
    plt.ylabel('T(x,ta)')
    plt.plot(x_dib, sol_exp[i][0], label=f"Sol. numèrica $\\gamma = {sol_exp[i][1]}$")
    plt.plot(x_dib, T_anal(x_dib,t_a), label="Sol. Analítica")
    plt.legend()
    plt.show()
    
    #Càlcul d'errors-----------
    error_exp=np.abs(sol_exp[i][0]-T_anal(x_dib,t_a))
    error_exp_list.append(error_exp)

#Plot error-----------------
plt.title("Errors Euler explícit")
plt.tick_params(axis='both', direction='in')
plt.xlabel("x/L")
plt.ylabel("Error absolut")
for j in range(3):
    plt.plot(x_dib,error_exp_list[j], label=f"$\\gamma = {sol_exp[j][1]}$")
plt.legend()
plt.show()


# EULER IMPLÍCIT----------------------------
Dx = 0.01
Dt = [Dx, 0.5*Dx]

sol_imp=[]
error_imp_list= []

for t in Dt:
  p_x=101
  p_t=int(t_a/(t))+1
  g=Dx/(t**2)

  Timp=np.zeros((p_x,p_t))


  Timp[0, :] = T_c   # CC: T(t,0)=T_c
  Timp[-1, :] = T_c  # CC: T(t,1)=T_c
  Timp[:,0] = T_c    # CI: T(0,x)=T_c


  for n in range(0,p_t-1,1): #itera sobre el temps
    for j in range(1,100,1): #itera sobre la posicio
         Timp[j,n]= ( g*(Timp[j+1,n] + Timp[j-1,n]) + Timp[j,n-1] + t) / (1+2*g)


  sol_imp.append((Timp[:,-1], t))

    

for i in range(2):
    
    plt.title("Euler implícit")
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('x/L')
    plt.ylabel('T(x,ta)')
    plt.plot(x_dib, sol_imp[i][0], label=f"Sol. numèrica $\\Delta t= {sol_imp[i][1]}$")
    plt.plot(x_dib, T_anal(x_dib,t_a), label="Sol. Analítica")
    plt.legend()
    plt.show()
    
    #Càlcul d'errors-----------
    error_imp=np.abs(sol_imp[i][0]-T_anal(x_dib,t_a))
    error_imp_list.append(error_imp)

#Plot error-----------------
plt.title("Errors Euler implícit")
plt.xlabel("x/L")
plt.ylabel("Error absolut")
for j in range(2):
    plt.plot(x_dib,error_imp_list[j], label=f"$\\Delta t= {sol_imp[j][1]}$")
plt.legend()
plt.show()