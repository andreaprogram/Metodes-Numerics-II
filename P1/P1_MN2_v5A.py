# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:08:18 2024

@author: andre
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 08:48:19 2024

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt



# Dades---------
T_c=36.5+273.15 #temperatura inicial i de contorn al cos
t_a = 0.025
P_ext=944000
k=0.56
l_0=0.02
T_0=P_ext*(l_0**2)/k #674

#SOLUCIÓ ANALÍTICA----------------------------
N_fourier= int(np.power(10,4))
def T_anal(x,t):
    sumatori= sum( (4*(1-np.exp(-((2*n-1)**2)*(np.pi**2)*t))*np.sin((2*n-1)*np.pi*x)) / ((2*n-1)*np.pi)**3 for n in range(1,N_fourier))
    
    return (T_c/T_0)+((sumatori))

x_dib=np.linspace(0,1,101)


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


  Texp[0, :] = T_c/T_0   # CC: T(t,0)=T_c
  Texp[-1, :] = T_c/T_0  # CC: T(t,1)=T_c
  Texp[:,0] = T_c/T_0    # CI: T(0,x)=T_c


  for n in range(0,p_t-1,1): #itera sobre el temps
    for j in range(1,100,1): #itera sobre la posicio
         Texp[j,n+1]= g*(Texp[j+1,n]-2*Texp[j,n]+Texp[j-1,n]) + g*(Dx**2) + Texp[j,n]

  sol_exp.append((Texp[:,-1], g))


for i in range(3):
    
    plt.title("Euler explícit")
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('x/L')
    plt.ylabel('T(x,ta)')
    plt.plot(x_dib, (sol_exp[i][0])*T_0, label=f"Sol. numèrica $\\gamma = {sol_exp[i][1]}$")
    plt.plot(x_dib, (T_anal(x_dib,t_a))*T_0, label="Sol. Analítica")

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


import numpy as np

# Parámetros
Dx = 0.01
Dt = [Dx, 0.5 * Dx]

sol_imp = []
error_imp_list=[]

for t in Dt:
    p_x = 101
    p_t = int(round(t_a / t)) + 1
    gamma = t / (Dx ** 2)

    # Inicializar matriz de temperaturas
    Timp = np.zeros((p_x, p_t))
    Timp[0, :] = T_c/T_0   # CC: T(t,0) = T_c
    Timp[-1, :] = T_c / T_0 # CC: T(t,1) = T_c
    Timp[:, 0] = T_c / T_0  # CI: T(0,x) = T_c

    # Iterar en el tiempo resolviendo sistemas tridiagonales
    for n in range(1, p_t):
        # Matriz tridiagonal y vector independiente
        A = np.zeros((p_x, p_x))
        b = np.zeros(p_x)

        # Rellenar matriz A y vector b
        for j in range(1, p_x - 1):
            A[j, j - 1] = -gamma
            A[j, j] = 1 + 2 * gamma
            A[j, j + 1] = -gamma
            b[j] = Timp[j, n - 1] + t  # Incluye término Δt

        # Condiciones de contorno
        A[0, 0] = A[-1, -1] = 1
        b[0] = b[-1] = T_c/T_0

        # Resolver sistema
        Timp[:, n] = np.linalg.solve(A, b)

    # Guardar resultados en 
    sol_imp.append((Timp[:, -1], t))


for i in range(2):
    
    plt.title("Euler implícit")
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('x/L')
    plt.ylabel('T(x,ta)')
    plt.plot(x_dib, (sol_imp[i][0])*T_0, label=f"Sol. numèrica $\\Delta t= {sol_imp[i][1]}$")

    #Càlcul d'errors-----------
    error_imp=np.abs((sol_imp[i][0])-(T_anal(x_dib,t_a)))
    error_imp_list.append(error_imp)
    

plt.plot(x_dib, (T_anal(x_dib,t_a))*T_0, label="Sol. Analítica")
plt.legend()
plt.show()


#Plot error-----------------
plt.title("Errors Euler implícit")
plt.xlabel("x/L")
plt.ylabel("Error absolut")
for j in range(2):
    plt.plot(x_dib,error_imp_list[j], label=f"$\\Delta t= {sol_imp[j][1]}$")


plt.legend()
plt.show()