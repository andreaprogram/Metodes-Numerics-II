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
t_a = 0.025 #temps al que volem trobar el perfil de temperatures
P_ext=944000
k=0.56
rho=1081
c_v=3686
l_0=0.02
T_0=P_ext*(l_0**2)/k #674 , Temperatura de normalització
t_0=(l_0**2)*rho*c_v/k #constant normalització temporal

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

#EULER EXPLÍCIT -------------------------------------------------------------------------------------

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


#Plotting-----------------------
#Colors triats:
    # gamma = 0.51 blau
    # gamma = 0.49 taronja
    # gamma = 0.25 verd

color_list = ["blue", "orange", "green", "magenta"]
for i in range(3):
    
    plt.title("Euler explícit")
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('x (cm)')
    plt.ylabel('T(x,ta) (K)')
    plt.plot(x_dib*2, (sol_exp[i][0])*T_0, label=f"Sol. numèrica $\\gamma = {sol_exp[i][1]}$", color=color_list[i])
    plt.plot(x_dib*2, (T_anal(x_dib,t_a))*T_0, label="Sol. Analítica", color="black")
    plt.axvline(0.75, linestyle='--', color="red", label="Teixit malalt")
    plt.axvline(1.25, linestyle='--', color="red")
    plt.axhline(323.15, linestyle='--', color="purple", label="50°C")
    plt.axvspan(0.75, 1.25, color='red', alpha=0.2)
    plt.legend(loc="upper right")
    plt.show()
    
    #Càlcul d'errors-----------
    error_exp=np.abs(sol_exp[i][0]-T_anal(x_dib,t_a))
    error_exp_list.append(error_exp)

#Plot error-----------------
plt.title("Errors Euler explícit")
plt.tick_params(axis='both', direction='in')
plt.xlabel("x(cm)")
plt.ylabel("Error absolut (K)")
plt.plot(x_dib*2,(error_exp_list[0])*T_0, label=f"$\\gamma = {sol_exp[0][1]}$", color=color_list[0])
plt.axvline(0.75, linestyle='--', color="red", label="Teixit malalt")
plt.axvline(1.25, linestyle='--', color="red")
plt.axvspan(0.75, 1.25, color='red', alpha=0.2)
plt.legend(loc="upper right")
plt.show()


for j in range(1,3,1):
    plt.plot(x_dib*2,(error_exp_list[j])*T_0, label=f"$\\gamma = {sol_exp[j][1]}$", color=color_list[j])

plt.axvline(0.75, linestyle='--', color="red", label="Teixit malalt")
plt.axvspan(0.75, 1.25, color='red', alpha=0.2)
plt.axvline(1.25, linestyle='--', color="red")
plt.title("Errors Euler explícit")
plt.tick_params(axis='both', direction='in')
plt.xlabel("x(cm)")
plt.ylabel("Error absolut (K)")
plt.legend(loc="upper right")
plt.show()



#EULER IMPLÍCIT--------------------------------------------------------------------------------------------
Dx = 0.01
Dt = [Dx, 0.5 * Dx]
t_a = 0.03

sol_imp = []
error_imp_list=[]

for t in Dt:
    p_x = 101
    p_t = int(round(t_a / t)) + 1
    gamma = t / (Dx ** 2)

    #Definim la nostra matriu de temperatures T_j^n
    Timp = np.zeros((p_x, p_t))
    Timp[0, :] = T_c/T_0   # CC: T(t,0) = T_c
    Timp[-1, :] = T_c / T_0 # CC: T(t,1) = T_c
    Timp[:, 0] = T_c / T_0  # CI: T(0,x) = T_c

    #Resolem el sistema d'equacions amb matrius
    #aquest proces ens donara el vector columna Tj^n pel qual després haurem de resoldre pel 
    for n in range(1, p_t):
        

            # MATRIU A
            A = np.zeros((p_x-2, p_x-2))    
        
            #Diagonals d'A
            upper_diag = (-1)*gamma*np.ones(p_x-3)
            diag = (1+2*gamma)*np.ones(p_x-2)
            lower_diag = (-1)*gamma*np.ones(p_x-3)

            A = np.diag(diag) + np.diag(lower_diag, -1) + np.diag(upper_diag, 1)
        
            #MATRIU b
            b = Timp[1:-1, n - 1] + t #b escull T_j^n-1 de j=1,...,N-1 
            b[0]+= gamma*Timp[0,n] #afegim les condicions de contorn al vector b
            b[-1] += gamma*Timp[-1,n]
        
            
            Timp[1:-1,n]= np.linalg.solve(A,b)
            
    sol_imp.append((Timp[:,-1], t))


#Plotting------------------------------
for i in range(2):
    
    plt.title("Euler implícit")
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('x(cm)')
    plt.ylabel('T(x,ta) (K)')
    plt.plot(x_dib*2, (sol_imp[i][0])*T_0, label=f"Sol. numèrica $\\Delta t= {sol_imp[i][1]}$", color=color_list[i])
    
    #Càlcul d'errors-----------
    error_imp=np.abs((sol_imp[i][0])-(T_anal(x_dib,t_a)))
    error_imp_list.append(error_imp)

plt.plot(x_dib*2, (T_anal(x_dib,t_a))*T_0, label="Sol. Analítica", color='black')
plt.axvline(0.75, linestyle='--', color="red", label="Teixit malalt")
plt.axvline(1.25, linestyle='--', color="red")
plt.axvspan(0.75, 1.25, color='red', alpha=0.2)
plt.axhline(323.15, linestyle='--', color="purple", label="50°C")
plt.legend(loc="lower right")
plt.show()


#Plot error-----------------
plt.title("Errors Euler implícit")
plt.xlabel("x(cm)")
plt.ylabel("Error absolut (K)")
for j in range(2):
    plt.plot(x_dib*2, (error_imp_list[j])*T_0, label=f"$\\Delta t= {sol_imp[j][1]}$", color=color_list[j])
    
plt.axvline(0.75, linestyle='--', color="red", label="Teixit malalt")
plt.axvline(1.25, linestyle='--', color="red")  
plt.axvspan(0.75, 1.25, color='red', alpha=0.2)  
plt.legend(loc="upper right")
plt.show()


#CRANK - NICOLSON ---------------------------------------------------------------------------------------
Dx = 0.01
Dt = [0.5*Dx, 0.25*Dx]
t_a=0.025

sol_cn = []
error_cn_list=[]

for t in Dt:
    p_x = 101
    p_t = int(round(t_a / t)) + 1
    gamma = (t / (Dx ** 2))
    
    #Ara, resoldrem com l'apartat anterior pero amb diferents A i b

    #Definim la nostra matriu de temperatures T_j^n+1
    Tcn = np.zeros((p_x, p_t))
    Tcn[0, :] = T_c/T_0   # CC: T(t,0) = T_c
    Tcn[-1, :] = T_c / T_0 # CC: T(t,1) = T_c
    Tcn[:, 0] = T_c / T_0  # CI: T(0,x) = T_c

    
    for n in range(1, p_t):
        

            # MATRIU A
            #És la mateixa que en EULEr IMPLÍCIT
            #Diagonals d'A
            upper_diag = (-0.5)*gamma*np.ones(p_x-3)
            diag = (1+gamma)*np.ones(p_x-2)
            lower_diag = (-0.5)*gamma*np.ones(p_x-3)

            A = np.diag(diag) + np.diag(lower_diag, -1) + np.diag(upper_diag, 1)
        
            #MATRIU b
            #   b = B*Tj^n-1 + terme CC + terme t
            
            diag_B = (1 - gamma) * np.ones(p_x - 2)
            upper_diag_B =  (0.5)*gamma * np.ones(p_x - 3)
            B = np.diag(diag_B) + np.diag(upper_diag_B, 1) + np.diag(upper_diag_B, -1)
            
            # Tj^n-1
            Tant = Tcn[1:-1, n-1] 
            
            # -A*Tj^n-1
            b = np.dot(Tant,B)
            
            #Terme t
            b += t*np.ones_like(b)
           
            #Terme CC
            b[0] += gamma* Tcn[0, n]
            b[-1] += gamma* Tcn[-1, n]
            
            Tcn[1:-1,n]= np.linalg.solve(A,b)
            
    sol_cn.append((Tcn[:,-1], t))


#Plotting------------------------------
for i in range(2):
    
    plt.title("Crank-Nicolson")
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('x(cm)')
    plt.ylabel('T(x,ta) (K)')
    plt.plot(x_dib*2, (sol_cn[i][0])*T_0, label=f"Sol. numèrica $\\Delta t = {sol_cn[i][1]}$", color=color_list[i])
    
    #Càlcul d'errors-----------
    error_cn=np.abs((sol_cn[i][0])-(T_anal(x_dib,t_a)))
    error_cn_list.append(error_cn)

plt.plot(x_dib*2, (T_anal(x_dib,t_a))*T_0, label="Sol. Analítica", color='black')
plt.axvline(0.75, linestyle='--', color="red", label="Teixit malalt")
plt.axvline(1.25, linestyle='--', color="red")
plt.axhline(323.15, linestyle='--', color="purple", label="50°C")
plt.axvspan(0.75, 1.25, color='red', alpha=0.2)
plt.legend(loc="lower right")
plt.show()


#Plot error-----------------
plt.title("Errors Crank-Nicolson")
plt.xlabel("x(cm)")
plt.ylabel("Error absolut (K)")
for j in range(2):
    plt.plot(x_dib*2, (error_cn_list[j])*T_0, label=f"$\\Delta t = {sol_cn[j][1]}$", color=color_list[j])
    
plt.axvline(0.75, linestyle='--', color="red", label="Teixit malalt")
plt.axvline(1.25, linestyle='--', color="red")    
plt.axvspan(0.75, 1.25, color='red', alpha=0.2)
plt.legend(loc="upper right")
plt.show()
