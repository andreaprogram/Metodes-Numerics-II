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




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Constants de normalització
def definir_constants():
    t_0 = 3600  # temps d'una hora en segons
    d_0 = 1.496e11  # distància mitja T-S en metres
    M_sol = 1.98840e30  # massa solar
    G_0 = np.power(d_0, 3) / ((t_0**2) * M_sol)
    G = 6.6743e-11 / G_0  # constant gravitacional normalitzada
    T_0 = 365 * 24 * 60 * 60  # període d'un any en segons
    return t_0, d_0, G, T_0 / t_0

t_0, d_0, G, T = definir_constants()
dt = 1  # pas temporal normalitzat (1 hora)

# Acceleració gravitatòria
def acceleracio(r):
    return -G * np.array(r) / (np.linalg.norm(r)**3)

# Mètode de Runge-Kutta (RK4)
def runge_kutta(r, v, dt):
    a1 = acceleracio(r)
    k1_r, k1_v = v, a1
    
    r2, v2 = r + 0.5 * k1_r * dt, v + 0.5 * k1_v * dt
    a2 = acceleracio(r2)
    k2_r, k2_v = v2, a2

    r3, v3 = r + 0.5 * k2_r * dt, v + 0.5 * k2_v * dt
    a3 = acceleracio(r3)
    k3_r, k3_v = v3, a3

    r4, v4 = r + k3_r * dt, v + k3_v * dt
    a4 = acceleracio(r4)
    k4_r, k4_v = v4, a4

    r_pas = r + (dt / 6) * (k1_r + 2*k2_r + 2*k3_r + k4_r)
    v_pas = v + (dt / 6) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
    return r_pas, v_pas

# Càlcul de trajectòria
def calcular_trajectoria(metode):
    passos = int(T / dt)
    r = np.zeros((passos, 2))
    v = np.zeros((passos, 2))
    r[0] = [-2.546e10 / d_0, 1.449e11 / d_0]
    v[0] = [-2.982e4 * t_0 / d_0, -5.280e3 * t_0 / d_0]

    for i in range(1, passos):
        if metode == "RK4":
            r[i], v[i] = runge_kutta(r[i-1], v[i-1], dt)
        elif metode == "Euler":
            a = acceleracio(r[i-1])
            v[i] = v[i-1] + a * dt
            r[i] = r[i-1] + v[i] * dt

    return r

# Dibuix trajectòries
def plot_trajectories(r_e, r_k):
    plt.figure(figsize=(8, 8))
    plt.plot(r_e[:, 0]*d_0, r_e[:, 1]*d_0, label="Euler")
    plt.plot(r_k[:, 0]*d_0, r_k[:, 1]*d_0, label="RK-4")
    plt.scatter(0, 0, s=280, color='orange', label="Sol")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.legend()
    
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    ax = plt.gca()
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    
    plt.show()

# Main: trajectòries
r_e = calcular_trajectoria("Euler")
r_k = calcular_trajectoria("RK4")
plot_trajectories(r_e, r_k)

# Trajectòria en 3D
def plot_trajectory_3D(rx, ry, rz):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(-rx, -ry, -rz, label="Trajectòria del Sol", linewidth=1)
    ax.scatter(0, 0, 0, s=200, color='green', label="Bigues i Riells del Fai")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.legend(loc='upper right')
    plt.show()

# Simulació 3D
passos = len(r_e)
rx = r_e[:, 0] * d_0
ry = r_e[:, 1] * d_0
rz = np.zeros(passos)  # Aproximació: pla XY
plot_trajectory_3D(rx, ry, rz)
