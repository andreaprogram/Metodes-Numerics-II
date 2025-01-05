#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 13:49:00 2024

@author: laurarodriguez
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants
latitud = 41 + 39 / 60 + 55.115 / 3600  # Latitud de Bigues i Riells en graus
latitud_rad = np.radians(latitud)       # Latitud en radians
gamma = np.radians(0)                   # Panell orientat cap al sud
dies_any = 365                          # Dies en un any
hores_dia = 24                          # Hores en un dia
I_0 = 6.320069973e7                     # Irradiància solar constant
R_S = 6.957e8                           # Radi del Sol
eta = 0.2                               # Eficiència del panell
A = 2                                   # Àrea d'un panell solar (m^2)
N1 = 1                                # Nombre de panells 


# Carregar les dades de coordenades sol-terra
data = np.loadtxt('coordenades_bigues_sol.txt', comments='#', delimiter=' ')
rx, ry, rz = data[:, 0], data[:, 1], data[:, 2]
modul = np.sqrt(rx**2 + ry**2 + rz**2)

# Funció per calcular la declinació solar
def declinacio(dia):
    return np.radians(23.45 * np.sin(np.radians(360 * (284 + dia) / 365)))

# Funció per calcular l'angle horari
def angle_horari(hora):
    return np.radians(15 * (hora - 12))

# Funció per calcular l'angle d'elevació solar
def elevacio_solar(delta, omega):
    sin_elevacio = np.sin(latitud_rad) * np.sin(delta) + np.cos(latitud_rad) * np.cos(delta) * np.cos(omega)
    return np.arcsin(sin_elevacio)  # En radians

# Funció per calcular el cos de l'angle d'incidència
def cos_theta_i(delta, omega, beta, gamma):
    term1 = np.sin(delta) * np.sin(latitud_rad) * np.cos(beta)
    term2 = -np.sin(delta) * np.cos(latitud_rad) * np.sin(beta) * np.cos(gamma)
    term3 = np.cos(delta) * np.cos(latitud_rad) * np.cos(beta) * np.cos(omega)
    term4 = np.cos(delta) * np.sin(latitud_rad) * np.sin(beta) * np.cos(gamma) * np.cos(omega)
    term5 = np.cos(delta) * np.sin(beta) * np.sin(gamma) * np.sin(omega)
    return term1 + term2 + term3 + term4 + term5

# Funció per calcular la potència
def calcular_potencia(angle_incident, N):
    P = np.zeros_like(angle_incident)
    for dia in range(P.shape[0]):
        for hora in range(P.shape[1]):
            cos_max = np.maximum(0, np.cos(angle_incident[dia, hora] * np.pi / 180))
            P[dia, hora] = A * I_0 * (R_S / modul[dia * 24 + hora])**2 * eta * cos_max * N
    return P

# Vector per provar diferents inclinacions (en graus)
betas = np.arange(38, 39, 0.01)  
potencies_totals = []

h=  815   #altura muntanya
L= 4.68e3    #distància muntanya
sigma=np.arctan(h/L)

for beta_deg in betas:
    beta = np.radians(beta_deg)
    angles_incidents = []

    # Simulació al llarg de l'any
    for dia in range(1, dies_any + 1):
        delta = declinacio(dia)
        angles_dia = []
        for hora in range(hores_dia):
            omega = angle_horari(hora)
            elevacio = elevacio_solar(delta, omega)
            if elevacio > sigma:  # Només si el Sol està per sobre de la muntanya
                cos_theta = cos_theta_i(delta, omega, beta, gamma)
                theta_i = np.degrees(np.arccos(cos_theta)) if -1 <= cos_theta <= 1 else np.nan
            else:
                theta_i = np.nan  # Sol sota l'horitzó
            angles_dia.append(theta_i)
        angles_incidents.append(angles_dia)

    angles_incidents = np.array(angles_incidents)
    angles_incidents = np.nan_to_num(angles_incidents, nan=90)

    # Calcul de potència per un panell
    P_T = calcular_potencia(angles_incidents, N1)

    # Potència total combinada
    
    potencies_totals.append(np.sum(P_T))

# Determinar la inclinació òptima
potencies_totals = np.array(potencies_totals)
beta_optima = betas[np.argmax(potencies_totals)]
potencia_maxima = potencies_totals.max()

print(f"La inclinació òptima (\u03B2) és: {beta_optima} graus")
print(f"La potència màxima anual és: {potencia_maxima:.2f} W")

# Representació gràfica
plt.figure(figsize=(10, 6))
plt.plot(betas, potencies_totals, marker='o')
plt.title("Potència total anual en funció de la inclinació (\u03B2)")
plt.xlabel("Inclinació (\u03B2) [graus]")
plt.ylabel("Potència total anual [W]")
plt.grid(True)
plt.show()
