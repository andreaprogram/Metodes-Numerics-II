import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Definim les constants
latitud = 41 + 39 / 60 + 55.115 / 3600  # Latitud de Bigues i Riells en graus
latitud_rad = np.radians(latitud)  # Latitud en radians
dies_any = 365  # Dies en un any
hores_dia = 24  # Hores en un dia

# Funció per calcular la declinació solar
def declinacio(dia):
    return np.radians(23.45 * np.sin(np.radians(360 * (284 + dia) / 365)))

# Funció per calcular l'angle horari
def angle_horari(hora):
    return np.radians(15 * (hora - 12))

# Funció per calcular el cosinus de l'angle d'incidència
def cos_theta_i(delta, omega, beta, gamma):
    term1 = np.sin(delta) * np.sin(latitud_rad) * np.cos(beta)
    term2 = -np.sin(delta) * np.cos(latitud_rad) * np.sin(beta) * np.cos(gamma)
    term3 = np.cos(delta) * np.cos(latitud_rad) * np.cos(beta) * np.cos(omega)
    term4 = np.cos(delta) * np.sin(latitud_rad) * np.sin(beta) * np.cos(gamma) * np.cos(omega)
    term5 = np.cos(delta) * np.sin(beta) * np.sin(gamma) * np.sin(omega)
    cos_theta = term1 + term2 + term3 + term4 + term5
    return np.clip(cos_theta, 0, 1)  # Només comptem la irradiació directa (cos > 0)

# Funció per calcular la potència total anual
def potencia_total(beta, gamma):
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)
    pot_total = 0
    for dia in range(1, dies_any + 1):
        delta = declinacio(dia)
        for hora in range(hores_dia):
            omega = angle_horari(hora)
            cos_theta = cos_theta_i(delta, omega, beta_rad, gamma_rad)
            pot_total += cos_theta
    return -pot_total  # Negatiu per a minimitzar amb scipy

# Cerca dels valors òptims
resultat = minimize(
    lambda x: potencia_total(x[0], x[1]),
    x0=[30, 180],  # Valors inicials (\beta=30^\circ, \gamma=180^\circ)
    bounds=[(0, 90), (0, 360)],  # Rangs per \beta i \gamma
    method='L-BFGS-B')

beta_opt, gamma_opt = resultat.x

print(f"Valor òptim de β (inclinació): {beta_opt:.2f} graus")
print(f"Valor òptim de γ (orientació): {gamma_opt:.2f} graus")
