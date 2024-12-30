import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
latitud = 41 + 39 / 60 + 55.115 / 3600  # Latitud de Bigues i Riells en graus
latitud_rad = np.radians(latitud)       # Latitud en radians
beta = np.radians(0)                    # Inclinació del panell en radians 
gamma = np.radians(180)                 # Panell orientat cap al sud
dies_any = 365                          # Dies en un any
hores_dia = 24                          # Hores en un dia

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
def cos_theta_i(delta, omega):
    term1 = np.sin(delta) * np.sin(latitud_rad) * np.cos(beta)
    term2 = -np.sin(delta) * np.cos(latitud_rad) * np.sin(beta) * np.cos(gamma)
    term3 = np.cos(delta) * np.cos(latitud_rad) * np.cos(beta) * np.cos(omega)
    term4 = np.cos(delta) * np.sin(latitud_rad) * np.sin(beta) * np.cos(gamma) * np.cos(omega)
    term5 = np.cos(delta) * np.sin(beta) * np.sin(gamma) * np.sin(omega)
    return term1 + term2 + term3 + term4 + term5

# Vector per emmagatzemar resultats
angles_incidents = []

# Simulació al llarg de l'any
for dia in range(1, dies_any + 1):
    delta = declinacio(dia)
    angles_dia = []
    for hora in range(hores_dia):
        omega = angle_horari(hora)
        elevacio = elevacio_solar(delta, omega)  # Calcula l'elevació solar
        if elevacio > 0:  # Només si el Sol està per sobre de l'horitzó
            cos_theta = cos_theta_i(delta, omega)
            theta_i = np.degrees(np.arccos(cos_theta)) if -1 <= cos_theta <= 1 else np.nan
        else:
            theta_i = np.nan  # Sol sota l'horitzó, angle d'incidència inexistent
        angles_dia.append(theta_i)
    angles_incidents.append(angles_dia)

angles_incidents = np.array(angles_incidents)

# Gràfic 3D de l'angle d'incidència
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')



# Eixos
dies = np.arange(1, dies_any + 1)
hores = np.arange(hores_dia)
X, Y = np.meshgrid(dies, hores)
Z = angles_incidents.T

# Superfície 3D
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
fig.colorbar(surf, ax=ax, label=r"$\theta_i$ (º)")
ax.set_title("Gràfic 3D de l'angle d'incidència corregit")
ax.set_xlabel("Dia de l'any")
ax.set_ylabel("Hora del dia")
ax.set_zlabel(r"$\theta_i$ (º)")
plt.show()

# Gràfic per a un dia concret
dia_concret = 1  # Podem anar posant el que volguem
angles_dia = angles_incidents[dia_concret - 1, :]

plt.figure(figsize=(10, 6))
plt.plot(range(hores_dia), angles_dia, marker='o')
plt.title(f"Angle d'incidència el dia {dia_concret}")
plt.xlabel("Hora del dia")
plt.ylabel(r"$\theta_i$ (º)")
plt.grid(True)
plt.show()

