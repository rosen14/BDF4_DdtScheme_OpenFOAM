import numpy as np
import matplotlib.pyplot as plt

# Parámetros
lamb = 1000.0
dt = 0.001/16/2
t_final = 0.01

# Tiempo
N = int(t_final / dt) + 3
t = np.linspace(-3*dt, t_final, N + 1)

# Solución exacta
phi_exact = np.cos(lamb * t)

# Inicialización
phi = np.zeros(N + 1)

# Inicialización con solución exacta (recomendado para evitar error de arranque)
phi[3] = np.cos(lamb * t[3])
phi[2] = np.cos(lamb * t[2])
phi[1] = np.cos(lamb * t[1])
phi[0] = np.cos(lamb * t[0])

print(phi[3])
print(phi[2])
print(phi[1])
print(phi[0])
# BDF4
for n in range(4, N + 1):
    rhs = 4*phi[n-1] - 3*phi[n-2] + (4/3)*phi[n-3] - (1/4)*phi[n-4]
    rhs -= dt * lamb * np.sin(lamb * t[n])
    phi[n] = (12 / 25) * rhs

# Error
error = np.abs(phi - phi_exact)

# Gráfico
plt.plot(t, phi, label="BDF4")
plt.plot(t, phi_exact, '--', label="Exacta")
plt.xlabel("t")
plt.ylabel("phi(t)")
plt.legend()
plt.title("Comparación BDF4 vs solución exacta")
plt.grid(True)
plt.show()

# Norma del error
print("Error L2:", np.linalg.norm(error, ord=2) * np.sqrt(dt))

'''
0.00100000
0.00050000
0.00025000
0.00012500
0.00006250
0.00003125
'''
