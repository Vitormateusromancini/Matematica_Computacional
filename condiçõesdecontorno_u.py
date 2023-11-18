import numpy as np
import matplotlib.pyplot as plt

def deflexao_placa_quadrada(q, sigma, Delta_z, E, L, delta_x, delta_y):
    # Parâmetros
    D = (E * Delta_z**3) / (12 * (1 - sigma**2))
    M, N = int(L / delta_x), int(L / delta_y)
    x = np.linspace(0, L, M)
    y = np.linspace(0, L, N)

    # Inicialização da matriz u
    u = np.zeros((M, N))

    # Condição de contorno u = 0 nas bordas
    u[:, 0] = 0
    u[:, -1] = 0
    u[0, :] = 0
    u[-1, :] = 0

    # Resolvendo a equação diferencial usando o método de diferenças finitas
    for _ in range(1000):  # Número de iterações
        for i in range(1, M - 1):
            for j in range(1, N - 1):
                u[i, j] = 0.25 * (u[i + 1, j] + u[i - 1, j] + u[i, j + 1] + u[i, j - 1]) - (q * delta_x**2) / D

    return u, x, y

# Parâmetros do problema
q = 33.6e3  # kN/m^2 para N/m^2
sigma = 0.3
Delta_z = 1e-2  # m
E = 2e11  # Pa
L = 2  # m
delta_x = 0.1  # m
delta_y = 0.1  # m

# Obtendo a deflexão
deflexao, x, y = deflexao_placa_quadrada(q, sigma, Delta_z, E, L, delta_x, delta_y)

# Visualização da deflexão
plt.imshow(deflexao, cmap='viridis', extent=[0, L, 0, L], origin='lower', interpolation='bilinear')
plt.colorbar(label='Deflexão (m)')
plt.title('Deflexão de uma placa quadrada')
plt.xlabel('Posição x (m)')
plt.ylabel('Posição y (m)')
plt.show()
