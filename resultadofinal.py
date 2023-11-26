import numpy as np
import matplotlib.pyplot as plt

def deflexao_placa_quadrada_u(q, sigma, Delta_z, E, L, delta_x, delta_y):
    # Parâmetros
    D = (E * Delta_z**3) / (12 * (1 - sigma**2))
    M, N = int(L / delta_x) + 1, int(L / delta_y) + 1
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
                u[i, j] = 0.25 * (u[i + 1, j] + u[i - 1, j] + u[i, j + 1] + u[i, j - 1] - delta_x**2 * q / D)

    return u, x, y


def deflexao_placa_quadrada_z(u, delta_x, x, y):
    # Parâmetros do sistema
    M, N = u.shape

    # Inicialização da matriz de deflexões z
    z = np.zeros((M, N))

    # Impor condição de contorno z = 0 nas bordas
    z[:, 0] = 0
    z[:, -1] = 0
    z[0, :] = 0
    z[-1, :] = 0

    # Resolução do sistema usando o método de diferenças finitas
    for _ in range(1000):  # Número de iterações
        for i in range(1, M-1):
            for j in range(1, N-1):
                z[i, j] = 0.25 * (z[i - 1, j] + z[i + 1, j]) + 0.25 * (z[i, j - 1] + z[i, j + 1]) - 0.25 * delta_x**2 * u[i, j]

    return z, x, y

    #     # Verificar convergência
    #     if np.linalg.norm(z - z_old) < 1e-6:
    #         break

    # return z

# Parâmetros do problema
q = float(input('Digite um valor para a carga distribuída (N/m²) :') ) # kN/m^2 para N/m^2
sigma = float(input('Digite um valor para \u03C3 :') )
Delta_z = float(input('Digite um valor para a espessura \u0394z (m) :') ) # m
E = float(input('Digite um valor para o Módulo de Elasticidade E (Pa) :') )  # Pa
L = float(input('Digite um valor para L (m) :') )  # m
delta_x = float(input('Digite um valor para \u0394x (m) :') ) # m
delta_y = float(input('Digite um valor para \u0394y (m) :') )  # m

# Obtendo a deflexão
deflexaoU, x, y = deflexao_placa_quadrada_u(q, sigma, Delta_z, E, L, delta_x, delta_y)
print('Deflexões u: ')
print(deflexaoU)

deflexaoZ, x, y = deflexao_placa_quadrada_z(deflexaoU, delta_x, x, y)
print('Deflexões z: ')
print(deflexaoZ)

# Visualização da deflexão
plt.imshow(deflexaoU, cmap='viridis', extent=[0, L, 0, L], origin='lower', interpolation='bilinear')
plt.colorbar(label='Deflexão (m)')
plt.title('Deflexão u de uma placa quadrada')
plt.xlabel('Posição x (m)')
plt.ylabel('Posição y (m)')
plt.show()

plt.imshow(deflexaoZ, cmap='viridis', extent=[0, L, 0, L], origin='lower', interpolation='bilinear')
plt.colorbar(label='Deflexão (m)')
plt.title('Deflexão z de uma placa quadrada')
plt.xlabel('Posição x (m)')
plt.ylabel('Posição y (m)')
plt.show()
