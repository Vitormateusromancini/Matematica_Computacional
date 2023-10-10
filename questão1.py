import numpy as np

# Parâmetros da placa
L = 2.0  # Comprimento das bordas da placa em metros
q = 33.6e3  # Carga de superfície em N/m^2 (33.6 kN/m^2)
sigma = 0.3  # Coeficiente de Poisson
delta_z = 1e-2  # Espessura da placa em metros
E = 2e11  # Módulo de elasticidade em Pa (2x10^11 Pa)
delta_x = 0.5  # Passo de discretização em x em metros
delta_y = 0.5  # Passo de discretização em y em metros

# Número de pontos discretos em x e y
n_x = int(L / delta_x)
n_y = int(L / delta_y)

# Inicialize uma matriz vazia para armazenar as deflexões
deflexoes = np.zeros((n_x, n_y))

# Calcule as deflexões
for i in range(n_x):
    for j in range(n_y):
        x = i * delta_x - L / 2
        y = j * delta_y - L / 2
        deflexoes[i, j] = (q * delta_x * delta_y / (16 * E * delta_z)) * (
            ((12 * (1 - sigma * 2) * x * 2) / L ** 4)
            + ((12 * (1 - sigma * 2) * y * 2) / L ** 4)
            - ((2 * (1 - sigma * 2) * x * y) / L * 4)
            - 1
        )

# Imprima as deflexões
print("Deflexões da placa:")
print(deflexoes)
