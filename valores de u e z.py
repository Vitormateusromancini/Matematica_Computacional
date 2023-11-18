import numpy as np

def deflection_plate_solver_u(size, q, sigma, delta_z, E, delta_x):
    # Parâmetros do sistema
    L = size  # Comprimento da placa
    W = size  # Largura da placa
    D = (E * delta_z**3) / (12 * (1 - sigma**2))  # Coeficiente de flexão

    # Número de nós em cada direção
    nx = int(L / delta_x) + 1
    ny = int(W / delta_x) + 1

    # Inicialização da matriz de deflexões
    u = np.zeros((nx, ny))

    # Aplicação das condições de contorno para u nas bordas
    u[:, 0] = 0  # Borda inferior fixa
    u[:, -1] = 0  # Borda superior fixa
    u[0, :] = 0  # Borda esquerda fixa
    u[-1, :] = 0  # Borda direita fixa

    # Resolução do sistema usando o método de diferenças finitas
    for _ in range(1000):  # Número de iterações
        u_old = np.copy(u)

        for i in range(1, nx-1):
            for j in range(1, ny-1):
                u[i, j] = 0.25 * (u[i + 1, j] + u[i - 1, j] + u[i, j + 1] + u[i, j - 1] - delta_x**2 * q / D)

        # Impor condição de contorno u = 0 nas bordas
        u[:, 0] = 0
        u[:, -1] = 0
        u[0, :] = 0
        u[-1, :] = 0

        # Verificar convergência
        if np.linalg.norm(u - u_old) < 1e-6:
            break

    return u

def deflection_plate_solver_z(u, delta_x):
    # Parâmetros do sistema
    nx, ny = u.shape

    # Inicialização da matriz de deflexões z
    z = np.zeros((nx, ny))

    # Aplicação das condições de contorno para z nas bordas
    z[:, 0] = 0  # Borda inferior fixa
    z[:, -1] = 0  # Borda superior fixa
    z[0, :] = 0  # Borda esquerda fixa
    z[-1, :] = 0  # Borda direita fixa

    # Resolução do sistema usando o método de diferenças finitas
    for _ in range(1000):  # Número de iterações
        z_old = np.copy(z)

        for i in range(1, nx-1):
            for j in range(1, ny-1):
                z[i, j] = 0.063 * (z[i - 1, j] + z[i + 1, j]) + 0.086 * (z[i, j - 1] + z[i, j + 1]) - 0.25 * delta_x**2 * u[i, j]

        # Impor condição de contorno z = 0 nas bordas
        z[:, 0] = 0
        z[:, -1] = 0
        z[0, :] = 0
        z[-1, :] = 0

        # Verificar convergência
        if np.linalg.norm(z - z_old) < 1e-6:
            break

    return z

def main():
    # Parâmetros do problema
    size = 2.0  # Tamanho da placa
    q = 33.6 * 10**3  # Carga de superfície constante (N/m^2)
    sigma = 0.3  # Coeficiente de Poisson
    delta_z = 1e-2  # Espessura da placa (m)
    E = 2e11  # Módulo de Elasticidade (Pa)
    delta_x = 0.5  # Passo de discretização (m)

    # Chamar a função do solver para u
    deflections_u = deflection_plate_solver_u(size, q, sigma, delta_z, E, delta_x)

    # Chamar a função do solver para z usando as deflexões u
    deflections_z = deflection_plate_solver_z(deflections_u, delta_x)

    # Exibir resultados
    print("Deflexões u na placa:")
    print(deflections_u)

    print("\nDeflexões z na placa:")
    print(deflections_z)

if __name__ == "__main__":
    main()
