# O método de diferenças finitas: EDPs elípticas
O método de diferenças finitas é uma abordagem numérica que aproxima as derivadas de uma equação diferencial parcial através de diferenças entre os valores discretos das variáveis ao longo de uma malha espacial, transformando a equação contínua em um sistema de equações algébricas solucionáveis.

A fórmula de diferenças finitas para aproximar numericamente a derivada de segunda ordem da função $u(x, y)$ na direção $x$ é:

$$ \frac{\partial^2 u}{\partial x^2} \approx \frac{u(x + h_x, y) - 2u(x, y) + u(x - h_x, y)}{h_x^2} $$

onde $h_x$ é o tamanho do passo na direção $x$ na malha discreta. Essa fórmula estima a derivada de segunda ordem utilizando os valores discretos da função $u$ em três pontos adjacentes na direção $x$. Na direção $y$ a fórmula é análoga.

$$ \frac{\partial^2 u}{\partial y^2} \approx \frac{u(x, y + h_y) - 2u(x, y) + u(x - h, y+h_y)}{h_y^2} $$

Se $h_x=h_y=h$ e, considerando que em cada ponto $(x_i,y_j)$ da malha, com $i=0,1,2,...,m$ e $j=0,1,2,...,n$, a função no ponto é dada por $u_{i,j}$ podemos escrever uma aproximação para a equação de Poisson no nó $i,j$ como

$$ \frac{u_{i+1,j} - 2u_{i,j} + u_{i-1,j}}{h^2} + \frac{u_{i,j+1} - 2u_{i,j} + u_{i,j-1}}{h^2} =f(x_i,y_j)$$

ou

$$ u_{i+1,j} +u_{i,j+1} - 4u_{i,j} + u_{i-1,j}+u_{i,j-1}=f(x_i,y_j)h^2$$


Desejamos agora resolver numericamente um problema de condução de calor

$$ \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

Para a solução numérica por diferenças finitas, usa-se uma grade de pontos onde as derivadas parciais da equação são substituidas por aproximações usando o esquema de diferenças finitas.

```python
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,2,21)
y = np.linspace(0,2,21)
xi, yj = np.meshgrid(x,y)

plt.plot(xi,yj,'b,')
plt.gca().set_aspect('equal')
plt.grid()

ni = 20
nj = 20

T = np.zeros([ni,nj])
T_ant = T.copy()
T[-1, :] = 50   # em cima
T[:, -1] = 100  # direita
T[0,  :] = 75   # embaixo
T[:,  0] = 0    # esquerda

err = 1000
n = 0
err_plot = []
int_plot = []

#print(T)
while err>0.0001:
    n=n+1
    int_plot.append(n)
    for i in range(1,len(T)-1):
        for j in range(1,len(T)-1):
            T[i,j] = (T[i-1,j]+T[i+1,j]+T[i,j-1]+T[i,j+1])/4
    err = np.linalg.norm(T-T_ant)/np.linalg.norm(T)
    T_ant = T.copy()

#plot
plt.pcolor(T[1:-1,1:-1],cmap='jet')
plt.gca().set_aspect('equal')
plt.colorbar()
plt.show()
```
# Projeto de Matematica Computacional

Projeto composto pelos discentes:

- Álvaro Augusto Lago Silva 
- Vitor Mateus Romancini do Amaral

Aulas ministradas e auxiliadas pelo docente:

- Prof: Tiago Martinuzzi Buriol

Neste projeto da disciplina vamos resolver numericamente equações diferenciais parciais (EDPs) envolve a discretização do domínio espacial e a implementação de métodos numéricos para aproximar as soluções. 
O trabalho escolhido foi as Deflexões em uma placa que dada o enunciado da questão pedia:  Desenvolva um programa computacional para determinar a deflexão de uma placa quadrada sujeita a uma carga de superfície constante. Teste seu programa para uma placa com bordas de 2 m de comprimento, q = 33,6 kN/m², σ = 0.3, Δz =  10−2  m e E =  2×1011  Pa. Use Δx = Δy = 0,5 m para seu teste de execução.

``` python
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
                z[i, j] = 0.25 * (z[i - 1, j] + z[i + 1, j]) + 0.25 * (z[i, j - 1] + z[i, j + 1]) - 0.25 * delta_x**2 * u[i, j]
                # z[i, j] = 0.25*((z[i - 1, j] + z[i + 1, j]) + (z[i, j - 1] + z[i, j + 1]) - 4*delta_x**2 * u[i, j])
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
```

Além da questão que pedia os mesmos cálculos da questão anterior, mas use 𝛥x = 𝛥y = 0.4 m
``` python
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
                z[i, j] = 0.25 * (z[i - 1, j] + z[i + 1, j]) + 0.25 * (z[i, j - 1] + z[i, j + 1]) - 0.25 * delta_x**2 * u[i, j]
                # z[i, j] = 0.25*((z[i - 1, j] + z[i + 1, j]) + (z[i, j - 1] + z[i, j + 1]) - 4*delta_x**2 * u[i, j])
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
    delta_x = 0.4  # Passo de discretização (m)

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
```
