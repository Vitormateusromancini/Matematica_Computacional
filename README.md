# Construção de Funções 

O primeiro conteúdo discutido em Matemática Computacional é sobre construção de funções que são ensinadas na matemática básica, mas agora abordado em métodos númericos. Como no exemplo abaixo para a função de terceiro grau. 

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-2, 2, 400)

y = x**3 - 2*x

plt.figure(figsize=(8, 6))
plt.plot(x, y, label='$y = x^3 - 2x$', color='blue')
plt.title('Gráfico da função $y = x^3 - 2x$')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()
plt.show()

```
No exemplo vemos a mportância de inserir as bibliotecas numpy para os cálculos matemáticos e o matplotlib.pyplot que inseri a plotagem e manipulação de gráficos. 
Definir os dados se utiliza np.linspace() que tem como finalidade na definição do intervalo de valores de x, além de que os valores de y são definidos pela função logo abaixo demostrado pelo código. 

Plotar os gráficos usa-se plt.plot() que traça o gráfico da função, especificando o estilo da linha, cor e a legenda. Os outros plot tem como finalidade enriquecer mais com detalhes o gráfico. E para exibir o gráfico se usa o plt.show() como demonstrado na imagem abaixo. 


![Figure python](https://github.com/Vitormateusromancini/Matematica_Computacional/assets/77472862/f612efd4-9df0-4abd-9ba0-fefa9397b5f1)

# Interpolação Polinomial

Uma função interpolada é aquela que passa em todos os n pontos dustintos de um intervalo [a,b], portanto uma função f(x) tal que ela seja igual a: 

                                    f(xi) = yi 

Com i = 0,1,2....n. 

Existem algumas formas de resolver uma interpolação e uma delas é usar a fórmula interpolatória de Lagrange 
                                    
Seja $f(x)$ definida em $x_0, x_1,...,x_n$, (n+1) pontos distintos de um intervalo $[a,b]$ e $y_i=f(x_i)$, $i=0,1,...,n$. O polinômio interpolador de Lagrange é dado por 
    $$ P(x) = y_0 l_0(x) + y_1 l_1(x) + ...+y_n l_n(x)$$
   
$l_k(x)$, $k=0,1,...,n$ são polinômios de grau $n$ obtidos pela fórmula

$$ l_k(x)=\\frac{(x-x_0)(x-x_1)...(x-x_{k-1})(x-x_{k+1})...(x-x_n)}{(x_k-x_0)(x_k-x_1)...(x_k-x_{k-1})(x_k-x_{k+1})...(x_k-x_n)}$$

ou, de forma compacta

$$P(x)=\\sum_{k=0}^{n} y_k l_k(x) $$
 
com

$$ l_k(x) = \\prod_{j=0 e j\\neq k}^{n} \\frac{(x-x_j)}{(x_k-x_j)}$$

onde o simbolo acima significa o produto. Um exemplo utilizando python está abaixo: 

```python
import matplotlib.pyplot as plt
import numpy as np 
x = [0.0, 0.5, 1.0]
y = [1.3, 2.5, 0.9]

def P(x): return  -5.6*x**2 + 5.2*x + 1.3

xnew = np.linspace(-1, 2, num=51)
print (xnew)
plt.plot(x, y, '.', xnew, P(xnew),'-', )
plt.grid()
plt.show()

print ('P(0.5)=', P(0.5))

```

# Erros absolutos e relativos

Erro absoluto: É a diferença entre o valor medido (ou aproximado) e o valor real de uma grandeza. Ele indica a magnitude do erro em unidades da própria grandeza.
Erro relativo: É a razão entre o erro absoluto e o valor real, expressa em porcentagem. Ele indica a significância do erro em relação ao valor real.

Um exemplo considerando uma sequência:

```Python
import math
x = 1
x_ant = x
x_exa = (1+math.sqrt(5))/2
    
for k in range(20):
        x = 1+1/x
        dif = abs(x-x_ant)/abs(x)
        err = abs(x-x_exa)/abs(x_exa)
        print (x, x_exa, err, dif)
        x_ant = x
```

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

``` python

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

```
