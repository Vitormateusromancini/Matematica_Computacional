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
Desejamos agora resolver numericamente um problema de condução de calor

$$ \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

Para a solução numérica por diferenças finitas, usa-se uma grade de pontos onde as derivadas parciais da equação são substituidas por aproximações usando o esquema de diferenças finitas.

# Projeto de Matematica Computacional

Neste projeto da disciplina vamos resolver numericamente equações diferenciais parciais (EDPs) envolve a discretização do domínio espacial e a implementação de métodos numéricos para aproximar as soluções. 
O trabalho escolhido foi as Deflexões em uma placa 
