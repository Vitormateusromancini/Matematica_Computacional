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
