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
