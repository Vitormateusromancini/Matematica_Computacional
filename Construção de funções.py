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

