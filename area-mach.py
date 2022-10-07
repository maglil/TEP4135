import numpy as np
import matplotlib.pyplot as plt

M = np.linspace(0,5,100)

k = 1.4

area_ratio_fun = 1/M*( (1 + 1/2*(k-1)*M**2)/(1/2*(k+1)) )**(k+1)/(2*(k-1))

area_ratio = 2.5

eq = area_ratio_fun-area_ratio

plt.plot(M,eq)
plt.grid(True)
plt.show()

