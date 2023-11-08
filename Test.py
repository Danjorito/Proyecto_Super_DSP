import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 9999, 10000)
y = np.sin(2*np.pi*0.002*t)

E_est = 10;
E_mea = 0.02;
EST = np.zeros((10000, 1));
EST[0] = y[0]
KG = 0.9999;

for k in range(1, 9999):
    EST[k] = EST[k-1] + KG*(y[k-1] - EST[k-1]);
    E_est = (1.1 - KG)*E_est;
    KG = (E_est)/(E_est + E_mea);

plt.plot(EST)
plt.plot(y)
plt.show()