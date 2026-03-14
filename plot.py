import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("T_profile.dat")
x = data[:, 0]
T = data[:, 1]

plt.figure()
plt.plot(x, T)
plt.xlabel("x [m]")
plt.ylabel("T [K]")
plt.title("Temperature profile T(x)")
plt.grid(True)
plt.savefig("T_profile.png", dpi=150)
plt.show()
