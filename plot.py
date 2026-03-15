import numpy as np
import matplotlib.pyplot as plt

num = np.loadtxt("T_profile.dat")
ana = np.loadtxt("T_analytical.dat")

plt.figure()
plt.plot(num[:, 0], num[:, 1], label="Численное")
plt.plot(ana[:, 0], ana[:, 1], '--', label="Аналитическое")
plt.xlabel("x [m]")
plt.ylabel("T [K]")
plt.title("Temperature profile T(x)")
plt.legend()
plt.grid(True)
plt.savefig("T_profile.png", dpi=150)
plt.show()
