import numpy as np
import matplotlib.pyplot as plt

# --- График 1: сравнение прогибов ---
mech   = np.loadtxt("deflection_mech.dat")
thermo = np.loadtxt("deflection_thermo.dat")

fig1, ax1 = plt.subplots()
ax1.plot(mech[:, 0],   mech[:, 1],   label="Механический (без температуры)")
ax1.plot(thermo[:, 0], thermo[:, 1], label="Термоупругий")
ax1.set_xlabel("x / L")
ax1.set_ylabel("w [м]")
ax1.set_title("Прогиб мембраны")
ax1.legend()
ax1.grid(True)
fig1.savefig("deflection_compare.png", dpi=150)

# --- График 2: профиль температуры ---
prof = np.loadtxt("T_profile.dat")

fig2, ax2 = plt.subplots()
ax2.plot(prof[:, 0], prof[:, 1], label="T(z)")
ax2.set_xlabel("z [м]")
ax2.set_ylabel("T [К]")
ax2.set_title("Профиль температуры T(z) при t = 20 с")
ax2.legend()
ax2.grid(True)
fig2.savefig("T_profile.png", dpi=150)
