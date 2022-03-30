import numpy as np
import matplotlib.pyplot as plt


R_exp = np.array([0.5, 1, 10, 50, 100])*1000
T_exp = np.array([108.5, 100, 94, 93.8, 93])
T_ = T_exp-90

R = np.linspace(0.1, 5000, 100000)
y = (9.25-1.5)*1/R + 3

rc = -0.002
rl = 10
w = 2*np.pi*15000
C = 92.2 * 10**(-9)
L = 1.3 * 10**(-3)

# type 1
Z_partT = (R*rc*rl + rl)**2+w**2*(R*rl*C+R*rc*L + L)**2
Z_partB = R**2*(rl**2 + w**2*L**2)
Z = 1/np.sqrt(Z_partT/Z_partB)

# type 2
# Z_partT = (rl-rc)**2 + w**2*(R*C*rc+R*C*rl+L+C*rl*rc)**2
# Z_partB = R**2*((rl-rc)**2+w**2*(L+C*rc*rl)**2)
# Z = 1/np.sqrt(Z_partT/Z_partB)


VT = 1
# VL = VT*(Z-rl)/Z
# IL = VL/R
PL = Z/R
#plt.plot(x, y)
#plt.scatter(R_exp, T_)
plt.figure(1)
plt.scatter(R, Z)
plt.figure(2)
plt.scatter(R, PL)

plt.show()
