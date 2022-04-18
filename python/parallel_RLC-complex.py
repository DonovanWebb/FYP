import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def sym():
    w, L, V, r_2, r_3, R = sp.symbols('w,L,V,r_2,r_3,R')
    Z, Z_b, P = sp.symbols('Z,Z_b,P')

    Z = w*L*1j + r_2 + (1/(r_3 - w*L*1j)+1/R)**-1
    Z_b = (Z*Z.conjugate())**(1/2)
    inner = (1-w*L/Z_b-r_2/Z_b)
    P = V**2*inner**2/R
    P = P.collect(R)
    R_max = sp.diff(P,R)
    return R_max


def find_P(w, R, L):
    I = 1
    M = L*10**-5
    V = w*I*M
    LT = 1*10**-3
    r_2 = 17
    r_3 = 1*10**-1


    Z = w*L*1j + r_2 + (1/(r_3 - w*L*1j)+1/R)**-1
    Z_b = (Z*Z.conjugate())**(1/2)
    P = V**2*(1-w*L/Z_b-r_2/Z_b)**2/R
    N = P/(w*LT*I**2)
    return N

R = np.linspace(50,1*10**4,10**4)


f_m = 30*10**3
fs = np.linspace(1000,f_m, 1000)

L_1 = 1*10**-3
N_1 = [max(find_P(w, R, L_1)) for w in 2*np.pi*fs]
N_1 = np.array(N_1)

L_2 = 1.2*10**-3
N_2 = [max(find_P(w, R, L_2)) for w in 2*np.pi*fs]
N_2 = np.array(N_2)

N = find_P(2*np.pi*1*10**3, R, L_1)


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
ax1.plot(R, N)
ax2.plot(fs, N_1)
ax3.plot(fs, N_2)
ax4.plot(fs, N_2/N_1)

plt.show()
