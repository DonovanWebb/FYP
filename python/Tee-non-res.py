"""
Calculating Power vs Load Resistance for parallel "Tee" RLC circuit
 with power supply in series with Inductance.
Non-perfect Inductance and Capacitance modelled by additional
 resistances rl and rc.
Vpp maintains constant.


"Tee" equivalent circuit:
 ┌────rl─L1──┬──L2─rl────┐
             │           │       ~V: Voltage supply
 V           M           │       R, rc, rl: Resistances
             │           R       L: Inductance
 r           │           │       C: Capacitance
 └───────────┴───────────┘

Consider the low freq case! hyper sensitive to any resistance in M as now r(M) > w*L(M)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# Constants
L1 = 10*10**-3  # Inductance L11 - M
L2 = 1*10**-3  # Inductance L22 - M
M = 10      # Inductance M
V = 20.0  # Voltage of supply

# varying load Resistance
R = np.linspace(0.001, 2*10**6, 10**3)

# Initial resistances and frequency
init_rl = 9
init_f = 2452


def calc_pow(I, V):  # I and V are PP values
    P = I.real*V.real #*np.cos(np.arctan(I.imag/I.real)-np.arctan(V.imag/V.real))/8
    return P


def power_tee(rl, L1, L2, M_, w, V, R):
    " Calculate Power dropped across load resistor RL "
    M = M_ * 10**-6

    ZL1 = w*L1*1j + rl
    ZL2 = w*L2*1j + rl

    ZM = w*M*1j #+ w*rl**3*10**-8 

    Ztot = ZL1 + (1/ZM + 1/ZL2)**-1

    Itot = V/Ztot

    IL1 = Itot
    VL1 = IL1*ZL1
    PL1 = calc_pow(IL1, VL1)

    VM = V - VL1
    IM = VM/ZM
    PM = calc_pow(IM, VM)

    IL2 = IL1 - IM
    VL2 = IL2 * ZL2
    PL2 = calc_pow(IL2, VL2)

    # IR = IL2
    # VR = IR*R
    # PR = calc_pow(IR, VR)

    # P = VR*VR.conjugate()/R
    # assert all([round(P[i].real, 5) == round(abs(P[i]), 5) for i in range(len(P))]) == True  # P should be real
    #P = P.real

    Ptot = PL1 + PM + PL2

    #if Ptot.all() != 0:
    if Ptot != 0:
        N = PL1/Ptot
    else:
        N = [1]



    # Checking Kirchoff rules:

    # 1st law
    # assert all([round(IL1[i], 5) == round((IM[i] + IL2[i]), 5) for i in range(len(IL1))]) == True  # yes
    # assert all([round(Itot[i], 5) == round((IR[i] + IM[i]), 5) for i in range(len(IL1))]) == True  # yes
    
    # 2nd law
    # assert all([round(V, 5) == round((VL1[i] + VM[i]), 5) for i in range(len(IL1))]) == True  # yes
    # assert all([round(VM[i], 5) == round((VL2[i] + VR[i]), 5) for i in range(len(IL1))]) == True  # yes

    return VL2.real, N, PL1


def find_freqs(V, R, rl):
    Ns_low = []
    freqs = range(1, 30000, 100)

    L2_ = 1.0*10**-3
    M_ = 1000     # Inductance M uH
    for f in freqs:
        w_ = 2*np.pi*f
        V_, N_, PL1_ = power_tee(rl, L1, L2_, M_, w_, V, R)
        max_V = max(V_)
        Ns_low.append(max_V)

    Ns_high = []
    L2_ = 1.2*10**-3
    M_ = 5000      # Inductance M uH
    for f in freqs:
        w_ = 2*np.pi*f
        V_, N_, PL1_ = power_tee(rl, L1, L2_, M_, w_, V, R)
        max_V = max(V_)
        Ns_high.append(max_V)

    N_ratio = [Ns_high[i]/Ns_low[i] for i in range(len(Ns_high))]

    return freqs, N_ratio, Ns_low, Ns_high

# --- Matplotlib plotting code --- #

# Create the figure and the line that we will manipulate
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3)
P, N, PL1 = power_tee(init_rl, L1, L2, M, init_f*2*np.pi, V, R)
lineP, = ax1.plot(R, P, lw=2)
lineN, = ax2.plot(R, N, lw=2)
ax1.set_xlabel('RL [Ohm]')
ax1.set_ylabel('Power [W]')
ax2.set_xlabel('RL [Ohm]')
ax2.set_ylabel('Efficiency')

# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control the Inductance Resistace
axrl = plt.axes([0.25, 0.15, 0.65, 0.03])
rl_slider = Slider(
    ax=axrl,
    label='rl [Ohm]',
    valmin=0,
    valmax=100,
    valinit=init_rl,
)

# Make a horizontal slider to control the difference to resonant frequency
axwd = plt.axes([0.25, 0, 0.65, 0.03])
wd_slider = Slider(
    ax=axwd,
    label='f [Hz]',
    valmin=1,
    valmax=30000,
    valinit=init_f,
)

freqs, N_ratio, Ns_low, Ns_high = find_freqs(V, R, init_rl)

lineNs, = ax4.plot(freqs, N_ratio, lw=2)
lineNlow, = ax5.plot(freqs, Ns_low, lw=2)
lineNhigh, = ax6.plot(freqs, Ns_high, lw=2)
ax4.set_ylim(-0.05*max(N_ratio), 1.05*max(N_ratio))
ax5.set_ylim(-0.05*max(Ns_low), 1.05*max(Ns_low))
ax6.set_ylim(-0.05*max(Ns_high), 1.05*max(Ns_high))

# register the update function with each slider

def update(val):
    w_ = 2*np.pi*wd_slider.val
    P, N, PL1 = power_tee(rl_slider.val, L1, L2, M, w_,V,R)
    max_P = 1.05*max(P)
    max_N = 1.05*max(N)

    lineP.set_ydata(P)
    lineN.set_ydata(N)
    ax1.set_ylim(-0.05*max_P, max_P)
    ax2.set_ylim(-0.05*max_N, max_N)


    freqs, N_ratio, Ns_low, Ns_high = find_freqs(V, R, rl_slider.val)

    lineNs.set_ydata(N_ratio)
    lineNlow.set_ydata(Ns_low)
    lineNhigh.set_ydata(Ns_high)
    ax4.set_ylim(-0.05*max(N_ratio), 1.05*max(N_ratio))
    ax5.set_ylim(-0.05*max(Ns_low), 1.05*max(Ns_low))
    ax6.set_ylim(-0.05*max(Ns_high), 1.05*max(Ns_high))
    fig.canvas.draw_idle()

rl_slider.on_changed(update)
wd_slider.on_changed(update)

plt.show()
