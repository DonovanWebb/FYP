"""
Calculating Power vs Load Resistance for parallel "Tee" RLC circuit
 with power supply in series with Inductance.
Non-perfect Inductance and Capacitance modelled by additional
 resistances rl and rc.
Vpp maintains constant.


"Tee" equivalent circuit:
 ┌───────L1──┬──L2──┬──────┐
             │      │      │       ~V: Voltage supply
 V           M      C      │       R, rc, rl: Resistances
             │      │      R       L: Inductance
             │      rc     │       C: Capacitance
             rl     │      │
 └───────────┴──────┴──────┘
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# Constants
C = 27*10**-9  # Capacitance
L1 = 1.0*10**-3  # Inductance L11 - M
L2 = 1.0*10**-3  # Inductance L22 - M
M = 1      # Inductance M uH
V = 20  # Voltage of supply

# varying load Resistance
R = np.linspace(0.001, 10**3, 10**3)

# Initial resistances and frequency
init_rl = 5
init_rc = 20
init_w = 1/(L2*C)**(1/2)  # At resonance
print(f"'Resonant' frequency: {init_w/(2*np.pi)}")


def calc_pow(I, V):  # I and V are PP values
    #P = abs(I)*abs(V)*np.cos(np.arctan(I.imag/I.real)-np.arctan(V.imag/V.real))/8
    P = I.real*V.real #*np.cos(np.arctan(I.imag/I.real)-np.arctan(V.imag/V.real))/8
    return P



def power_tee(rl, rc, C, L1, L2, M_, w, V, R):
    " Calculate Power dropped across load resistor RL "
    M = M_ * 10**-6

    ZL1 = w*L1*1j + rc
    ZL2 = w*L2*1j + rl

    ZM = w*M*1j

    ZC = (w*C*1j)**-1

    zCR = (1/R + 1/ZC)**-1
    Ztot = ZL1 + (1/ZM + 1/(ZL2+zCR))**-1

    Itot = V/Ztot #abs(Ztot)

    IL1 = Itot
    VL1 = IL1*ZL1
    PL1 = calc_pow(IL1, VL1)

    VM = V - VL1
    IM = VM/ZM
    PM = calc_pow(IM, VM)

    IL2 = Itot - IM
    VL2 = IL2 * ZL2
    PL2 = calc_pow(IL2, VL2)

    VC = VM - VL2
    IC = VC/ZC
    PC = calc_pow(IC, VC)

    VR = VC
    IR = VR/R
    # PR = VR*VR.conjugate()/(8*R)
    PR = calc_pow(IR, VR)
    # print(PR[3125])

    # assert all([round(PR[i].real, 5) == round(abs(PR[i]), 5) for i in range(len(PR))]) == True  # P should be real

    Ptot = PL1 + PM + PL2 + PC + PR

    N = PR/PL1  # should strictly be PR/Ptot
    # plt.plot(R,PL1/Ptot)
    # plt.show()


    # Checking Kirchoff rules:

    # 1st law
    # assert all([round(IL1[i], 5) == round((IM[i] + IL2[i]), 5) for i in range(len(IL1))]) == True  # yes
    # assert all([round(Itot[i], 5) == round((IR[i] + IC[i] +IM[i]), 5) for i in range(len(IL1))]) == True  # yes
    # assert all([round(IL2[i], 5) == round((IR[i] + IC[i]), 5) for i in range(len(IL1))]) == True  # yes
    
    # 2nd law
    # assert all([round(V, 5) == round((VL1[i] + VM[i]), 5) for i in range(len(IL1))]) == True  # yes
    # assert all([round(V, 5) == round((VL1[i] + VL2[i] + VC[i]), 5) for i in range(len(IL1))]) == True  # yes
    # assert all([round(VM[i], 5) == round((VL2[i] + VC[i]), 5) for i in range(len(IL1))]) == True  # yes

    return PR, N, PL1



# --- Matplotlib plotting code --- #

# Create the figure and the line that we will manipulate
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3)
P, N, PL1 = power_tee(init_rl, init_rc, C, L1, L2, M, init_w, V, R)
lineP, = ax1.plot(R, P, lw=2)
lineN, = ax2.plot(R, N, lw=2)
linePL1, = ax3.plot(R, PL1, lw=2)
ax1.set_xlabel('RL [Ohm]')
ax1.set_ylabel('Power [W]')
ax2.set_xlabel('RL [Ohm]')
ax2.set_ylabel('Efficiency')
ax3.set_xlabel('RL [Ohm]')
ax3.set_ylabel('Power WT [W]')

# adjust the main plot to make room for the sliders
plt.subplots_adjust(bottom=0.25)

# Make a horizontal slider to control the Inductance Resistace
axrl = plt.axes([0.25, 0.15, 0.65, 0.03])
rl_slider = Slider(
    ax=axrl,
    label='rl [Ohm]',
    valmin=0,
    valmax=100,
    valinit=init_rl,
)

# Make a vertically oriented slider to control the Capacitance Resistance
axrc = plt.axes([0.25, 0.1, 0.65, 0.03])
rc_slider = Slider(
    ax=axrc,
    label="rc [Ohm]",
    valmin=0,
    valmax=20,
    valinit=init_rc,
    # orientation="vertical"
)

# Make a horizontal slider to control the difference to resonant frequency
axwd = plt.axes([0.25, 0.05, 0.65, 0.03])
wd_slider = Slider(
    ax=axwd,
    label='f diff [Hz]',
    valmin=-1*init_w/(2*np.pi),
    valmax=init_w/(2*np.pi),
    valinit=0,
)

axm = plt.axes([0.25, 0, 0.65, 0.03])
m_slider = Slider(
    ax=axm,
    label='Mutual inductance [uH]',
    valmin=0.01,
    valmax=10,
    valinit=M,
)

# The function to be called anytime a slider's value changes
def update(val):
    w_ = init_w+2*np.pi*wd_slider.val
    P, N, PL1 = power_tee(rl_slider.val,rc_slider.val,C,L1, L2, m_slider.val, w_,V,R)
    max_P = 1.05*max(P)
    max_N = 1.05*max(N)
    min_PL1 = 0.99*min(PL1)
    max_PL1 = 1.01*max(PL1)
    lineP.set_ydata(P)
    lineN.set_ydata(N)
    linePL1.set_ydata(PL1)
    ax1.set_ylim(0, max_P)
    ax2.set_ylim(0, max_N)
    ax3.set_ylim(min_PL1, max_PL1)
    fig.canvas.draw_idle()


# register the update function with each slider
rl_slider.on_changed(update)
rc_slider.on_changed(update)
wd_slider.on_changed(update)
m_slider.on_changed(update)


Ns_low = []
freqs = range(1, 30000,10)

L2_ = 1.0*10**-3
M_ = 5.0     # Inductance M uH
for f in freqs:
    C_ = 1/(L2_*(2*np.pi*f)**2)  # Capacitance
    w_ = 2*np.pi*f
    P, N, PL1 = power_tee(init_rl, init_rc, C_, L1, L2_, M_, w_, V, R)
    max_N = max(N)
    Ns_low.append(max_N)

Ns_high = []
L2_ = 1.2*10**-3
M_ = 10.0      # Inductance M uH
for f in freqs:
    C_ = 1/(L2_*(2*np.pi*f)**2)  # Capacitance
    w_ = 2*np.pi*f
    P, N, PL1 = power_tee(init_rl, init_rc, C_, L1, L2_, M_, w_, V, R)
    max_N = max(N)
    Ns_high.append(max_N)

N_ratio = [Ns_high[i]/Ns_low[i] for i in range(len(Ns_high))]
ax4.plot(freqs, N_ratio)
ax5.plot(freqs, Ns_low)
ax6.plot(freqs, Ns_high)
plt.show()
