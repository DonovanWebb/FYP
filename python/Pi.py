"""
Calculating Power vs Load Resistance for parallel "Pi" RLC circuit
 with power supply in series with Inductance.
Non-perfect Inductance and Capacitance modelled by additional
 resistances rl and rc.
Vpp maintains constant.


"Pi" equivalent circuit:
 ┌──r───┬─L12─┬─────┬──────┐
        │     │     │      │       ~V: Voltage supply
 V      L1    L2    C      │       R, rc, rl: Resistances
        │     │     │      R       L: Inductance
        │     │     rc     │       C: Capacitance
        rl    rl    │      │
 └──────┴─────┴─────┴──────┘
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# Constants
C = 100*10**-9  # Capacitance
L1 = 1.5*10**-6  # Inductance L11
L2 = 1.5*10**-6  # Inductance L22
L12 = 1.13*10**-3      # Inductance L12
V = 1.0  # Voltage of supply

# varying load Resistance
R = np.linspace(0.001, 10**5, 10**6)

# Initial resistances and frequency
init_rl = 0
init_rc = 0
init_w = 1/(L2*C)**(1/2)  # At resonance
print(f"'Resonant' frequency: {init_w/(2*np.pi)}")


def power_tee(rl, rc, C, L1, L2, L12, w, V, R):
    " Calculate Power dropped across load resistor RL "
    r = 1.0
    ZL1 = w*L1*1j + rl
    ZL2 = w*L2*1j + rl
    ZL12 = w*L12*1j
    ZC = 1/(w*C*1j) + rc

    zCRL2 = (1/ZL2 + 1/ZC + 1/R)**-1
    Ztot = r + (1/ZL1 + 1/(ZL12+zCRL2))**-1

    Itot = V/Ztot

    Ir = Itot
    Vr = Ir*r
    Pr = Ir*Ir.conjugate() * r.real

    VL1 = V - Vr
    IL1 = VL1/ZL1
    PL1 = IL1*IL1.conjugate() * ZL1.real

    IL12 = Itot-IL1
    VL12 = IL12*ZL12
    PL12 = IL12*IL12.conjugate() * ZL12.real

    VL2 = VL1-VL12
    IL2 = VL2/ZL2
    PL2 = IL2*IL2.conjugate() * ZL2.real

    VC = VL2
    IC = VL2/ZC
    PC = IC*IC.conjugate() * ZC.real

    VR = VC
    IR = VR/R
    PR = VR*VR.conjugate()/R

    # assert all([round(PR[i].real, 5) == round(abs(PR[i]), 5) for i in range(len(PR))]) == True  # P should be real

    Zreal = np.sqrt(Ztot*Ztot.conjugate()) # useful

    Ptot = Pr + PL1 + PL12 + PL2 + PC + PR

    N = PR/Ptot


    # Checking Kirchoff rules:

    # 1st law
    # assert all([round(Ir[i], 5) == round((IL1[i] + IL12[i]), 5) for i in range(len(IL12))]) == True  # yes
    # assert all([round(IL12[i], 5) == round((IL2[i] + IR[i] + IC[i]), 5) for i in range(len(IL12))]) == True  # yes
    # assert all([round(Itot[i], 5) == round((IL1[i] + IL2[i] +IR[i] + IC[i]), 5) for i in range(len(IL12))]) == True  # yes
    
    # 2nd law
    # assert all([round(V, 5) == round((Vr[i] + VL12[i] + VL2[i]), 5) for i in range(len(IL12))]) == True  # yes
    # assert all([round(VL1[i], 5) == round((VL12[i] + VL2[i]), 5) for i in range(len(IL12))]) == True  # yes

    return PR, N



# --- Matplotlib plotting code --- #

# Create the figure and the line that we will manipulate
fig, (ax1, ax2) = plt.subplots(1,2)
P, N = power_tee(init_rl, init_rc, C, L1, L2, L12, init_w, V, R)
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

# Make a vertically oriented slider to control the Capacitance Resistance
axrc = plt.axes([0.25, 0.07, 0.65, 0.03])
rc_slider = Slider(
    ax=axrc,
    label="rc [Ohm]",
    valmin=0,
    valmax=100,
    valinit=init_rc,
    # orientation="vertical"
)

# Make a horizontal slider to control the difference to resonant frequency
axwd = plt.axes([0.25, 0, 0.65, 0.03])
wd_slider = Slider(
    ax=axwd,
    label='f diff [Hz]',
    valmin=-50000,
    valmax=50000,
    valinit=0,
)


# The function to be called anytime a slider's value changes
def update(val):
    w_ = init_w+2*np.pi*wd_slider.val
    P, N = power_tee(rl_slider.val,rc_slider.val,C,L1, L2, L12, w_,V,R)
    max_P = 1.05*max(P)
    max_N = 1.05*max(N)
    lineP.set_ydata(P)
    lineN.set_ydata(N)
    ax1.set_ylim(0, max_P)
    ax2.set_ylim(0, max_N)
    fig.canvas.draw_idle()


# register the update function with each slider
rl_slider.on_changed(update)
rc_slider.on_changed(update)
wd_slider.on_changed(update)

plt.show()
