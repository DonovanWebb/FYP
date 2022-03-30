"""
Calculating Power vs Load Resistance for parallel "Tee" RLC circuit
 with power supply in series with Inductance.
Non-perfect Inductance and Capacitance modelled by additional
 resistances rl and rc.
Vpp maintains constant.


"Tee" equivalent circuit:
 ┌───────L1──┬──L2───────┐
             │           │       ~V: Voltage supply
 V           M           │       R, rc, rl: Resistances
             │           R       L: Inductance
             │           │       C: Capacitance
             rl          │
 └───────────┴───────────┘
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# Constants
L1 = 1.13*10**-3  # Inductance L11 - M
L2 = 1.13*10**-3  # Inductance L22 - M
M = 1.5*10**-6      # Inductance M
V = 1.0  # Voltage of supply

# varying load Resistance
R = np.linspace(0.001, 10**4, 10**4)

# Initial resistances and frequency
init_rl = 0
init_w = 0.01  # At resonance
print(f"frequency: {init_w/(2*np.pi)}")


def power_tee(rl, L1, L2, M, w, V, R):
    " Calculate Power dropped across load resistor RL "
    ZL1 = w*L1*1j + rl
    ZL2 = w*L2*1j + rl
    ZM = w*M*1j

    Ztot = ZL1 + (1/ZM + 1/(ZL2+R))**-1

    Itot = V/Ztot

    IL1 = Itot
    VL1 = IL1*ZL1

    VM = V - VL1
    IM = VM/ZM

    IL2 = IL1 - IM
    VL2 = IL2 * ZL2

    IR = IL2
    VR = IL2*R

    P = VR*VR.conjugate()/R
    # assert all([round(P[i].real, 5) == round(abs(P[i]), 5) for i in range(len(P))]) == True  # P should be real
    #P = P.real

    Zreal = np.sqrt(Ztot*Ztot.conjugate())
    Ptot = (V**2/Zreal)

    N = P/Ptot


    # Checking Kirchoff rules:

    # 1st law
    # assert all([round(IL1[i], 5) == round((IM[i] + IL2[i]), 5) for i in range(len(IL1))]) == True  # yes
    # assert all([round(Itot[i], 5) == round((IR[i] + IM[i]), 5) for i in range(len(IL1))]) == True  # yes
    
    # 2nd law
    # assert all([round(V, 5) == round((VL1[i] + VM[i]), 5) for i in range(len(IL1))]) == True  # yes
    # assert all([round(VM[i], 5) == round((VL2[i] + VR[i]), 5) for i in range(len(IL1))]) == True  # yes

    return P, N



# --- Matplotlib plotting code --- #

# Create the figure and the line that we will manipulate
fig, (ax1, ax2) = plt.subplots(1,2)
P, N = power_tee(init_rl, L1, L2, M, init_w, V, R)
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
    label='f diff [Hz]',
    valmin=0.01,
    valmax=30000,
    valinit=0,
)


# The function to be called anytime a slider's value changes
def update(val):
    w_ = init_w+2*np.pi*wd_slider.val
    P, N = power_tee(rl_slider.val, L1, L2, M, w_,V,R)
    max_P = 1.05*max(P)
    max_N = 1.05*max(N)
    lineP.set_ydata(P)
    lineN.set_ydata(N)
    ax1.set_ylim(0, max_P)
    ax2.set_ylim(0, max_N)
    fig.canvas.draw_idle()


# register the update function with each slider
rl_slider.on_changed(update)
wd_slider.on_changed(update)

plt.show()
