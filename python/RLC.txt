"""
Calculating Power vs Load Resistance for parallel RLC circuit
 with power supply in series with Inductance.
Non-perfect Inductance and Capacitance modelled by additional
 resistances rl and rc.
Vpp maintains constant.


  ┌──────┬──────┐
  ~V     │      │       ~V: Voltage supply
  │      C      │       RL, rc, rl: Resistances
  L      │      RL      L: Inductance
  │      rc     │       C: Capacitance
  rl     │      │
  └──────┴──────┘
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# Constants
C = 22*10**-9  # Capacitance
L = 1.13*10**-3  # Inductance
V = 1  # Voltage of supply

# varying load Resistance
RL = np.linspace(0.001, 10**5, 10**6)

# Initial resistances and frequency
init_rl = 0
init_rc = 0
init_w = 1/(L*C)**(1/2)  # At resonance
print(f"'Resonant' frequency: {init_w/(2*np.pi)}")


def power_(rl, rc, C, L, w, V, RL):
    " Calculate Power dropped across load resistor RL "
    ZI = w*L*1j + rl
    ZC = 1/(w*C*1j) + rc
    Z = ZI + (1/ZC + 1/RL)**-1
    VL = V*(1-ZI/Z)

    # Take Re(VL)**2 or Re(VL**2) or |VL**2| or |VL|**2??
    # P = abs(VL)**2/RL
    P = abs(VL**2)/RL
    # P = VL.real**2/RL
    # P = (VL**2).real/RL # Doesnt work
    return P


# --- Matplotlib plotting code --- #

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = plt.plot(RL, power_(init_rl, init_rc, C, L, init_w, V, RL), lw=2)
ax.set_xlabel('RL [Ohm]')
ax.set_ylabel('Power [W]')

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
    valmin=-50,
    valmax=50,
    valinit=0,
)


# The function to be called anytime a slider's value changes
def update(val):
    w_ = init_w+2*np.pi*wd_slider.val
    P = power_(rl_slider.val,rc_slider.val,C,L,w_,V,RL)
    max_y = 1.05*max(P)
    line.set_ydata(P)
    ax.set_ylim(0, max_y)
    fig.canvas.draw_idle()


# register the update function with each slider
rl_slider.on_changed(update)
rc_slider.on_changed(update)
wd_slider.on_changed(update)

plt.show()
