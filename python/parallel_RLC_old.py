import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button


init_rl = 0
init_rc = 0
C = 22*10**-9
L = 1.13*10**-3
init_w = 1/(L*C)**(1/2)
print(f"'Resonant' frequency: {init_w/(2*np.pi)}")
V = 1

RL = np.linspace(0.001,10**5,10**6)

def power_(rl,rc,C,L,w,V,RL):
    ZI = w*L*1j + rl
    ZC = 1/(w*C*1j) + rc
    Z = ZI + (1/ZC + 1/RL)**-1
    # Z = abs(Z)
    # Z2 = Z**2
    VL = V*(1-ZI/Z)
    # Take only real VL or take real VL**2 or abs or |VL|??

    # P = abs(VL)**2/RL
    P = abs(VL**2/RL)
    # P = VL.real**2/RL
    # P = (VL**2).real/RL # Doesnt work
    return P

def power(rl,rc,C,L,w,V,RL):

    # Real part of Z numerator
    A = rl - w**2*L*C**rl*(RL+rc)+RL+w**2*L*C*(RL+rc)+w**2*C**2*rl*(RL+rc)**2+w**2*C**2*RL*rc*(RL+rc)
    # Im part of Z numerator
    B = w*(w**2*L*C**2*rl*(RL+rc)**2-C*RL*(RL+rc)+L+C*RL*rc)
    # Denominator
    denom = 1+w**2*C**2*(RL+rc)**2

    #Z = (A**2 - B**2)**(1/2)/C
    Z2 = (A**2 - B**2)/denom**2
    Z = Z2**0.5

    P = V*RL/Z2
    Z_ = ((RL**2+(w*C*RL**2+w*L)**2)**0.5)/(1+w**2*C**2*RL**2)
    Z__ = ((RL**2+(w*L-RL**2*w*C+RL**2*w)**2)/(1+RL**2*w**2*C**2))**0.5
    return Z2


# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = plt.plot(RL, power_(init_rl,init_rc,C,L,init_w,V,RL), lw=2)
ax.set_xlabel('RL [Ohm]')

# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control the frequency.
axrl = plt.axes([0.25, 0.15, 0.65, 0.03])
rl_slider = Slider(
    ax=axrl,
    label='rl [Ohm]',
    valmin=0,
    valmax=100,
    valinit=init_rl,
)

axwd = plt.axes([0.25, 0, 0.65, 0.03])
wd_slider = Slider(
    ax=axwd,
    label='f diff [Hz]',
    valmin=-50,
    valmax=50,
    valinit=0,
)

# Make a vertically oriented slider to control the amplitude
#axrc = plt.axes([0.1, 0.25, 0.0225, 0.63])
axrc = plt.axes([0.25, 0.07, 0.65, 0.03])
rc_slider = Slider(
    ax=axrc,
    label="rc [Ohm]",
    valmin=0,
    valmax=100,
    valinit=init_rc,
    # orientation="vertical"
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





#plt.plot(RL,P)
#plt.plot(RL,A)
plt.show()
