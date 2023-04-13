import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

global lam,E0,sigma

lam=532e-9
eps0=8.85e-12
c=3e8
I0=1e14
E0=np.sqrt(2*I0/(c*eps0))
sigma=5*lam


# The parametrized function to be plotted
def f(z,t):
    return E0*np.exp(-0.5*((z)/sigma)**2)*np.cos(2*np.pi*(z)/lam)+E0*np.exp(-0.5*((-z-c*t)/sigma)**2)*np.cos(2*np.pi*(-z-2*c*t)/lam)
    
def I(z,t):
    return (E0*np.exp(-0.5*((z-c*t)/sigma)**2)*np.cos(2*np.pi*(z-2*c*t)/lam)*E0*np.exp(-0.5*((-z-c*t)/sigma)**2)*np.cos(2*np.pi*(-z-2*c*t)/lam))

z=np.linspace(-30e-6,30e-6,100000)

# Define initial parameters
init_t=-150e-15

# Create the figure and the line that we will manipulate
fig,(ax1,ax2)=plt.subplots(2)
line1,=ax1.plot(z,f(z,init_t),lw=2)
ax1.set_xlabel('Position [micro m]')
ax1.set(xlim=(min(z),max(z)),ylim=(-2*E0,2*E0))
line2,=ax2.plot(z,I(z,init_t),lw=2)
ax2.set_xlabel('Position [micro m]')
ax2.set(xlim=(min(z),max(z)),ylim=(-E0**2,E0**2))

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control the frequency.
axt=fig.add_axes([0.25, 0.1, 0.65, 0.03])
time_slider=Slider(
    ax=axt,
    label='Time [s]',
    valmin=-150e-15,
    valmax=150e-15,
    valinit=init_t,
)

# The function to be called anytime a slider's value changes
def update(val):
    line1.set_ydata(f(z,time_slider.val))
    line2.set_ydata(I(z,time_slider.val))
    fig.canvas.draw_idle()

# register the update function with each slider
time_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(event):
    time_slider.reset()
button.on_clicked(reset)

plt.show()

t=np.linspace(-0.6e-13,0.6e-13,10000)

E=E0*np.exp(-0.5*((c*t)/sigma)**2)*np.cos(2*np.pi*(2*c*t)/lam)

I0=0*t+2*np.trapz(E**2,t)
I0+=0*t+2*np.convolve(E,E,mode='same')
fig,ax3=plt.subplots()
line3,=ax3.plot(t,I0,lw=2)

plt.show()

I1=0*t+2*np.trapz(E**4,t)
I1+=0*t+8*np.convolve(E**3,E,mode='same')
I1+=0*t+6*np.convolve(E**2,E**2,mode='same')
fig,ax3=plt.subplots()
line3,=ax3.plot(t,I1,lw=2)

plt.show()
