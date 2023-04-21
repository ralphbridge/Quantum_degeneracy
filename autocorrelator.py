import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

plt.rcParams['text.usetex'] = True

global lam,E0,sigma

eps0=8.85e-12
c=3e8
lam=532e-9 # laser frequency
I0=1e14 # laser intensity
sigma=5*lam # laser variance
E0=np.sqrt(2*I0/(c*eps0))

# Definitions for Electric fields and Intensity as functions of z and t
def f(z,t):
    return E0*np.exp(-0.5*((z)/sigma)**2)*np.cos(2*np.pi*(z)/lam)+E0*np.exp(-0.5*((-z-c*t)/sigma)**2)*np.cos(2*np.pi*(-z-2*c*t)/lam)
    
def I(z,t):
    return (E0*np.exp(-0.5*((z)/sigma)**2)*np.cos(2*np.pi*(z)/lam)+E0*np.exp(-0.5*((-z-c*t)/sigma)**2)*np.cos(2*np.pi*(-z-2*c*t)/lam))**2

z=np.linspace(-30e-6,30e-6,100000)

# Define initial parameters
#init_t=-150e-15
init_t=-0.6e-13

# Figure properties for the slider
fig,(ax1,ax2)=plt.subplots(2)
line1,=ax1.plot(z,f(z,init_t),lw=2)
ax1.set_xlabel(r'$Position\ \mu m$')
ax1.set(xlim=(min(z),max(z)),ylim=(-2*E0,2*E0))
line2,=ax2.plot(z,I(z,init_t),lw=2)
ax2.set_xlabel(r'$Position\ \mu m$')
ax2.set(xlim=(min(z),max(z)),ylim=(-5*E0**2,5*E0**2))

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
I_pulses=np.linspace(-0.6e-13,0.6e-13,10000)
S_l=np.zeros(len(t))
S_q=np.zeros(len(t))
j=0

for tau in t:
	S_l[j]=np.trapz(I(z,tau),z)
	S_q[j]=np.trapz((I(z,tau))**2,z)
	j+=1
	
fig,(ax3,ax4)=plt.subplots(2,tight_layout=True)
ax3.plot(t,S_l,lw=2)

ax3.set_xlabel(r'Time difference $\tau\ s$', fontsize=16)
ax3.set_ylabel('$S_{linear}\ J/m^2$', fontsize=16)
ax3.set_title(r'Linear autocorrelator', fontsize=16, color='r')

ax4.plot(t,S_q,lw=2)

ax4.set_xlabel(r'Time difference $\tau\ s$', fontsize=16)
ax4.set_ylabel('$S_{quadratic}\ J/m^2$', fontsize=16)
ax4.set_title(r'Quadratic autocorrelator', fontsize=16, color='r')

plt.show()

#E=E0*np.exp(-0.5*((c*t)/sigma)**2)*np.cos(2*np.pi*(2*c*t)/lam) # Pulse electric field

#I0=2*np.convolve(E,E,mode='same') # Interference term (documentation, eq. (4.7) page 54)
#I0+=2*np.trapz(E**2,t) # Offset 

#fig,(ax3,ax4)=plt.subplots(2,tight_layout=True)
#ax3.plot(t,I0,lw=2)

#ax3.set_xlabel(r'Time difference $\tau\ s$', fontsize=16)
#ax3.set_ylabel('$S_{linear}\ J/m^2$', fontsize=16)
#ax3.set_title(r'Linear autocorrelator', fontsize=16, color='r')

#I1=8*np.convolve(E**3,E,mode='same') # Interference terms (documentation, eq. (4.8) page 54)
#I1+=6*np.convolve(E**2,E**2,mode='same')
#I1+=2*np.trapz(E**4,t) # Offset
#ax4.plot(t,I1,lw=2)

#ax4.set_xlabel(r'Time difference $\tau\ s$', fontsize=16)
#ax4.set_ylabel('$S_{quadratic}\ J/m^2$', fontsize=16)
#ax4.set_title(r'Quadratic autocorrelator', fontsize=16, color='r')

#plt.show()
