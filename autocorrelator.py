import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

global lam,E0,sigma

eps0=8.85e-12
c=3e8 # group velocity
vp=c # phase velocity
lam=800e-9 # laser frequency
k=2*np.pi/lam
w=3*k*vp
#alpha=1e29 # Chirp correction
I0=1e14 # laser intensity
#sigma=5*lam # laser variance
FWHM=8e-15
sigma=c*FWHM/(2*np.sqrt(2*np.log(2)))
E0=np.sqrt(2*I0/(c*eps0))
tmax=0.5e-13

# Definitions for Electric fields and Intensity as functions of z and t
def E(t,alpha):
    Ef=np.zeros(len(t))
    z=0
    Ef=E0*np.exp(-((-z-c*t)/(2*sigma))**2)*np.cos(-k*z-(w+alpha*t)*t)
    return Ef

def S_l(t,alpha):
    S=np.zeros(len(t))
    j=0
    for tau in t:
	    S[j]=np.trapz((E(t,alpha)+E(t+tau,alpha))**2,t)
	    j+=1
    return S

def S_q(t,alpha):
    S=np.zeros(len(t))
    j=0
    for tau in t:
	    S[j]=np.trapz(((E(t,alpha)+E(t+tau,alpha))**2)**2,t)
	    j+=1
    return S

t=np.linspace(-tmax,tmax,2000)

##################################################################################################################
############################################################################################### Slider

init_alpha=1.3e29

# Create the figure and the line that we will manipulate
fig,(ax1,ax2,ax3)=plt.subplots(3,1,tight_layout=True)

line1,=ax1.plot(t,E(t,init_alpha),lw=1)
ax1.set_xlabel(r'Time $t\ s$', fontsize=16)
ax1.set_ylabel('$E\ V/m$', fontsize=16)
ax1.set_title(r'Electric field', fontsize=16, color='r')
ax1.set_xlim([-tmax,tmax])

line2,=ax2.plot(t,S_l(t,init_alpha),lw=1)
ax2.set_xlabel(r'Time difference $\tau\ s$', fontsize=16)
ax2.set_ylabel('$S_{linear}\ J/m^2$', fontsize=16)
ax2.set_title(r'Linear autocorrelator', fontsize=16, color='r')
ax2.set_xlim([-0.8*tmax,0.8*tmax])

line3,=ax3.plot(t,S_q(t,init_alpha),lw=1)
ax3.set_xlabel(r'Time difference $\tau\ s$', fontsize=16)
ax3.set_ylabel('$S_{quadratic}\ J/m^2$', fontsize=16)
ax3.set_title(r'Quadratic autocorrelator', fontsize=16, color='r')
ax3.set_xlim([-0.8*tmax,0.8*tmax])

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.4, bottom=0.4)

# Make a horizontal slider to control the frequency.
axt=fig.add_axes([0.25, 0, 0.65, 0.03])
alpha_slider=Slider(
    ax=axt,
    label='Chirp',
    valmin=init_alpha,
    valmax=1e30,
    valinit=init_alpha,
)

#fig.savefig("chirp.pdf",bbox_inches='tight')

# The function to be called anytime a slider's value changes
def update(val):
    line1.set_ydata(E(t,alpha_slider.val))
    line2.set_ydata(S_l(t,alpha_slider.val))
    line3.set_ydata(S_q(t,alpha_slider.val))
    fig.canvas.draw_idle()

# register the update function with each slider
alpha_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values
resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(event):
    time_slider.reset()
button.on_clicked(reset)

plt.show()

############################################################################################### End of slider
##################################################################################################################

# Estimate dispersion for our laser: a) due to air, b) due to optics elements
# Start at a transform limit pulse (assume the pulse is good from the laser)
# Ask for the FROG at Kees' lab
# July 18
# Read Keese's DAMOP2005 tutorial to understand their parametrization, and look for the real parameters to give physical meaning to the chirp code I have implemented
# Added ssh key
