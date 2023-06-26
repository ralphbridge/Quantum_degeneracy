import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

global lam,E0,sigma

eps0=8.85e-12
c=3e8 # group velocity
vp=0.6*c # phase velocity
lam=532e-9 # laser frequency
k=2*np.pi/lam
w=5*k*vp
alpha=5e28 # Chirp correction
I0=1e14 # laser intensity
sigma=5*lam # laser variance
E0=np.sqrt(2*I0/(c*eps0))

# Definitions for Electric fields and Intensity as functions of z and t
def E(t):
    z=0
    return E0*np.exp(-0.5*((-z-c*t)/sigma)**2)*np.cos(-k*z-(w+alpha*t)*t)

z=np.linspace(-30e-6,30e-6,100000)

t=np.linspace(-0.6e-13,0.6e-13,5000)

S_l=np.zeros(len(t))
S_q=np.zeros(len(t))
j=0

for tau in t:
	S_l[j]=np.trapz((E(t)+E(t+tau))**2,t)
	S_q[j]=np.trapz(((E(t)+E(t+tau))**2)**2,t)
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
