import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

eps0=8.85e-12
c=3e8
lamlaser=756e-9 # laser central wavelength
k=2*np.pi/lamlaser
w=2*np.pi*c/lamlaser
FWHM=8e-15
# T0=FWHM/(2*np.sqrt(2*np.log(2)))
T0=FWHM
E00=np.sqrt(2*1e14/(c*eps0))

gdd=50e-30 # Total GDD [(Air+mirrors+chirped mirrors+FSAC)+glass] in fs^2
# gdd=0
# alpha=1e28 # complex part of the beam parameter Gamma
alpha=w*gdd/T0**3

def E0shift(E0,j):
    Ef=np.zeros(len(E0))
    Eftemp=np.append(np.append(E0,E0),E0)
    for i in range(len(E0)):
        idx=int((len(E0)-1)/2)+i+j
        Ef[i]=Eftemp[idx]
    del Eftemp
    return Ef

# def S_l(t,alpha,E0):
#     S=np.zeros(len(t))
#     S=2*np.trapz((E(t,alpha,E0))**2,t)
#     S+=2*np.convolve(E(t,alpha,E0),E(t,alpha,E0),mode='same')
#     return S

def S_q(t,alpha,E0):
    S=np.zeros(len(t))
    j=0
    for tau in t:            
        S[j]=np.trapz((E0+E0shift(E0,j))**4,t)
        j+=1
    return S

n=10000
t=np.linspace(-1e-13,1e-13,n)
T1=np.sqrt(1+(gdd/T0**2)**2)*T0
Efield=(E00/(np.sqrt(2*np.pi)*T1))*np.exp(-t**2/(2*T1**2))*np.cos(w*t+alpha*t**2)
#Slin=S_l(t,alpha,E0)
Squad=S_q(t,alpha,Efield)

# Create the figure and the line that we will manipulate
fig,(ax1,ax3)=plt.subplots(2,1,tight_layout=True)

line1,=ax1.plot(t*1e15,Efield,lw=1)
ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax1.set_ylabel('$E\ V/m$', fontsize=16)
ax1.set_title(r'Electric field', fontsize=16, color='r')
ax1.set_xlim([-100,100])

# line2,=ax2.plot(t*1e15,Slin,lw=1)
# ax2.set_xlabel(r'Time difference $\tau\ fs$', fontsize=16)
# ax2.set_ylabel('$S_{linear}\ W/m^2$', fontsize=16)
# ax2.set_title(r'Linear detector', fontsize=16, color='r')
# ax2.set_xlim([-300,300])

line3,=ax3.plot(t*1e15 ,Squad,lw=1)
ax3.set_xlabel(r'Time difference $\tau\ fs$', fontsize=16)
ax3.set_ylabel('$S_{quadratic}\ W/m^2$', fontsize=16)
ax3.set_title(r'Quadratic detector', fontsize=16, color='r')
ax3.set_xlim([-100,100])

plt.show()