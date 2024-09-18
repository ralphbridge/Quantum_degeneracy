import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate as inter

eps0=8.85e-12
c=299792458
lambda0=756e-9 # laser central wavelength
w0=2*np.pi*c/lambda0

def InverseFourier(Fw,w,t):
    ft=np.zeros(np.size(w),dtype=np.complex128)
    for j in range(len(t)):
        for i in range(len(w)):
            ft[j]+=Fw[i]*np.exp(1j*w[i]*t[j])+np.conj(Fw[i])*np.exp(-1j*w[i]*t[j])
    return np.real(ft)

#################################### Getting measured spectra

S=pd.read_excel("spectrum.xlsx").to_numpy()
n=np.size(S,0)
lam=np.zeros(np.size(S,0))
Il=np.zeros(np.size(S,0))

for j in range(0,np.size(S,1)):
    for i in range (0,np.size(S,0)):
        #print(S[i][j])
        if j==0:
            lam[i]=(S[i][j])*1e-9
        elif j==1:
            Il[i]=S[i][j]

plt.plot(lam*1e9,Il)
plt.xlabel('wavelength nm')
plt.grid()
plt.show()

f=np.zeros(n)
spectrumf=np.zeros(n)

for i in range(n):
    f[n-i-1]=c/(lam[i])
    spectrumf[n-i-1]=(lam[i]**2)*Il[i]/c

plt.plot(f,spectrumf)
plt.xlabel('frequency')
plt.title('Original frequency spectrum')
plt.grid()
plt.show()

######################## Increasing time resolution (by increasing frequency range)

Ef0=np.sqrt(2*spectrumf/(c*eps0))
t0=np.linspace(-50e-15,50e-15,len(f))
Et0=InverseFourier(Ef0,2*np.pi*f,t0)

########################## Trimming spectrum low and high frequencies and interpolating it to get equally spaced frequencies

NN=799 # Length of actual spectrum (obtained from matlab)
f_trim=np.zeros(NN)
spectrumf_trim=np.zeros(NN)

for i in range(0,NN):
    f_trim[i]=f[i+266] # 266 is the position of the original spectrum where it starts being bigger than 0
    spectrumf_trim[i]=spectrumf[i+266] # New angular frequency spectrum (trimmed)

plt.plot(f_trim,spectrumf_trim)
plt.xlabel('frequency')
plt.title('Trimmed frequency spectrum')
plt.grid()
plt.show()

spectrumf_interp_function=inter.CubicSpline(f_trim,spectrumf_trim)

N=799
f_interp=np.linspace(min(f_trim),max(f_trim),N)
spectrumf_interp=np.zeros(N)
df=f_interp[1]-f_interp[0]

w=2*np.pi*f_interp

for i in range(N):
    if i==0:
        spectrumf_interp[i]=spectrumf_trim[i]
    elif i==N-1:
        spectrumf_interp[i]=spectrumf_trim[NN-1]
    else:
        spectrumf_interp[i]=spectrumf_interp_function(f_interp[i])
        if spectrumf_interp[i]<0:
            spectrumf_interp[i]=0

Ef0_interp=np.sqrt(2*spectrumf_interp/(c*eps0))*(1+0j) # Initial electric field spectrum (trimmed)

plt.plot(f_interp,spectrumf_interp)
plt.xlabel('frequency')
plt.title('Interpolated frequency spectrum')
plt.grid()
plt.show()

######################################

t0_interp=np.linspace(-50e-15,50e-15,len(f_interp))
Et0_interp=InverseFourier(Ef0_interp,2*np.pi*f_interp,t0_interp)

fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)
plt.subplot(2,1,1)
line1,=plt.plot(t0*1e15,Et0)
# ax1.set_xlim([min(t0)*1e15/10,max(t0)*1e15/10])
ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax1.set_ylabel(r'$E_0\ V/m$', fontsize=16)
ax1.grid(which='both')

plt.subplot(2,1,2)
line2,=ax2.plot((t0_interp)*1e15,Et0_interp)
ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax2.set_ylabel('$E_0^{interp}\ V/m$', fontsize=16)
# ax2.set_xlim([-50,50])
ax2.grid(which='both')
plt.show()