import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import sympy as sym
from scipy import interpolate as inter

eps0=8.85e-12
c=299792458
lambda0=756e-9 # laser central wavelength
w0=2*np.pi*c/lambda0
E00=np.sqrt(2*1e14/(c*eps0))

zd=2
zw=3e-3

x=sym.Symbol('x')

n_air=1+0.05792105/(238.0185-(1e-12)*(x/(2*np.pi*c))**2)+0.00167917/(57.362-(1e-12)*(x/(2*np.pi*c))**2)
n_bk7=(1+1.03961212*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-0.00600069867)+0.231792344*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-0.0200179144)+1.01046945*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-103.560653))**0.5

n0_air=n_air.subs(x,w0)
np0_air=sym.diff(n_air,x).subs(x,w0)
npp0_air=sym.diff(sym.diff(n_air,x),x).subs(x,w0)
nppp0_air=sym.diff(sym.diff(sym.diff(n_air,x),x)).subs(x,w0)

k0_air=float(n0_air*w0/c)
kp0_air=float((n0_air+w0*np0_air)/c)
kpp0_air=float((2*np0_air+w0*npp0_air)/c)
kppp0_air=float((3*npp0_air+w0*nppp0_air)/c)

n0_bk7=n_bk7.subs(x,w0)
np0_bk7=sym.diff(n_bk7,x).subs(x,w0)
npp0_bk7=sym.diff(sym.diff(n_bk7,x),x).subs(x,w0)
nppp0_bk7=sym.diff(sym.diff(sym.diff(n_bk7,x),x)).subs(x,w0)

k0_bk7=float(n0_bk7*w0/c)
kp0_bk7=float((n0_bk7+w0*np0_bk7)/c)
kpp0_bk7=float((2*np0_bk7+w0*npp0_bk7)/c)
kppp0_bk7=float((3*npp0_bk7+w0*nppp0_bk7)/c)

def E0shift(E0,j):
    Ef=np.zeros(len(E0))
    Eftemp=np.append(np.append(E0,E0),E0)
    for i in range(len(E0)):
        idx=int((len(E0)-1)/2)+i+j
        Ef[i]=Eftemp[idx]
    del Eftemp
    return Ef

def S_l(t,E0):
    S=np.zeros(len(t))
    j=0
    for tau in t:            
        S[j]=np.trapz((E0+E0shift(E0,j))**2,t)
        j+=1
    return S

def S_q(t,E0):
    S=np.zeros(len(t))
    j=0
    for tau in t:            
        S[j]=np.trapz((E0+E0shift(E0,j))**4,t)
        j+=1
    return S

def roundup(x):
    return int(math.ceil(x / 100.0)) * 100

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

Il=(Il-min(Il))/max(Il)

spectruml=Il # Initial measured spectrum

f=np.zeros(n)
spectrum=np.zeros(n)

for i in range(n):
    f[n-i-1]=c/(lam[i])
    spectrum[n-i-1]=(lam[i]**2)*spectruml[i]/c
    
df=(max(f)-min(f))/n

# f=np.append(f,np.linspace(max(f)+df,3*(max(f)+max(f)-min(f)),3*n))
# spectrum=np.append(spectrum,np.zeros(3*n))
# n=len(f)
# df=(max(f)-min(f))/n

Ef=np.sqrt(2*spectrum/(c*eps0))*(1+0j) # Initial electric field

# Initial phases for original measured spectrum (up to third order)
# i=0
# for w in 2*np.pi*f:
#     Ef[i]=Ef[i]*np.exp(1j*GDD_laser[i]*(w-w0)**2/math.factorial(2))
#     Ef[i]=Ef[i]*np.exp(1j*TOD_laser[i]*(w-w0)**3/math.factorial(3))
    
#     i+=1
    
t=np.arange(-n/2,n/2)/(n*df)
pulse=np.fft.ifftshift(np.fft.ifft(Ef))

plt.plot(t*1e15,pulse)
plt.xlim([-20,20])
plt.show()

########################## Getting GDD data from Thorlabs chirped mirrors data

CM=pd.read_excel("UMxx-15FS_data.xlsx").to_numpy()
n_cm=np.size(CM,0)
lam_cm=np.zeros(np.size(CM,0))
w_cm=np.zeros(np.size(CM,0))
GD_cm=np.zeros(np.size(CM,0))

for j in range(0,np.size(CM,1)):
    for i in range (0,np.size(CM,0)):
        if j==0:
            lam_cm[i]=(CM[i][j])*1e-9
            w_cm[i]=2*np.pi*c/((CM[n_cm-1-i][j])*1e-9)
        else:
            GD_cm[i]=(CM[i][j])*1e-15

GDD_cm_data=np.zeros(n_cm)

for i in range(0,n_cm):
    if i==0:
        m2=(GD_cm[i+1]-GD_cm[i])/(w_cm[i+1]-w_cm[i])
        GDD_cm_data[i]=m2
    elif i==n_cm-1:
        m1=(GD_cm[i]-GD_cm[i-1])/(w_cm[i]-w_cm[i-1])
        GDD_cm_data[i]=m1
    else:
        m1=(GD_cm[i]-GD_cm[i-1])/(w_cm[i]-w_cm[i-1])
        m2=(GD_cm[i+1]-GD_cm[i])/(w_cm[i+1]-w_cm[i])
        GDD_cm_data[i]=(m1+m2)/2

GDD_cm_interp=inter.interp1d(w_cm,GDD_cm_data)
GDD_cm=np.zeros(n)

# i=0
# for w in 2*np.pi*f:
#     GDD_cm[i]=GDD_cm_interp(w)
#     i+=1

plt.plot(w_cm,GDD_cm_data*1e30)
plt.xlim([2e15,3e15])
plt.ylim([-100,200])
plt.grid()
plt.show()

print(GDD_cm_interp(2*np.pi*c/800e-9)*1e30) # Check this line, it does not interpolate correctly

############################## Computing the phases due to each element in the layout

Em=np.zeros(np.size(Ef),dtype=np.complex_)

i=0
for w in 2*np.pi*f:
    # Dispersion phases introduced by air up to the third order
    Em[i]=Ef[i]*np.exp(1j*kp0_air*zd)*np.exp(1j*kp0_air*(w-w0)*zd)
    Em[i]=Em[i]*np.exp(1j*kpp0_air*(w-w0)**2*zd/math.factorial(2))
    Em[i]=Em[i]*np.exp(1j*kppp0_air*(w-w0)**3*zd/math.factorial(3))
    
    # Dispersion phases introduced by BK7 Fused Silica glass window up to the third order
    Em[i]=Em[i]*np.exp(1j*kp0_bk7*zd)*np.exp(1j*kp0_bk7*(w-w0)*zw)
    Em[i]=Em[i]*np.exp(1j*kpp0_bk7*(w-w0)**2*zw/math.factorial(2))
    Em[i]=Em[i]*np.exp(1j*kppp0_bk7*(w-w0)**3*zw/math.factorial(3))
    
    # Dispersion phases introduced by bounces off P01 Silver mirrors up to the third order
    Em[i]=Em[i]*np.exp(1j*GDD_p01m*(w-w0)**2/math.factorial(2))
    Em[i]=Em[i]*np.exp(1j*TOD_p01m*(w-w0)**3/math.factorial(3))
    
    # Dispersion phases introduced by bounces off Chirped mirrors up to the third order
    Em[i]=Em[i]*np.exp(1j*GDD_cm[i]*(w-w0)**2/math.factorial(2))
    Em[i]=Em[i]*np.exp(1j*TOD_cm[i]*(w-w0)**3/math.factorial(3))
    
    i+=1

spectrum_final=1

####################################

######## Creating equally spaced frequency domain and spectrum

# dfmax=0
# dfmin=df

# for i in range(n-1):
#     if f[i+1]-f[i]>dfmax:
#         dfmax=f[i+1]-f[i]
#     elif f[i+1]-f[i]<dfmin:
#         dfmin=f[i+1]-f[i]

# df=df

# ftemp=np.arange(min(f),max(f),df)

# n=len(ftemp)
# spectrumtemp=np.zeros(n)

# spectrum_interp=inter.interp1d(f,spectrum)
# # spectrum_interp=inter.UnivariateSpline(f,spectrum)

# for i in range(n):
#     spectrumtemp[i]=spectrum_interp(ftemp[i])

# f=ftemp
# spectrum=spectrumtemp

# del ftemp, spectrumtemp

# ########## Increasing time resolution (by increasing frequency range) by 3^N

# N=1

# for i in range(N):
#     f=np.append(np.append(np.arange(min(f)-n*df,min(f),df),f),np.arange(max(f)+df,max(f)+(n+1)*df,df))
#     spectrum=np.append(np.append(np.zeros(n),spectrum),np.zeros(n))
#     n=len(f)

###########

# pulse=np.fft.ifftshift(np.fft.ifft(spectrum))

# It=abs(pulse)

fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)
plt.subplot(2,1,1)
line1,=ax1.plot(lam*1e9,spectruml)
ax1.set_xlabel(r'Wavelength $\lambda\ nm$', fontsize=16)
ax1.set_ylabel(r'$I\ W/m^2$', fontsize=16)
ax1.set_xlim([min(lam)*1e9,max(lam)*1e9])

major_tick = np.arange(roundup(min(lam)*1e9), roundup(max(lam)*1e9),200)#[200, 400, 600, 800, 1000]
minor_tick = np.arange(roundup(min(lam)*1e9)+100, roundup(max(lam)*1e9),200)#[300, 500, 700, 900]
ax1.set_xticks(major_tick) # Grid
ax1.set_xticks(minor_tick, minor=True)
ax1.grid(which='both')

plt.subplot(2,1,2)
line2,=ax2.plot(t*1e15,It)
ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax2.set_ylabel('$I\ W/m^2$', fontsize=16)
ax2.set_xlim([-100,100])

major_tick = np.arange(-100,100,20)
minor_tick = np.arange(-100,100,10)
ax2.set_xticks(major_tick)
ax2.set_xticks(minor_tick, minor=True)
ax2.grid(which='both')

# fig.savefig("Frequency_Time.pdf",bbox_inches='tight')
plt.show()

############### Gaussian fit (3 Gaussian functions)

N=5000

tt=t[len(t)//2-60:len(t)//2+60]
Itt=It[len(t)//2-60:len(t)//2+60]

plt.plot(t*1e15,It)
# plt.xlim(-20,20)
plt.xlabel(r'Time $t\ fs$', fontsize=16)
plt.ylabel('Intensity $I\ W/m^2$', fontsize=16)
# major_tick = np.arange(-100,100,20)
# minor_tick = np.arange(-100,100,10)
# plt.xticks(major_tick)
# plt.xticks(minor_tick, minor=True)
plt.grid(which='both')
# plt.legend(["Gaussian fit","Inverse Fourier transform"],loc="upper right")
# plt.savefig("GaussianFitting.pdf",bbox_inches='tight')
plt.show()

###############

Efield=np.sqrt(2*It/(c*eps0))*np.cos(w*t+alpha*t**2)
Slinear=S_l(t,Efield)
Squad=S_q(t,Efield)

fig,(ax1,ax3)=plt.subplots(2,1,tight_layout=True)

line1,=ax1.plot(t*1e15,Efield,lw=1)
ax1.plot(t*1e15,np.sqrt(2*It/(c*eps0)),'.',alpha=0.01)
ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax1.set_ylabel(r'$E\ V/m$', fontsize=16)
ax1.set_xlim([-50,50])
ax1.grid()

# line2,=ax2.plot(t*1e15 ,Slinear,lw=1)
# ax2.set_xlabel(r'Time delay $\tau\ fs$', fontsize=16)
# ax2.set_ylabel(r'$S_{linear}\ W/m^2$', fontsize=16)
# ax2.set_xlim([-20,20])
# ax2.grid()

line3,=ax3.plot(t*1e15 ,Squad,lw=1)
ax3.set_xlabel(r'Time delay $\tau\ fs$', fontsize=16)
ax3.set_ylabel(r'$S_{quadratic}\ W/m^2$', fontsize=16)
ax3.set_xlim([-50,50])
ax3.grid()
# plt.savefig("FieldTrace_chirp.pdf",bbox_inches='tight')
plt.show()

# plt.plot(t_fit*1e15 ,Squad,lw=1)
# plt.xlabel(r'Time delay $\tau\ fs$', fontsize=16)
# plt.ylabel(r'$S_{quadratic}\ W/m^2$', fontsize=16)
# plt.xlim([-100,100])
# plt.grid()
# plt.show()

# Check difference between calculated and experimental total GDD
# Discuss Agrawal's expressions with Herman
# Get theoretical expression for Fourier transform (for simple Gaussian spectrum)
