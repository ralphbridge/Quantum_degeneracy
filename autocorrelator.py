import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import sympy as sym

eps0=8.85e-12
c=3e8
lambda0=756e-9 # laser central wavelength
w0=2*np.pi*c/lambda0
E00=np.sqrt(2*1e14/(c*eps0))

x=sym.Symbol('x')

n_air=1+0.05792105/(238.0185-(1e-12)*(x/(2*np.pi*c))**2)+0.00167917/(57.362-(1e-12)*(x/(2*np.pi*c))**2)
n_bk7=(1+1.03961212*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-0.00600069867)+0.231792344*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-0.0200179144)+1.01046945*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-103.560653))**0.5

n0_air=n_air.subs(x,w0)
np0_air=sym.diff(n_air,x).subs(x,w0)
npp0_air=sym.diff(sym.diff(n_air,x),x).subs(x,w0)
nppp0_air=sym.diff(sym.diff(sym.diff(n_air,x),x)).subs(x,w0)

k0_air=n0_air*w0/c
kp0_air=(n0_air+w0*np0_air)/c
kpp0_air=(2*np0_air+w0*npp0_air)/c
kppp0_air=(3*npp0_air+w0*nppp0_air)/c

n0_bk7=n_bk7.subs(x,w0)
np0_bk7=sym.diff(n_bk7,x).subs(x,w0)
npp0_bk7=sym.diff(sym.diff(n_bk7,x),x).subs(x,w0)
nppp0_bk7=sym.diff(sym.diff(sym.diff(n_bk7,x),x)).subs(x,w0)

k0_bk7=n0_bk7*w0/c
kp0_bk7=(n0_bk7+w0*np0_bk7)/c
kpp0_bk7=(2*np0_bk7+w0*npp0_bk7)/c
kppp0_bk7=(3*npp0_bk7+w0*nppp0_bk7)/c

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
###############

S=pd.read_excel("spectrum.xlsx").to_numpy()
n=np.size(S,0)
lam=np.zeros(np.size(S,0))
Il=np.zeros(np.size(S,0))
Iluv=np.zeros(np.size(S,0))
Ilir=np.zeros(np.size(S,0))

for j in range(0,np.size(S,1)):
    for i in range (0,np.size(S,0)):
        #print(S[i][j])
        if j==0:
            lam[i]=(S[i][j])*1e-9
        elif j==1:
            Il[i]=S[i][j]
        elif j==2:
            Iluv[i]=S[i][j]
        else:
            Ilir[i]=S[i][j]

Il=(Il-min(Il))/max(Il)
Iluv=(Iluv-min(Iluv))/max(Iluv)
Ilir=(Ilir-min(Ilir))/max(Ilir)

spectruml=Il

f=np.zeros(n)
spectrum=np.zeros(n)

for i in range(n):
    f[n-i-1]=c/(lam[i])
    spectrum[n-i-1]=(lam[i]**2)*spectruml[i]/c

df=(max(f)-min(f))/n

####################################

Ef=np.sqrt(2*spectrum/(c*eps0))

i=0

Em=np.zeros(size(Ef))

for w in 2*np.pi*f:
    Em[i]=Ef[i]*exp()

spectrum_final=

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

t=np.arange(-n/2,n/2)/(n*df)

pulse=np.fft.ifftshift(np.fft.ifft(spectrum))

It=abs(pulse)

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

fit=GaussianFit3(tt,Itt)

amp2=0.23*fit[3]
cen2=1.5*fit[4]
sigma2=0.5*fit[5]

amp3=0.23*fit[6]
cen3=1.5*fit[7]
sigma3=0.5*fit[8]

amp1=2.5*fit[9]
cen1=0
sigma1=1.5*fit[11] # Had to pick this value so I only use three Gaussians

t_fit=np.linspace(5*min(tt),5*max(tt),N)

sigma1=chirping(sigma1,gdd*1e-30)
sigma2=chirping(sigma2,gdd*1e-30)
sigma3=chirping(sigma3,gdd*1e-30)

It_fit=_1gaussian(t_fit,amp1,cen1,sigma1)+\
    _1gaussian(t_fit,amp2,cen2,sigma2)+\
    _1gaussian(t_fit,amp3,cen3,sigma3)
    # _1gaussian(t_fit,amp5,cen5,sigma5)
    
t=t_fit
It=It_fit

plt.plot(t_fit*1e15,It_fit,'.')
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