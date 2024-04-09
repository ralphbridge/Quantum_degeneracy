import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as scipy
from scipy.interpolate import interp1d
from scipy import optimize

eps0=8.85e-12
c=3e8
lamlaser=756e-9 # laser central wavelength
k=2*np.pi/lamlaser
w=2*np.pi*c/lamlaser
FWHM=8e-15
# T0=FWHM/(2*np.sqrt(2*np.log(2)))
T0=FWHM
E00=np.sqrt(2*1e14/(c*eps0))

# gdd=10e-30 # Total GDD [(Air+mirrors+chirped mirrors+FSAC)+glass] in fs^2
gdd=0
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

def _3gaussian(x, p):
    return p[0]*(1/(p[1]*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((t-p[2])/p[1])**2))) + \
            p[3]*(1/(p[4]*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((t-p[5])/p[4])**2)))+ \
            p[6]*(1/(p[7]*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((t-p[8])/p[7])**2)))

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

Iltest1=np.exp(-((lam-756e-9)/(2*20e-9))**2) # Single Gaussian
Iltest2=np.exp(-((lam-756e-9)/(2*10e-9))**2)+0.75*np.exp(-((lam-800e-9)/(2*10e-9))**2)-0.4*np.exp(-((lam-790e-9)/(2*15e-9))**2) # Double Gaussian

spectruml=Iltest2

f=np.zeros(n)
spectrum=np.zeros(n)

for i in range(n):
    f[n-i-1]=c/(lam[i])
    spectrum[n-i-1]=(lam[i]**2)*spectruml[i]/c

df=(max(f)-min(f))/n

######## Creating equally spaced frequency domain and spectrum

# dfmax=0
# dfmin=df

# for i in range(n-1):
#     if f[i+1]-f[i]>dfmax:
#         dfmax=f[i+1]-f[i]
#     elif f[i+1]-f[i]<dfmin:
#         dfmin=f[i+1]-f[i]

# df=dfmax

# ftemp=np.arange(min(f),max(f),df)

# n=len(ftemp)
# spectrumtemp=np.zeros(n)

# spectrum_interp=interp1d(f,spectrum)
# #spectrum_interp=UnivariateSpline(f,spectrum)

# for i in range(n):
#     spectrumtemp[i]=spectrum_interp(ftemp[i])

# f=ftemp
# spectrum=spectrumtemp

# del ftemp, spectrumtemp

########### Increasing time resolution (by increasing frequency range) by 3^N

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
ax1.set_ylabel('Amplitude', fontsize=16)
ax1.set_xlim([min(lam)*1e9,max(lam)*1e9])

plt.subplot(2,1,2)
line2,=ax2.plot(t*1e15,It)
ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax2.set_ylabel('Amplitude', fontsize=16)
ax2.set_xlim([-100,100])

plt.show()

E0=np.sqrt(2*It/(c*eps0)) # check this two lines
E0=E00*E0/max(E0)

############### Gaussian fit (3 Gaussian functions)

pars=[max(E0),10e-15,t[np.argmin(E0)],max(E0)/5,20e-15,t[np.argmin(E0)]//4,max(E0)/5,20e-15,3*t[np.argmin(E0)]//4] # Amp, sigma, center (per function) in form of list
popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_3gaussian, t, E0, p0=pars)

E0fit=_3gaussian(t, *popt_2gauss)

fig = plt.figure(figsize=(4,3))

ax1.plot(t,E0, "ro")
ax1.plot(t,E0fit, 'k--')#,\
         #label="y= %0.2f$e^{%0.2fx}$ + %0.2f" % (popt_exponential[0], popt_exponential[1], popt_exponential[2]))
    
ax1.set_xlim([-100,100])
ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax1.set_ylabel('Amplitude', fontsize=16)
ax1.legend(loc="best")
ax1.tick_params(axis='both',which='major', direction="out", top="on", right="on", bottom="on", length=8, labelsize=8)
ax1.tick_params(axis='both',which='minor', direction="out", top="on", right="on", bottom="on", length=5, labelsize=8)

fig.tight_layout()
plt.show()

############### 

# Efield_exp=np.multiply(E0,np.cos((w+alpha*t)*t))
n=10000
tt=np.linspace(min(t),max(t),n)
t=tt
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
# fig.savefig("twogaussians_nochirp.pdf",bbox_inches='tight')

# fig = plt.figure()
# plt.plot(t*1e15,Squad,lw=1)
# plt.xlabel(r'Time $t\ fs$', fontsize=16)
# plt.ylabel(r'Intensity $I(t)\ W/m^2$', fontsize=16)
# plt.xlim([-100,100])

# plt.show()
# fig.savefig("10fstrace_nochirp.pdf",bbox_inches='tight')

# Check difference between calculated and experimental total GDD
# Discuss Agrawal's expressions with Herman
# Get theoretical expression for Fourier transform (for simple Gaussian spectrum)
