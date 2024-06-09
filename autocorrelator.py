import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as scipy

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

def S_q(t,E0):
    S=np.zeros(len(t))
    j=0
    for tau in t:            
        S[j]=np.trapz((E0+E0shift(E0,j))**4,t)
        j+=1
    return S

def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))

def _3gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2, amp3,cen3,sigma3, amp4,cen4,sigma4, amp5,cen5,sigma5):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))+ \
            amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen3)/sigma3)**2)))+ \
            amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen4)/sigma4)**2)))+ \
            amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen5)/sigma5)**2)))

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

spectruml=Il

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

sigma1=1e-15
amp1=max(It)*sigma1
cen1=t[np.argmax(It)]

sigma2=5e-15
amp2=(max(It)/5)*sigma2
cen2=-2e-14

sigma3=5e-15
amp3=(max(It)/5)*sigma3
cen3=2e-14

sigma4=5e-15
amp4=(max(It)/5)*sigma4
cen4=-1e-14

sigma5=5e-15
amp5=(max(It)/5)*sigma5
cen5=1e-14

popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_3gaussian, t, It, p0=[amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5])

perr_2gauss = np.sqrt(np.diag(pcov_2gauss))

pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
pars_3 = popt_2gauss[6:9]
pars_4 = popt_2gauss[9:12]

pars_1=pars_4

del pars_4

amp1=pars_1[0]
cen1=pars_1[1]
sigma1=pars_1[2]

amp2=pars_2[0]
cen2=pars_2[1]
sigma2=pars_2[2]

amp3=pars_3[0]
cen3=pars_3[1]
sigma3=pars_3[2]

del pars_1,pars_2,pars_3

It_fit=2*_1gaussian(t,amp1,cen1,sigma1)+_1gaussian(t,amp2,cen2,sigma2)+_1gaussian(t,amp3,cen3,sigma3)

plt.plot(t*1e15,It_fit,'.k')
plt.plot(t*1e15,It)

plt.xlabel("t fs",family="serif",  fontsize=12)
plt.ylabel("Intensity",family="serif",  fontsize=12)

plt.xlim([-100,100])

plt.show()

############### 

# Efield_exp=np.multiply(E0,np.cos((w+alpha*t)*t))
n=10000
tt=np.linspace(min(t),max(t),n)
t=tt
T1=np.sqrt(1+(gdd/T0**2)**2)*T0
Efield=(E00/(np.sqrt(2*np.pi)*T1))*np.exp(-t**2/(2*T1**2))*np.cos(w*t+alpha*t**2)
Squad=S_q(t,Efield)

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
