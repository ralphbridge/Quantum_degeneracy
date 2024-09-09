import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as scipy
import math

eps0=8.85e-12
c=3e8
lamlaser=756e-9 # laser central wavelength
k=2*np.pi/lamlaser
w=2*np.pi*c/lamlaser
FWHM=8e-15
T0=FWHM/(2*np.sqrt(2*np.log(2)))
# T0=FWHM
E00=np.sqrt(2*1e14/(c*eps0))

# gdd=2 # Total GDD [(Air+mirrors+chirped mirrors+FSAC)+glass] in fs^2
gdd=1
# alpha=1e28 # complex part of the beam parameter Gamma
alpha=w*gdd*1e-30/T0**3

case=1

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

def _1gaussian(x, amp,cen,sigma):
    return amp*(1/(sigma*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen)/sigma)**2)))

def _3gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2, amp3,cen3,sigma3, amp4,cen4,sigma4, amp5,cen5,sigma5):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))+ \
            amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen3)/sigma3)**2)))+ \
            amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen4)/sigma4)**2)))+ \
            amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen5)/sigma5)**2)))

def GaussianFit(x_array,y_array):

    sigma=6.5e-15
    amp=max(y_array)
    cen=0

    popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_1gaussian, x_array, y_array, p0=[amp, cen, sigma])

    # perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
    
    plt.plot(x_array*1e15,y_array)
    plt.xlim([-20,20])
    plt.grid(which='both')
    plt.show()
    
    y_fit=_1gaussian(x_array, popt_2gauss[0], popt_2gauss[1], popt_2gauss[2])
    
    plt.plot(x_array*1e15,y_fit,'.')
    plt.show()
    
    return popt_2gauss

def GaussianFit3(x_array,y_array):

    sigma1=1e-15
    amp1=max(y_array)*sigma1
    cen1=x_array[np.argmax(y_array)]

    sigma2=5e-15
    amp2=(max(y_array)/5)*sigma2
    cen2=-2e-14

    sigma3=5e-15
    amp3=(max(y_array)/5)*sigma3
    cen3=2e-14

    sigma4=5e-15
    amp4=(max(y_array)/5)*sigma4
    cen4=-1e-14

    sigma5=5e-15
    amp5=(max(y_array)/5)*sigma5
    cen5=1e-14

    popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_3gaussian, x_array, y_array, p0=[amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5])

    # perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
    
    return popt_2gauss

def chirping(T0,GDD):
    T1=np.sqrt(1+(GDD/(T0**2))**2)*T0
    return T1

def roundup(x):
    return int(math.ceil(x / 100.0)) * 100
###############

if case!=1 and case!=2:
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

# Iltest1=np.exp(-((lam-756e-9)/(np.sqrt(2)*15e-9))**2) # Single Gaussian
# Iltest2=np.exp(-((lam-756e-9)/(np.sqrt(2)*20e-9))**2)+0.75*np.exp(-((lam-800e-9)/(np.sqrt(2)*10e-9))**2)-0.4*np.exp(-((lam-790e-9)/(np.sqrt(2)*15e-9))**2) # Triple Gaussian

if case==3:
    spectruml=Il
elif case==4:
    spectruml=Iluv
elif case==5:
    spectruml=Ilir

if case!=1 and case!=2:
    f=np.zeros(n)
    spectrum=np.zeros(n)

    for i in range(n):
        f[n-i-1]=c/(lam[i])
        spectrum[n-i-1]=(lam[i]**2)*spectruml[i]/c

    df=(max(f)-min(f))/n
    
    Ef=np.sqrt(2*spectrum/(c*eps0))

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

if case==3:

    tt=t[len(t)//2-60:len(t)//2+60]
    Itt=It[len(t)//2-60:len(t)//2+60]

    fit=GaussianFit3(tt,Itt)

    amp2=0.15*fit[3]
    cen2=2.5*fit[4]
    sigma2=0.8*fit[5]

    amp3=0.15*fit[6]
    cen3=2.5*fit[7]
    sigma3=0.8*fit[8]

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
elif case!=1 and case!=2:
    tt=t[len(t)//2-200:len(t)//2+200]
    Itt=It[len(t)//2-200:len(t)//2+200]

    fit=GaussianFit(tt,Itt)

    # amp=fit[0]
    cen=fit[0]
    sigma=fit[1]
    
    # print(amp)
    print(cen)
    print(sigma)

    t_fit=np.linspace(5*min(tt),5*max(tt),N)

    sigma=chirping(sigma,gdd*1e-30)

    It_fit=1/(sigma*(np.sqrt(2*np.pi)))*(np.exp((-1.0/2.0)*(((t_fit-cen)/sigma)**2)))
    
    It_fit=max(Itt)*It_fit/max(It_fit)
else:
    sigma=T0
    sigma=chirping(sigma,gdd*1e-30)
    t_fit=np.linspace(-50,50,N)*1e-15
    It_fit=1/(sigma*(np.sqrt(2*np.pi)))*(np.exp((-1.0/2.0)*(((t_fit)/sigma)**2)))
    
t=t_fit
It=It_fit

if case==3:
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
ax1.set_xlim([-20,20])
ax1.grid()

# line2,=ax2.plot(t*1e15 ,Slinear,lw=1)
# ax2.set_xlabel(r'Time delay $\tau\ fs$', fontsize=16)
# ax2.set_ylabel(r'$S_{linear}\ W/m^2$', fontsize=16)
# ax2.set_xlim([-20,20])
# ax2.grid()

line3,=ax3.plot(t*1e15 ,Squad,lw=1)
ax3.set_xlabel(r'Time delay $\tau\ fs$', fontsize=16)
ax3.set_ylabel(r'$S_{quadratic}\ W/m^2$', fontsize=16)
ax3.set_xlim([-20,20])
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