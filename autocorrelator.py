import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

global lam,sigma,eps0,c,sigmat,E00

eps0=8.85e-12
c=3e8 # group velocity
vp=c # phase velocity
lamlaser=800e-9 # laser frequency
k=2*np.pi/lamlaser
w=k*vp
FWHM=8e-15
sigmat=FWHM/(2*np.sqrt(2*np.log(2)))
E00=np.sqrt(2*1e14/(c*eps0))
#alpha=0.3e27
alpha=0

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

def pulse_profile():
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
    
    Iltest1=np.exp(-((lam-756e-9)/(2*10e-9))**2) # Single Gaussian
    Iltest2=0.925*np.exp(-((lam-756e-9)/(2*8e-9))**2)+0.5*np.exp(-((lam-813e-9)/(2*20e-9))**2) # Double Gaussian

    spectruml=Il
    
    ######## Interpolating section
    
    nint=0 # Number of interpolated iterations: N added points between experimental points is N=2^n-1

    for i in range(nint):
        lamtemp=np.zeros(2*n-1)
        spectrumltemp=np.zeros(2*n-1)
    
        spectruml_interp=interp1d(lam,spectruml)

        for i in range(n):
            lamtemp[2*i]=lam[i]
            spectrumltemp[2*i]=spectruml[i]
            if i<n-1:
                lamtemp[2*i+1]=(lam[i+1]+lam[i])/2
                spectrumltemp[2*i+1]=spectruml_interp(lamtemp[2*i+1])

        n=2*n-1
        lam=lamtemp
    
        del spectruml, spectruml_interp
    
        spectruml=spectrumltemp
    
        del spectrumltemp
    
    #################

    f=np.zeros(n)
    spectrum=np.zeros(n)

    for i in range(n):
        f[n-i-1]=c/(lam[i])
        spectrum[n-i-1]=(lam[i]**2)*spectruml[i]/c

    df=(max(f)-min(f))/n
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
    #fig.savefig("10fstotime.pdf",bbox_inches='tight')

    Et=np.sqrt(2*It/(c*eps0))
    Et=E00*Et/max(Et)

    return (lam,spectruml,Et,t)

##################################################################################################################
############################################################################################### Slider

lam,spectruml,E0,t=pulse_profile()

Efield=np.multiply(E0,np.cos((w+alpha*t)*t))
#Slin=S_l(t,alpha,E0)
Squad=S_q(t,alpha,Efield)

# Create the figure and the line that we will manipulate
#fig,(ax1,ax2,ax3)=plt.subplots(3,1,tight_layout=True)
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
#fig.savefig("chirp.pdf",bbox_inches='tight')

plt.plot(t*1e15,Squad,lw=1)
plt.xlabel(r'Time $t\ fs$', fontsize=16)
plt.ylabel(r'Intensity $I(t)\ W/m^2$', fontsize=16)
plt.title(r'Autocorrelation trace (nonlinear detector)', fontsize=16, color='r')
plt.xlim([-100,100])

plt.show()

############################################################################################### End of slider
##################################################################################################################

# Estimate dispersion for our laser: a) due to air, b) due to optics elements
# Check if number of fringes is consistent with tau=8fs
# Get theoretical expression for Fourier transform (for simple Gaussian spectrum)
# Check convergence for larger n and for larger t ranges(?)
# Change Intensity vs frequency for Intensity vs lambda in slides 16 out of 25
# Add w+alpha*t discussion to justify why this was not an issue for Sam
# Add title on slide(s?) 17 of 25
# IMPORTANT: Ask Herman how does time axis re-scale in this process
