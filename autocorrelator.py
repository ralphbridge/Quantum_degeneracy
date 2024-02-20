import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import pandas as pd

global lam,E0,sigma,eps0,c

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
#E0=np.sqrt(2*I0/(c*eps0))
tmax=0.5e-13

# Definitions for Electric fields and Intensity as functions of z and t
def E(t,alpha,E0):
    Ef=np.zeros(len(t))
    for j in range(len(t)):
        Ef=E0[j]*np.exp(-((c*t[j])/(2*sigma))**2)*np.cos((w+alpha*t[j])*t[j])
        j+=1
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

def pulse_profile():
    S=pd.read_excel("spectrum.xlsx",skiprows=1).to_numpy()
    n=np.size(S,0)

    lam=np.zeros(np.size(S,0))
    Ilam10fs=np.zeros(np.size(S,0))
    Ilam100fsuv=np.zeros(np.size(S,0))
    Ilam100fsir=np.zeros(np.size(S,0))

    for j in range(0,np.size(S,1)):
	    for i in range (0,np.size(S,0)):
		    #print(S[i][j])
		    if j==0:
    			lam[i]=(S[i][j])*1e-9
		    elif j==1:
    			Ilam10fs[i]=S[i][j]
		    elif j==2:
    			Ilam100fsuv[i]=S[i][j]
		    else:
			    Ilam100fsir[i]=S[i][j]

    Ilam10fs=(Ilam10fs-min(Ilam10fs))/max(Ilam10fs)
    Ilam100fsuv=(Ilam100fsuv-min(Ilam100fsuv))/max(Ilam100fsuv)
    Ilam100fsir=(Ilam100fsir-min(Ilam100fsir))/max(Ilam100fsir)

    f=np.zeros(n)
    If10fs=np.zeros(n)
    If100fsuv=np.zeros(n)
    If100fsir=np.zeros(n)

    for i in range(n):
        f[n-i-1]=c/lam[i]
        If10fs[n-i-1]=(lam[i]**2)*Ilam10fs[i]/c
        If100fsuv[n-i-1]=(lam[i]**2)*Ilam100fsuv[i]/c
        If100fsir[n-i-1]=(lam[i]**2)*Ilam100fsir[i]/c

    #f=np.linspace(0.5e13,10e14,n)
    df=(max(f)-min(f))/n
    t=np.arange(0,n)/(n*df)

    Iftest=np.exp(-(f-4e14)**2/(4*(0.5e13)**2))

    spectrum=If10fs

    for i in range(n): # I use this section to "clean up" the measured spectrum in lambda
        if spectrum[i]<=0.1*max(spectrum):
            spectrum[i]=0

    pulse=np.fft.ifftshift(np.fft.ifft(spectrum))

    #fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)
    #plt.subplot(2,1,1)
    #line1,=ax1.plot(f,spectrum)
    #ax1.set_xlabel(r'Frequency $f\ Hz$', fontsize=16)
    #ax1.set_ylabel('Amplitude', fontsize=16)
    #ax1.set_xlim([min(f),max(f)])

    #Itmp=abs(pulse)
    #plt.subplot(2,1,2)
    #line2,=ax2.plot(t,Itmp)
    #ax2.set_xlabel(r'Time $t\ s$', fontsize=16)
    #ax2.set_ylabel('Amplitude', fontsize=16)
    #ax2.set_xlim([min(t),max(t)])

    #plt.show()
    #fig.savefig("10fstotime.pdf",bbox_inches='tight')

    It=np.zeros(n)

    j=1
    for i in range(n):
	    if Itmp[i]>=0.005*max(Itmp):
		    It[j]=Itmp[i]
		    j+=1
		
    It=np.trim_zeros(It)
    t=t[:len(It)]

    plt.plot(t,It)
    plt.show()

    Et=np.sqrt(2*I/(c*eps0))

    return Et,t

    #fig.savefig("100fsuvtotime_nbg.pdf",bbox_inches='tight')

#t=np.linspace(-tmax,tmax,2000)

##################################################################################################################
############################################################################################### Slider

init_alpha=1.3e29

E0,t=pulse_profile()

Efield=E(t,init_alpha,E0)

# Create the figure and the line that we will manipulate
fig,(ax1,ax2)=plt.subplots(3,1,tight_layout=True)

line1,=ax1.plot(t,Efield,lw=1)
ax1.set_xlabel(r'Time $t\ s$', fontsize=16)
ax1.set_ylabel('$E\ V/m$', fontsize=16)
ax1.set_title(r'Electric field', fontsize=16, color='r')
ax1.set_xlim([-tmax,tmax])

line2,=ax2.plot(t,S_q(t,init_alpha),lw=1)
ax2.set_xlabel(r'Time difference $\tau\ s$', fontsize=16)
ax2.set_ylabel('$S_{quadratic}\ J/m^2$', fontsize=16)
ax2.set_title(r'Quadratic autocorrelator', fontsize=16, color='r')
ax2.set_xlim([-0.8*tmax,0.8*tmax])

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
    line2.set_ydata(S_q(t,alpha_slider.val))
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
# Check if number of fringes is consistent with tau=8fs
# Get theoretical expression for Fourier transform (for simple Gaussian spectrum)
# Use matlab fit function to get sine and cosine expressions(?)
# Check convergence for larger n and for larger t ranges
# Change Intensity vs frequency for Intensity vs lambda in slides 16 out of 25
# Add w+alpha*t discussion to justify why this was not an issue for Sam
# Add title on slide(s?) 17 of 25
