import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import pandas as pd
from scipy.interpolate import interp1d

global lam,E0,sigma,eps0,c,sigmat,E00

eps0=8.85e-12
c=3e8 # group velocity
vp=c # phase velocity
lam=800e-9 # laser frequency
k=2*np.pi/lam
w=k*vp
FWHM=8e-15
sigmat=FWHM/(2*np.sqrt(2*np.log(2)))
E00=np.sqrt(2*1e14/(c*eps0))

def E(t,alpha,E0):
    Ef=np.zeros(len(t))
    for j in range(len(t)):
        Ef[j]=E00*E0[j]*np.cos((w+alpha*t[j])*t[j])
        j+=1
    return Ef

def S_l(t,alpha,E0):
    S=np.zeros(len(t))
    S=2*np.trapz((E(t,alpha,E0))**2,t)
    S+=2*np.convolve(E(t,alpha,E0),E(t,alpha,E0),mode='same')
    return S

def S_q(t,alpha,E0):
    S=np.zeros(len(t))
    S=2*np.trapz((E(t,alpha,E0))**4,t)
    S+=8*np.convolve((E(t,alpha,E0))**3,E(t,alpha,E0),mode='same')
    S+=6*np.convolve((E(t,alpha,E0))**2,(E(t,alpha,E0))**2,mode='same')
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
    Iltest2=np.exp(-((lam-756e-9)/(2*8e-9))**2)+0.6*np.exp(-((lam-813e-9)/(2*20e-9))**2) # Double Gaussian

    spectruml=Il
    
    ######## Interpolating section

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
    
    del spectruml
    
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

    # for i in range(n): # I use this section to "clean up" the measured spectrum in lambda
    #     if spectrum[i]<=0.5*max(spectrum):
    #         spectrum[i]=0

    pulse=np.fft.ifftshift(np.fft.ifft(spectrum))
    
    Itmp=abs(pulse)

    fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)
    plt.subplot(2,1,1)
    line1,=ax1.plot(lam*1e9,spectruml)
    ax1.set_xlabel(r'Wavelength $\lambda\ nm$', fontsize=16)
    ax1.set_ylabel('Amplitude', fontsize=16)
    ax1.set_xlim([min(lam)*1e9,max(lam)*1e9])

    plt.subplot(2,1,2)
    line2,=ax2.plot(t*1e15,Itmp)
    ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
    ax2.set_ylabel('Amplitude', fontsize=16)
    ax2.set_xlim([-200,200])

    plt.show()
    #fig.savefig("10fstotime.pdf",bbox_inches='tight')

    It=np.zeros(n)

#     j=1
#     for i in range(n):
# 	    if Itmp[i]>=0.005*max(Itmp):
# 		    It[j]=Itmp[i]
# 		    j+=1
# 		
#     It=np.trim_zeros(It)
#     t=t[:len(It)]
    
    Itfinal=Itmp

    # plt.plot(t*1e15,Itfinal,lw=1)
    # plt.xlabel(r'Time $t\ fs$', fontsize=16)
    # plt.ylabel(r'Intensity $I(t)\ J/m^2$', fontsize=16)
    # plt.title(r'Intensity profile as function of time', fontsize=16, color='r')
    # plt.xlim([-20,20])

    # plt.show()
    
    del It, Itmp

    Et=np.sqrt(2*Itfinal/(c*eps0))
    Et=Et/max(Et)
    print(len(t))

    return (Et,t)

    #fig.savefig("100fsuvtotime_nbg.pdf",bbox_inches='tight')

##################################################################################################################
############################################################################################### Slider

init_alpha=0

E0,t=pulse_profile()

Efield=E(t,init_alpha,E0)
Slin=S_l(t,init_alpha,E0)
Squad=S_q(t,init_alpha,E0)

# Create the figure and the line that we will manipulate
#fig,(ax1,ax2,ax3)=plt.subplots(3,1,tight_layout=True)
fig,(ax1,ax3)=plt.subplots(2,1,tight_layout=True)

line1,=ax1.plot(t*1e15,Efield,lw=1)
ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax1.set_ylabel('$E\ V/m$', fontsize=16)
ax1.set_title(r'Electric field', fontsize=16, color='r')
ax1.set_xlim([-300,300])

# line2,=ax2.plot(t*1e15,Slin,lw=1)
# ax2.set_xlabel(r'Time difference $\tau\ fs$', fontsize=16)
# ax2.set_ylabel('$S_{linear}\ W/m^2$', fontsize=16)
# ax2.set_title(r'Linear detector', fontsize=16, color='r')
# ax2.set_xlim([-300,300])

line3,=ax3.plot(t*1e15 ,Squad,lw=1)
ax3.set_xlabel(r'Time difference $\tau\ fs$', fontsize=16)
ax3.set_ylabel('$S_{quadratic}\ W/m^2$', fontsize=16)
ax3.set_title(r'Quadratic detector', fontsize=16, color='r')
ax3.set_xlim([-300,300])

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

plt.show()
#fig.savefig("chirp.pdf",bbox_inches='tight')

plt.plot(t*1e15,Squad,lw=1)
plt.xlabel(r'Time $t\ fs$', fontsize=16)
plt.ylabel(r'Intensity $I(t)\ W/m^2$', fontsize=16)
plt.title(r'Autocorrelation trace (nonlinear detector)', fontsize=16, color='r')
plt.xlim([-200,200])

plt.show()

# The function to be called anytime a slider's value changes
def update(val):
    line1.set_ydata(E(t,alpha_slider.val,E0))
    line2.set_ydata(S_q(t,alpha_slider.val,E0))
    line3.set_ydata(S_q(t,alpha_slider.val,E0))
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
# IMPORTANT: Ask Herman how does time axis re-scale in this process
