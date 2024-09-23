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

GDD_laser=124e-30 # Check FSAC manual for correct values
TOD_laser=0
GDD_fsac=230e-30
TOD_fsac=345e-45

n_bounces=2

zd=2.06+0.05*(n_bounces-1) # Distance from laser pinhole to FSAC pinhole
zw=0e-3 # Measure this again

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

def InverseFourier(Fw,w,t):
    ft=np.zeros(np.size(w),dtype=np.complex128)
    for j in range(len(t)):
        for i in range(len(w)):
            ft[j]+=Fw[i]*np.exp(1j*w[i]*t[j])+np.conj(Fw[i])*np.exp(-1j*w[i]*t[j])
    return np.real(ft)

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
        S[j]=np.trapezoid((E0+E0shift(E0,j))**2,t)
        j+=1
    return S

def S_q(t,E0):
    S=np.zeros(len(t))
    j=0
    for tau in t:            
        S[j]=np.trapezoid((E0+E0shift(E0,j))**4,t)
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

Il=(Il-min(Il))/max(Il) # Initial measured spectrum

f=np.zeros(n)
spectrumf=np.zeros(n)

for i in range(n):
    f[n-i-1]=c/(lam[i])
    spectrumf[n-i-1]=(lam[i]**2)*Il[i]/c

######################## Increasing time resolution (by increasing frequency range)

Ef0=np.sqrt(2*spectrumf/(c*eps0))
t0=np.linspace(0,50e-15,len(f))
Et0=InverseFourier(Ef0,2*np.pi*f,t0)

########################## Trimming spectrum low and high frequencies and interpolating it to get equally spaced frequencies

NN=799 # Length of actual spectrum (obtained from matlab)
f_trim=np.zeros(NN)
spectrumf_trim=np.zeros(NN)

for i in range(0,NN):
    f_trim[i]=f[i+266] # 266 is the position of the original spectrum where it starts being bigger than 0
    spectrumf_trim[i]=spectrumf[i+266] # New angular frequency spectrum (trimmed)

spectrumf_interp_function=inter.CubicSpline(f_trim,spectrumf_trim)

N=6900
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

# Initial phases for original measured spectrum (up to third order)

for i in range(N):
    # Dispersion phases introduced by the laser itself to the pulse
    Ef0_interp[i]=Ef0_interp[i]*np.exp(1j*GDD_laser*(w[i]-w0)**2/math.factorial(2))
    Ef0_interp[i]=Ef0_interp[i]*np.exp(1j*TOD_laser*(w[i]-w0)**3/math.factorial(3))

######################################

t0_interp=np.linspace(0,50e-15,len(f_interp))
Et0_interp=InverseFourier(Ef0_interp,2*np.pi*f_interp,t0_interp)

# fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)
# plt.subplot(2,1,1)
# line1,=plt.plot(t0*1e15,Et0)
# ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
# ax1.set_ylabel(r'$E_0\ V/m$', fontsize=16)
# ax1.grid(which='both')

# plt.subplot(2,1,2)
# line2,=ax2.plot((t0_interp)*1e15,Et0_interp)
# ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
# ax2.set_ylabel('$E_0^{interp}\ V/m$', fontsize=16)
# ax2.grid(which='both')
# plt.show()

########################## Getting GDD data from Thorlabs chirped mirrors data

CM=pd.read_excel("UMxx-15FS_data.xlsx").to_numpy()
n_cm=np.size(CM,0)
lam_cm=np.zeros(n_cm)
GD_cm_lam=np.zeros(n_cm)

for j in range(np.size(CM,1)):
    for i in range (np.size(CM,0)):
        if j==0:
            lam_cm[i]=(CM[i][j])*1e-9
        else:
            GD_cm_lam[i]=(CM[i][j])*1e-15

w_cm=np.zeros(n_cm)
GD_cm=np.zeros(n_cm)
for i in range(n_cm):
    w_cm[n_cm-1-i]=2*np.pi*c/lam_cm[i]
    GD_cm[n_cm-1-i]=GD_cm_lam[i]

GD_cm_interp_function=inter.CubicSpline(w_cm,GD_cm)
GD_cm_interp=GD_cm_interp_function(w)
# GD_cm_interp=np.interp(w,w_cm,GD_cm)

GDD_cm=np.zeros(N)
TOD_cm=np.zeros(N)
for i in range(N): # Higher accuracy order derivatives from Fornberg 1988
    if i==0:
        m2=(GD_cm_interp[i+1]-GD_cm_interp[i])/(w[i+1]-w[i])
        GDD_cm[i]=m2
    elif i==N-1:
        m1=(GD_cm_interp[i]-GD_cm_interp[i-1])/(w[i]-w[i-1])
        GDD_cm[i]=m1
    elif i==1 or i==N-2:
        GDD_cm[i]=(-GD_cm_interp[i-1]+GD_cm_interp[i+1])/(2*2*np.pi*df)
        TOD_cm[i]=(GD_cm_interp[i-1]-2*GD_cm_interp[i]+GD_cm_interp[i+1])/((2*np.pi*df)**2)
    else:
        GDD_cm[i]=(GD_cm_interp[i-2]-8*GD_cm_interp[i-1]+8*GD_cm_interp[i+1]-GD_cm_interp[i+2])/(12*2*np.pi*df)
        TOD_cm[i]=(-GD_cm_interp[i-2]+16*GD_cm_interp[i-1]-30*GD_cm_interp[i]+16*GD_cm_interp[i+1]-GD_cm_interp[i+2])/(12*(2*np.pi*df)**2)

########################## Getting GDD data from Thorlabs P01 mirrors data

P01=pd.read_excel("P01_data.xlsx").to_numpy()
n_p01=np.size(P01,0)
lam_p01=np.zeros(n_p01)
GDD_p01_data_lam=np.zeros(n_p01)

for j in range(0,np.size(P01,1)):
    for i in range (0,np.size(P01,0)):
        if j==0:
            lam_p01[i]=(P01[i][j])*1e-9
        else:
            GDD_p01_data_lam[i]=(P01[i][j])*1e-30

w_p01=np.zeros(n_p01)
GDD_p01_data=np.zeros(n_p01)
for i in range(n_p01):
    w_p01[n_p01-1-i]=2*np.pi*c/lam_p01[i]
    GDD_p01_data[n_p01-1-i]=GDD_p01_data_lam[i]

GDD_p01_interp_function=inter.CubicSpline(w_p01,GDD_p01_data)
GDD_p01_interp=GDD_p01_interp_function(w)

TOD_p01_data=np.zeros(n_p01)
for i in range(0,n_p01):
    if i==0:
        m2=(GDD_p01_data[i+1]-GDD_p01_data[i])/(w_p01[i+1]-w_p01[i])
        TOD_p01_data[i]=m2
    elif i==n_p01-1:
        m1=(GDD_p01_data[i]-GDD_p01_data[i-1])/(w_p01[i]-w_p01[i-1])
        TOD_p01_data[i]=m1
    else:
        m1=(GDD_p01_data[i]-GDD_p01_data[i-1])/(w_p01[i]-w_p01[i-1])
        m2=(GDD_p01_data[i+1]-GDD_p01_data[i])/(w_p01[i+1]-w_p01[i])
        TOD_p01_data[i]=(m1+m2)/2

TOD_p01_interp_function=inter.CubicSpline(w_p01,GDD_p01_data)
TOD_p01_interp=TOD_p01_interp_function(w)

#print(GDD_p01_interp(2*np.pi*c/800e-9)*1e30) # Check this line, it does not interpolate correctly

############################## Computing the phases due to each element in the layout

Ef=np.zeros(np.size(Ef0_interp),dtype=np.complex128)

for i in range(N):
    # Dispersion phases introduced by air up to the third order
    # Ef[i]=Ef0_interp[i]*np.exp(1j*k0_air*zd)
    # Ef[i]=Ef[i]*np.exp(1j*kp0_air*(w[i]-w0)*zd)
    Ef[i]=Ef0_interp[i]*np.exp(1j*kpp0_air*(w[i]-w0)**2*zd/math.factorial(2))
    # Ef[i]=Ef[i]*np.exp(1j*kppp0_air*(w[i]-w0)**3*zd/math.factorial(3))
    
    # Dispersion phases introduced by BK7 Fused Silica glass window up to the third order
    # Ef[i]=Ef[i]*np.exp(1j*k0_bk7*zw)
    # Ef[i]=Ef[i]*np.exp(1j*kp0_bk7*(w[i]-w0)*zw)
    # Ef[i]=Ef[i]*np.exp(1j*kpp0_bk7*(w[i]-w0)**2*zw/math.factorial(2))
    # Ef[i]=Ef[i]*np.exp(1j*kppp0_bk7*(w[i]-w0)**3*zw/math.factorial(3))
    
    # Dispersion phases introduced by bounces off Chirped mirrors
    # Ef[i]=Ef[i]*np.exp(1j*GDD_cm[i]*(w[i]-w0)**2/math.factorial(2))
    # Ef[i]=Ef[i]*np.exp(1j*TOD_cm[i]*(w[i]-w0)**3/math.factorial(3))
    Ef[i]=Ef[i]*np.exp(1j*n_bounces*GD_cm_interp[i]*(w[i]-w0))
    
    # # Dispersion phases introduced by bounces off P01 Silver mirrors
    # Ef[i]=Ef[i]*np.exp(1j*GDD_p01[i]*(w[i]-w0)**2/math.factorial(2))
    # Ef[i]=Ef[i]*np.exp(1j*TOD_p01[i]*(w[i]-w0)**3/math.factorial(3))
    Ef[i]=Ef[i]*np.exp(1j*GDD_p01_interp[i]*(w[i]-w0)**2/math.factorial(2))
    
    # # Dispersion phases introduced by the FSAC up to the third order
    Ef[i]=Ef[i]*np.exp(1j*GDD_fsac*(w[i]-w0)**2/math.factorial(2))
    Ef[i]=Ef[i]*np.exp(1j*TOD_fsac*(w[i]-w0)**3/math.factorial(3))

##################################################

# fs=max(f_interp)-min(f_interp) # Sampling frequency (equals the bandwidth or maximum frequency)
t=np.linspace(-3e-12,3e-12,len(f_interp))
Et=InverseFourier(Ef,2*np.pi*f_interp,t)

# fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)
# plt.subplot(2,1,1)
# line1,=plt.plot(t0*1e15,Et0)
# ax1.set_xlim([-20,20])
# ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
# ax1.set_ylabel(r'$E_0\ V/m$', fontsize=16)
# # major_tick = np.arange(roundup(min(lam)*1e9), roundup(max(lam)*1e9),200)#[200, 400, 600, 800, 1000]
# # minor_tick = np.arange(roundup(min(lam)*1e9)+100, roundup(max(lam)*1e9),200)#[300, 500, 700, 900]
# # ax1.set_xticks(major_tick) # Grid
# # ax1.set_xticks(minor_tick, minor=True)
# ax1.grid(which='both')

# plt.subplot(2,1,2)
# line2,=ax2.plot(t*1e15,Et)
# ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
# ax2.set_ylabel('$E\ V/m$', fontsize=16)
# # ax2.set_xlim([-100,100])
# # major_tick = np.arange(-100,100,20)
# # minor_tick = np.arange(-100,100,10)
# # ax2.set_xticks(major_tick)
# # ax2.set_xticks(minor_tick, minor=True)
# ax2.grid(which='both')

# # fig.savefig("Frequency_Time.pdf",bbox_inches='tight')
# plt.show()

###############

# Slinear=S_l(t,Et)
Squad=S_q(t,Et)

fig,(ax3)=plt.subplots(1,1,tight_layout=True)

# line1,=ax1.plot(t*1e15,Et,lw=1)
# ax1.plot(t*1e15,Et,'.',alpha=0.01)
# ax1.set_xlabel(r'Time $t\ fs$', fontsize=16)
# ax1.set_ylabel(r'$E\ V/m$', fontsize=16)
# ax1.set_xlim([-50,50])
# ax1.grid()

# line2,=ax2.plot(t*1e15 ,Slinear,lw=1)
# ax2.set_xlabel(r'Time delay $\tau\ fs$', fontsize=16)
# ax2.set_ylabel(r'$S_{linear}\ W/m^2$', fontsize=16)
# ax2.set_xlim([-20,20])
# ax2.grid()

line3,=ax3.plot((t-t[np.argmax(Squad)])*1e15,Squad,lw=1)
ax3.set_xlabel(r'Time delay $\tau\ fs$', fontsize=16)
ax3.set_ylabel(r'$S_{quadratic}\ W/m^2$', fontsize=16)
ax3.set_xlim([-250,250])
ax3.grid()
plt.savefig('FieldTrace_'+str(n_bounces)+'b_far.pdf',bbox_inches='tight')
plt.show()

# Check first, second and third order phases added due to air and BK7
