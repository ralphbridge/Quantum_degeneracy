import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

c=299792458

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

fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)
plt.subplot(2,1,1)
line1,=ax1.plot(f,spectrum)
ax1.set_xlabel(r'Frequency $f\ Hz$', fontsize=16)
ax1.set_ylabel('Amplitude', fontsize=16)
ax1.set_xlim([min(f),max(f)])

Itmp=abs(pulse)
plt.subplot(2,1,2)
line2,=ax2.plot(t,Itmp)
ax2.set_xlabel(r'Time $t\ s$', fontsize=16)
ax2.set_ylabel('Amplitude', fontsize=16)
ax2.set_xlim([min(t),max(t)])

plt.show()

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

#fig.savefig("100fsuvtotime_nbg.pdf",bbox_inches='tight')
