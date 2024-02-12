import numpy as np
import matplotlib.pyplot as plt

df=0.1e13
f=np.arange(5e13,10e14,df)
n=len(f)
t=(1/(5*df*n))*np.arange(n)
If=np.exp(-(f-4e14)**2/(4*(5e13)**2))
It=np.fft.fftshift(If)

fig,(ax1,ax2)=plt.subplots(2,1,tight_layout=True)

line1,=ax1.plot(f,If)
ax1.set_xlabel(r'Frequency $f\ Hz$', fontsize=16)
ax1.set_ylabel('Amplitude', fontsize=16)
ax1.set_xlim([min(f),max(f)])

line2,=ax2.plot(t,abs(It))
ax2.set_xlabel(r'Time $t\ s$', fontsize=16)
ax2.set_ylabel('Amplitude$', fontsize=16)
ax2.set_xlim([min(t),max(t)])

plt.show()
