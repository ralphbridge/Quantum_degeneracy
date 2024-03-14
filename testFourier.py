import numpy as np
import matplotlib.pyplot as plt

# df=6.44e11
df=4.04e12
# df=1.02e11

#f=np.arange(2.72973e14,1.59162e15,df)
f=np.arange(2.72973e14,0.59162e15,df)
If=np.exp(-((f-4e14)/(2*1e13))**2)

n=len(f)
t=np.arange(-n/2,n/2)/(n*df)
It=abs(np.fft.ifftshift(np.fft.ifft(If)))

fig,(ax1,ax2,ax3)=plt.subplots(3,1,tight_layout=True)
plt.subplot(3,1,1)
line1,=ax1.plot(f,If)
ax1.set_xlabel(r'Frequency $f\ Hz$', fontsize=16)
ax1.set_ylabel('Amplitude in f', fontsize=16)
# ax1.set_xlim([min(lam)*1e9,max(lam)*1e9])

plt.subplot(3,1,2)
line2,=ax2.plot(t*1e15,It)
ax2.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax2.set_ylabel('Amplitude in t', fontsize=16)
#ax2.set_xlim([-100,100])

plt.subplot(3,1,3)
line3,=ax3.plot(t*1e15,It)
ax3.stem(t*1e15,It)
ax3.set_xlabel(r'Time $t\ fs$', fontsize=16)
ax3.set_ylabel('Amplitude in t', fontsize=16)
ax3.set_xlim([-100,100])

plt.show()
