import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
import pandas as pd

def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2)))

def _3gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2, amp3,cen3,sigma3, amp4,cen4,sigma4, amp5,cen5,sigma5):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))+ \
            amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen3)/sigma3)**2)))+ \
            amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen4)/sigma4)**2)))+ \
            amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen5)/sigma5)**2)))
            
S=pd.read_excel("pulse.xlsx").to_numpy()

n=np.size(S,0)

x_array=S[:,0]
y_array_2gauss=S[:,1]

sigma1=1e-15
amp1=max(y_array_2gauss)*sigma1
cen1=x_array[np.argmax(y_array_2gauss)]

sigma2=5e-15
amp2=(max(y_array_2gauss)/5)*sigma2
cen2=-2e-14

sigma3=5e-15
amp3=(max(y_array_2gauss)/5)*sigma3
cen3=2e-14

sigma4=5e-15
amp4=(max(y_array_2gauss)/5)*sigma4
cen4=-1e-14

sigma5=5e-15
amp5=(max(y_array_2gauss)/5)*sigma5
cen5=1e-14

popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_3gaussian, x_array, y_array_2gauss, p0=[amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5])

perr_2gauss = np.sqrt(np.diag(pcov_2gauss))

pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
pars_3 = popt_2gauss[6:9]
pars_4 = popt_2gauss[9:12]

pars_1=pars_4

del pars_4

fig,(ax1,ax2)= plt.subplots(2,1,tight_layout=True)
plt.subplot(2,1,1)

myfit=2*_1gaussian(x_array,pars_1[0],pars_1[1],pars_1[2])+\
    _1gaussian(x_array,pars_2[0],pars_2[1],pars_2[2])+\
    _1gaussian(x_array,pars_3[0],pars_3[1],pars_3[2])

ax1.plot(x_array*1e15,y_array_2gauss)

ax1.set_xlabel("x_array",family="serif",  fontsize=12)
ax1.set_ylabel("y_array",family="serif",  fontsize=12)

plt.subplot(2,1,2)

ax2.plot(x_array*1e15,myfit,'.k')
ax2.plot(x_array*1e15,y_array_2gauss)

ax2.set_xlabel("x_array",family="serif",  fontsize=12)
ax2.set_ylabel("y_array",family="serif",  fontsize=12)

fig.tight_layout()
plt.show()
