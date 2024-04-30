import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker
import pandas as pd

def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2)))

def _3gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2, amp3,cen3,sigma3, amp4,cen4,sigma4, amp5,cen5,sigma5):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen2)/sigma2)**2)))+ \
            amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen3)/sigma3)**2)))+ \
            amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen4)/sigma4)**2)))+ \
            amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen5)/sigma5)**2)))
            
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

# t=np.zeros(np.size(S,0))
# I=np.zeros(np.size(S,0))

# linearly spaced x-axis of 10 values between 1 and 10
# n=200
# x_array = np.linspace(1,100,n)

# amp1 = 100
# sigma1 = 10
# cen1 = 50
# y_array_gauss = amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2)))

# # creating some noise to add the the y-axis data
# y_noise_gauss = (np.exp((np.random.ranf(50))))/5
# y_array_gauss += y_noise_gauss

# amp1 = 100
# sigma1 = 10
# cen1 = 40

# amp2 = 75
# sigma2 = 5
# cen2 = 65

# amp3 = 10
# sigma3 = 2
# cen3 = 15

# y_array_2gauss = amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + \
#                 amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen2)/sigma2)**2))) + \
#                 amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen3)/sigma3)**2)))

# creating some noise to add the the y-axis data
# y_noise_2gauss = (np.exp((np.random.ranf(n))))/5
# y_array_2gauss += y_noise_2gauss

popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_3gaussian, x_array, y_array_2gauss, p0=[amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5])

perr_2gauss = np.sqrt(np.diag(pcov_2gauss))

pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
pars_3 = popt_2gauss[6:9]
pars_4 = popt_2gauss[9:12]
# gauss_peak_1 = _1gaussian(x_array, *pars_1)
# gauss_peak_2 = _1gaussian(x_array, *pars_2)
# gauss_peak_3 = _1gaussian(x_array, *pars_3)
# gauss_peak_4 = _1gaussian(x_array, *pars_4)
# gauss_peak_5 = _1gaussian(x_array, *pars_5)

pars_1=pars_4

del pars_4

fig,(ax1,ax2)= plt.subplots(2,1,tight_layout=True)
plt.subplot(2,1,1)

myfit=2*_1gaussian(x_array,pars_1[0],pars_1[1],pars_1[2])+\
    _1gaussian(x_array,pars_2[0],pars_2[1],pars_2[2])+\
    _1gaussian(x_array,pars_3[0],pars_3[1],pars_3[2])

ax1.plot(x_array*1e15,y_array_2gauss)

         #label="y= %0.2f$e^{%0.2fx}$ + %0.2f" % (popt_exponential[0], popt_exponential[1], popt_exponential[2]))
    
# ax1.set_xlim(-5,105)
# ax1.set_ylim(-0.5,8)

ax1.set_xlabel("x_array",family="serif",  fontsize=12)
ax1.set_ylabel("y_array",family="serif",  fontsize=12)

plt.subplot(2,1,2)

ax2.plot(x_array*1e15,_3gaussian(x_array, *popt_2gauss),'.k')#,\
ax2.plot(x_array*1e15,y_array_2gauss)
         #label="y= %0.2f$e^{%0.2fx}$ + %0.2f" % (popt_exponential[0], popt_exponential[1], popt_exponential[2]))
    
# ax1.set_xlim(-5,105)
# ax1.set_ylim(-0.5,8)

ax2.set_xlabel("x_array",family="serif",  fontsize=12)
ax2.set_ylabel("y_array",family="serif",  fontsize=12)


# ax1.xaxis.set_major_locator(ticker.MultipleLocator(20))
#ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))

# ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

# ax1.tick_params(axis='both',which='major', direction="out", top="on", right="on", bottom="on", length=8, labelsize=8)
# ax1.tick_params(axis='both',which='minor', direction="out", top="on", right="on", bottom="on", length=5, labelsize=8)

fig.tight_layout()
plt.show()
# fig.savefig("fit2Gaussian.png", format="png",dpi=1000)

# this cell prints the fitting parameters with their errors
# print("-------------Peak 1-------------")
# print("amplitude = %0.2f (+/-) %0.2f" % (pars_1[0], perr_2gauss[0]))
# print("center = %0.2f (+/-) %0.2f" % (pars_1[1], perr_2gauss[1]))
# print("sigma = %0.2f (+/-) %0.2f" % (pars_1[2], perr_2gauss[2]))
# print("area = %0.2f" % np.trapz(gauss_peak_1))
# print("--------------------------------")
# print("-------------Peak 2-------------")
# print("amplitude = %0.2f (+/-) %0.2f" % (pars_2[0], perr_2gauss[3]))
# print("center = %0.2f (+/-) %0.2f" % (pars_2[1], perr_2gauss[4]))
# print("sigma = %0.2f (+/-) %0.2f" % (pars_2[2], perr_2gauss[5]))
# print("area = %0.2f" % np.trapz(gauss_peak_2))
# print("--------------------------------")
# print("-------------Peak 3-------------")
# print("amplitude = %0.2f (+/-) %0.2f" % (pars_3[0], perr_2gauss[6]))
# print("center = %0.2f (+/-) %0.2f" % (pars_3[1], perr_2gauss[7]))
# print("sigma = %0.2f (+/-) %0.2f" % (pars_3[2], perr_2gauss[8]))
# print("area = %0.2f" % np.trapz(gauss_peak_3))
# print("--------------------------------")
