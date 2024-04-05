import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker

def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2)))

def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen2)/sigma2)**2)))

# linearly spaced x-axis of 10 values between 1 and 10
x_array = np.linspace(1,100,50)

# amp1 = 100
# sigma1 = 10
# cen1 = 50
# y_array_gauss = amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2)))

# # creating some noise to add the the y-axis data
# y_noise_gauss = (np.exp((np.random.ranf(50))))/5
# y_array_gauss += y_noise_gauss

amp1 = 100
sigma1 = 10
cen1 = 40

amp2 = 75
sigma2 = 5
cen2 = 65

y_array_2gauss = amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sigma1)**2))) + \
                amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen2)/sigma2)**2)))

# creating some noise to add the the y-axis data
y_noise_2gauss = (np.exp((np.random.ranf(50))))/5
y_array_2gauss += y_noise_2gauss

fig = plt.figure(figsize=(4,3))
gs = gridspec.GridSpec(1,1)
ax1 = fig.add_subplot(gs[0])

ax1.plot(x_array, y_array_2gauss, "ro")

ax1.set_xlim(-5,105)
ax1.set_ylim(-0.5,8)

ax1.set_xlabel("x_array",family="serif",  fontsize=12)
ax1.set_ylabel("y_array",family="serif",  fontsize=12)

ax1.xaxis.set_major_locator(ticker.MultipleLocator(20))
#ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))

ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

ax1.tick_params(axis='both',which='major', direction="out", top="on", right="on", bottom="on", length=8, labelsize=8)
ax1.tick_params(axis='both',which='minor', direction="out", top="on", right="on", bottom="on", length=5, labelsize=8)

fig.tight_layout()
fig.savefig("raw2Gaussian.png", format="png",dpi=1000)

popt_2gauss, pcov_2gauss = scipy.optimize.curve_fit(_2gaussian, x_array, y_array_2gauss, p0=[amp1, cen1, sigma1, \
                                                                                          amp2, cen2, sigma2])

perr_2gauss = np.sqrt(np.diag(pcov_2gauss))

pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
gauss_peak_1 = _1gaussian(x_array, *pars_1)
gauss_peak_2 = _1gaussian(x_array, *pars_2)

fig = plt.figure(figsize=(4,3))
gs = gridspec.GridSpec(1,1)
ax1 = fig.add_subplot(gs[0])

ax1.plot(x_array, y_array_2gauss, "ro")
ax1.plot(x_array, _2gaussian(x_array, *popt_2gauss), 'k--')#,\
         #label="y= %0.2f$e^{%0.2fx}$ + %0.2f" % (popt_exponential[0], popt_exponential[1], popt_exponential[2]))
    
ax1.set_xlim(-5,105)
ax1.set_ylim(-0.5,8)

ax1.set_xlabel("x_array",family="serif",  fontsize=12)
ax1.set_ylabel("y_array",family="serif",  fontsize=12)

ax1.legend(loc="best")

ax1.xaxis.set_major_locator(ticker.MultipleLocator(20))
#ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))

ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))

ax1.tick_params(axis='both',which='major', direction="out", top="on", right="on", bottom="on", length=8, labelsize=8)
ax1.tick_params(axis='both',which='minor', direction="out", top="on", right="on", bottom="on", length=5, labelsize=8)

fig.tight_layout()
fig.savefig("fit2Gaussian.png", format="png",dpi=1000)

# this cell prints the fitting parameters with their errors
print "-------------Peak 1-------------"
print "amplitude = %0.2f (+/-) %0.2f" % (pars_1[0], perr_2gauss[0])
print "center = %0.2f (+/-) %0.2f" % (pars_1[1], perr_2gauss[1])
print "sigma = %0.2f (+/-) %0.2f" % (pars_1[2], perr_2gauss[2])
print "area = %0.2f" % np.trapz(gauss_peak_1)
print "--------------------------------"
print "-------------Peak 2-------------"
print "amplitude = %0.2f (+/-) %0.2f" % (pars_2[0], perr_2gauss[3])
print "center = %0.2f (+/-) %0.2f" % (pars_2[1], perr_2gauss[4])
print "sigma = %0.2f (+/-) %0.2f" % (pars_2[2], perr_2gauss[5])
print "area = %0.2f" % np.trapz(gauss_peak_2)
print "--------------------------------"