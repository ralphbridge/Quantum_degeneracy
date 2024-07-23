'''
Intensity FWHM

Reads in interferometric autocorrelation data and calculates the FWHM
 of the pulse assuming a Gaussian envelope.

Central wavelength should be set manually below.
 
Marshall Scott (mscott@thorlabs.com)
20170417 - Initial version
'''

import numpy as np
from numpy.fft import fft, ifft, fftshift, ifftshift
import matplotlib.pyplot as plt
import csv
from scipy import signal

def format_plot(xlabel, ylabel):
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.rc('text', usetex=True)  # Allow Latex formatting
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=14)
    return ax


def read_data(fname, skip):
    """
    read data from a csv file
    fname - file name
    skip - number of header lines to skip
    """
    with open(fname, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for i in range(skip):
            next(spamreader, None)  # skip line
        lam = []
        curr = []
        for row in spamreader:
            lam.append(float(row[3]))
            curr.append(float(row[4]))
        return np.array(lam), np.array(curr)


def FWHM(t, s):
    s = s - np.max(s) / 2.  # put half max at 0
    sgn = np.sign(s)  # convert to 1s and 0s
    d = sgn[1:] - sgn[:-1]  # take the derivative
    # find the left and right most indexes
    left_idx = np.where(d > 0)
    left_idx = left_idx[0][0]
    right_idx = np.where(d < 0)
    right_idx = right_idx[0][-1]
    return t[right_idx] - t[left_idx]  # return FWHM

lam = 0.8  # [um] Central wavelength
fname = 'sample_signal.csv'  # Autocorrelation trace
bgname = 'sample_background.csv'  # Background trace (laser blocked)

# Import data
t, bg = read_data(bgname, 18)
t, s = read_data(fname, 18)

# Subtract background
s = s - np.mean(bg)
s_max = np.max(s)

# Shift time coordinates
t = t - t[np.argmax(s)]

# Fourier transform
dt = t[1] - t[0]
f = np.linspace(-1/2./dt, 1/2./dt, len(t))
S = ifftshift(fft(fftshift(s)))*dt

# Find the fundamental frequency
half = int(len(S)/2)
Schop = np.abs(S[half:])
peakind = signal.find_peaks_cwt(Schop, np.arange(5, 30))
mxind = np.argmax(Schop[peakind[1:]])
mxind = peakind[mxind + 1] + half
mx_f = f[mxind]  # [Hz] fundamental frequency
T = 1./mx_f  # [s] period
delay = t * lam / T  # [um] optical delay
delay = delay / 0.3  # [fs] optical delay

# Filter the original signal.
S[np.abs(f) > mx_f/2.] = 0
s = fftshift(ifft(ifftshift(S)))/dt
s = np.abs(s)

s = s - np.min(s)
t_ac = FWHM(delay, s)
t_g = t_ac / np.sqrt(2)
t_s = t_ac / 1.543

print('AC FWHM: ', t_ac, ' fs')
print('Gaussian Pulse FWHM: ', t_g, ' fs')
print('Sech^2 Pulse FWHM: ', t_s, ' fs')
# print 'Max signal: ', s_max, ' V'

# Plot
ax = format_plot(r'Optical delay (fs)', r'Amplitude (AU)')
ax.plot(delay, s, color='k', ls='-', linewidth=1.5, label='x')
plt.savefig('1.png')
ax = format_plot(r'frequency (Hz)', r'Amplitude (AU)')
ax.plot(f, np.abs(S), color='k', ls='-', linewidth=1.5, label='x')
# plt.axis([0, 50, 0, 1.1*np.max(np.abs(S))])
plt.savefig('2.png')
# plt.close()
