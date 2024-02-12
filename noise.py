import numpy as np
import matplotlib.pyplot as plt

import scipy
from scipy import signal

# Create synthetic signal
dt = 0.001
t = np.arange(0, 1, dt)
signal = np.sin(2 * np.pi * 50 * t) + np.sin(2 * np.pi * 120 * t) # Sum of 2 Sequencies
signal_clean = signal
signal = signal + 2.5 * np.random.randn(len(t)) # Add some noise
min_signal, max_signal = signal.min(), signal.max()

plt.plot(t, signal, color='c', linewidth=1.5, label='Noisy')
plt.plot(t, signal_clean, color='k', linewidth=2, label='Clean')
plt.xlim(t[0], t[-1])
plt.xlabel('t axis')
plt.ylabel('Vals')
plt.legend()
plt.show()

# Compute the Fast Fourier Transform (FFT)
n = len(t)
fhat = np.fft.fft(signal, n)                 # Compute the FFT
psd = fhat * np.conj(fhat) / n          
freq = (1 / (dt * n)) * np.arange(n)    # frequency array
idxs_half = np.arange(1, np.floor(n / 2), dtype=np.int32)  # first half index

fig, axs = plt.subplots(2, 1)

plt.sca(axs[0])
plt.plot(t, signal, color='c', linewidth=1.5, label='Noisy')
plt.plot(t, signal_clean, color='k', linewidth=2, label='Clean')
plt.xlim(t[0], t[-1])
plt.xlabel('t axis')
plt.ylabel('Vals')
plt.legend()

plt.sca(axs[1])
plt.plot(freq[idxs_half], psd[idxs_half], color='c', linewidth=2, label='PSD Noisy')
plt.xlim(freq[idxs_half[0]], freq[idxs_half[-1]])
plt.xlabel('t axis')
plt.ylabel('Vals')
plt.legend()

plt.tight_layout()
plt.show()

threshold = 100
psd_idxs = psd > threshold # array of 0 and 1
psd_clean = psd * psd_idxs # zero out all the unnecessary powers
fhat_clean = psd_idxs * fhat # used to retreive the signal

signal_filtered = np.fft.ifft(fhat_clean) # inverse fourier transform

# plt.rcParams['figure.figsize'] = [8,10]
fig, axs = plt.subplots(4, 1)

plt.sca(axs[0])
plt.plot(t, signal, color='b', linewidth=0.5, label='Noisy')
plt.plot(t, signal_clean, color='r', linewidth=1, label='Clean')
plt.ylim(min_signal, max_signal)
plt.xlabel('t axis')
plt.ylabel('Vals')
plt.legend()

plt.sca(axs[1])
plt.plot(freq[idxs_half], np.abs(psd[idxs_half]), color='b', linewidth=0.5, label='PSD noisy')
plt.xlabel('Frequencies in Hz')
plt.ylabel('Amplitude')
plt.legend()

plt.sca(axs[2])
plt.plot(freq[idxs_half], np.abs(psd_clean[idxs_half]), color='r', linewidth=1, label='PSD clean')
plt.xlabel('Frequencies in Hz')
plt.ylabel('Amplitude')
plt.legend()

plt.sca(axs[3])
plt.plot(t, signal_filtered, color='r', linewidth=1, label='Clean Signal Retrieved')
plt.xlim(t[0], t[-1])
plt.ylim(min_signal, max_signal)
plt.xlabel('t axis')
plt.ylabel('Vals')
plt.legend()

plt.tight_layout()

plt.show()
