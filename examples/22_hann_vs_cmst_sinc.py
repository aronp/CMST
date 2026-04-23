import numpy as np
import matplotlib.pyplot as plt
import os

def cmst_window_p2(t):
    mask = (np.abs(t) < 1.0)
    res = np.zeros_like(t)
    res[mask] = np.exp((t[mask]**2) / (t[mask]**2 - 1))
    return res

def hann_window(t):
    # Hann is defined as 0.5 + 0.5*cos(pi*t) on [-1, 1]
    mask = (np.abs(t) < 1.0)
    res = np.zeros_like(t)
    res[mask] = 0.5 + 0.5 * np.cos(np.pi * t[mask])
    return res

# Parameters
N = 1024 * 16
t = np.linspace(-1.5, 1.5, N)
dt = t[1] - t[0]

# Kernels
sinc_kernel = np.sinc(10 * t)
cmst_filt = sinc_kernel * cmst_window_p2(t)
hann_filt = sinc_kernel * hann_window(t)

# FFTs
def get_mag(x):
    fft_vals = np.fft.fft(x)
    return 20 * np.log10(np.abs(fft_vals) / np.max(np.abs(fft_vals)) + 1e-20)

freqs = np.fft.fftshift(np.fft.fftfreq(N, dt))
mag_cmst = np.fft.fftshift(get_mag(cmst_filt))
mag_hann = np.fft.fftshift(get_mag(hann_filt))

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(freqs, mag_hann, 'r--', label='Hann-Sinc (Standard)')
plt.plot(freqs, mag_cmst, 'b', label='CMST-Sinc ', linewidth=2)
plt.title("CMST vs. Hann Sinc Window")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.xlim(-100, 100)
plt.ylim(-200, 10)
plt.grid(True, which='both', alpha=0.3)
plt.legend()



plt.savefig("Hann_vs_CMST_sinc.png")
print("Saved Hann_vs_CMST_sinc.png")
    
# Only show if NOT running on CI
if not os.environ.get('CI'):
    plt.show()

