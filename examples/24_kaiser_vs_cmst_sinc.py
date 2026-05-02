import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal.windows import kaiser, blackmanharris, hann

def cmst_window_p2(t):
    mask = (np.abs(t) < 1.0)
    res = np.zeros_like(t)
    res[mask] = np.exp((1) / (t[mask]**2 - 1))
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

kaiser_val = 5
w_kaiser = kaiser(N, kaiser_val)
# Kernels
sinc_kernel = np.sinc(10 * t)
cmst_filt = sinc_kernel * cmst_window_p2(t)
hann_filt = sinc_kernel * w_kaiser


# FFTs
def get_mag(x):
    fft_vals = np.fft.fft(x)
    return 20 * np.log10(np.abs(fft_vals) / np.max(np.abs(fft_vals)) + 1e-20)

freqs = np.fft.fftshift(np.fft.fftfreq(N, dt))
mag_cmst = np.fft.fftshift(get_mag(cmst_filt))
mag_hann = np.fft.fftshift(get_mag(hann_filt))

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(freqs, mag_hann, 'r--', label=f'Kaiser({kaiser_val})-Sinc (Standard)')
plt.plot(freqs, mag_cmst, 'b', label='CMST-Sinc ', linewidth=2)
plt.title(f"CMST vs. Kaiser({kaiser_val}) Sinc Window")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.xlim(-100, 100)
plt.ylim(-400, 10)
plt.grid(True, which='both', alpha=0.3)
plt.legend()



plt.savefig("Kaiser_vs_CMST_sinc.png")
print("Saved Kaiser_vs_CMST_sinc.png")
    
# Only show if NOT running on CI
if not os.environ.get('CI'):
    plt.show()
    
pass_mask = (freqs >= -4.0) & (freqs <= 4.0)
stop_mask = (freqs >= 15.0) | (freqs <= -15.0)

ripple_cmst = np.max(mag_cmst[pass_mask]) - np.min(mag_cmst[pass_mask])
ripple_kaiser = np.max(mag_hann[pass_mask]) - np.min(mag_hann[pass_mask])

rejection_cmst = np.max(mag_cmst[stop_mask])
rejection_kaiser = np.max(mag_hann[stop_mask])

# Display the results
print("\n" + "="*35)
print(f"{'Metric':<20} | {'CMST':<8} | {'Kaiser':<8}")
print("-" * 35)
print(f"{'Passband Ripple (dB)':<20} | {ripple_cmst:.4f} | {ripple_kaiser:.4f}")
print(f"{'Stopband Rejection (dB)':<20} | {rejection_cmst:.2f} | {rejection_kaiser:.2f}")
print("="*35 + "\n")    

power_cmst = 10**(mag_cmst / 10)
power_kaiser = 10**(mag_hann / 10)

total_leakage_cmst = np.trapezoid(power_cmst[stop_mask], x=freqs[stop_mask])
total_leakage_kaiser = np.trapezoid(power_kaiser[stop_mask], x=freqs[stop_mask])

total_p_cmst = np.trapezoid(power_cmst, freqs)
total_p_kaiser = np.trapezoid(power_kaiser, freqs)

leakage_ratio_cmst = 10 * np.log10(total_leakage_cmst / total_p_cmst)
leakage_ratio_kaiser = 10 * np.log10(total_leakage_kaiser / total_p_kaiser)

print(f"Total Stopband Leakage (CMST):   {leakage_ratio_cmst:.2f} dB")
print(f"Total Stopband Leakage (Kaiser): {leakage_ratio_kaiser:.2f} dB")

