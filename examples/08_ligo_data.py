import numpy as np
import matplotlib.pyplot as plt
from gwpy.timeseries import TimeSeries

# Fetch a long, stable chunk (64 seconds)
t0 = 1126259445
duration = 64
data = TimeSeries.fetch_open_data('H1', t0, t0 + duration, sample_rate=16384)

# Parameters for the Waterfall
win_len = 16 * 16384  # 16-second window for ultra-high resolution
step = 8 * 16384     # Keep 50% overlap (8-second jumps)
num_steps = (len(data.value) - win_len) // step

# Frequency Range of Interest (The 108.8 Hz Pulsar)
f_min, f_max = 330,333
freqs = np.fft.rfftfreq(win_len, d=1/16384)
mask = (freqs >= f_min) & (freqs <= f_max)

waterfall_data = []

# Generate the CMST window for this slice size
def get_cmst_p2(N):
    t = np.linspace(-1, 1, N)
    w = np.zeros(N)
    mask_win = (t > -1) & (t < 1)
    w[mask_win] = np.exp(1 + t[mask_win]**2 - 1/(1 - t[mask_win]**2))
    return w / np.max(w)

w_c2 = get_cmst_p2(win_len)
# w_c2 = get_planck_safe(win_len, eps=0.1) 

#  Slide and FFT
for i in range(num_steps):
    start = i * step
    end = start + win_len
    slice_data = data.value[start:end] * w_c2
    fft_mag = np.abs(np.fft.rfft(slice_data))
    waterfall_data.append(fft_mag[mask])

# Plotting the Waterfall
plt.figure(figsize=(12, 8))
plt.imshow(np.log10(waterfall_data), aspect='auto', 
           extent=[f_min, f_max, duration, 0], 
           cmap='magma')
plt.colorbar(label='Log Magnitude')
plt.title("CMST Archaeology: Pulsar Injection Waterfall (H1)")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Time (Seconds from t0)")
plt.show()

# Recombine by taking the median across the time axis
deep_stack_cmst = np.median(waterfall_data, axis=0)

# Plotting the Deep-Stack
plt.figure(figsize=(15, 6))
plt.semilogy(freqs[mask], deep_stack_cmst, color='darkred', lw=1.5)

plt.title("CMST Deep-Stack: 64-Second Integrated")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Integrated Magnitude")
plt.grid(True, which='both', alpha=0.3)
plt.show()


