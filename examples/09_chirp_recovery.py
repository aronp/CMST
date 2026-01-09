import numpy as np
import matplotlib.pyplot as plt

# Setup Simulation Parameters
fs = 4096*4  
T = 1.0    # 1 second of data
t = np.linspace(0, T, int(fs * T))

# Generate a Gravitational Wave "Chirp"
# Frequency starts at 30Hz and sweeps to 250Hz (typical BBH merger)
f_start = 30
f_end = 250
chirp = np.sin(2 * np.pi * (f_start * t + (f_end - f_start) * t**2 / (2 * T)))

# Apply a quick ramp-up and a ringdown-like decay to make it "real"
envelope = np.exp(-8 * (t - T)**2) 
signal = chirp * envelope

# Add "Real-World" Noise
# Instead of white noise, let's add a "Seismic Wall" (low frequency dominance)
# and a few "Spectral Lines" (60Hz hum)
white_noise = np.random.normal(0, 1.5, len(t))
seismic_noise = 5.0 * np.sin(2 * np.pi * 10 * t) # Heavy 10Hz noise
power_line = 2.0 * np.sin(2 * np.pi * 60 * t)   # 60Hz noise
data = signal + 0.5*(white_noise + seismic_noise + power_line)

# Define Windows
def cmst_p2(N):
    tau = np.linspace(-1, 1, N)
    w = np.zeros(N)
    mask = np.abs(tau) < 1
    w[mask] = np.exp(1 + tau[mask]**2 - 1/(1 - tau[mask]**2))
    return w / np.max(w)

w_kaiser = np.kaiser(len(t), 14)
w_cmst = cmst_p2(len(t))

#  Analysis
fft_kaiser = np.abs(np.fft.rfft(data * w_kaiser))
fft_cmst = np.abs(np.fft.rfft(data * w_cmst))
fft_pure = np.abs(np.fft.rfft(signal))
freqs = np.fft.rfftfreq(len(t), 1/fs)

#  Visualization
plt.figure(figsize=(12, 6))
plt.semilogy(freqs, fft_kaiser, label='Kaiser 14', color='gray', alpha=0.7)
plt.semilogy(freqs, fft_cmst, label='CMST (p=2)', color='red', linewidth=1.5)
plt.semilogy(freqs, fft_pure, label='Pure', color='blue', linewidth=1.5)

plt.axvspan(30, 250, color='green', alpha=0.1, label='Target Chirp Band')
plt.title("Chirp Recovery in High-Noise Environment (Seismic + Line Interference)")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (Log)")
plt.xlim(0, 500)
plt.legend()
plt.grid(True, which='both', alpha=0.3)
plt.savefig('chirp_simulation.png')
plt.show()
print("Saved chirp_simulation.png")
# Calculate SNR-like metric in the band
band_mask = (freqs >= 30) & (freqs <= 250)
snr_kaiser = np.mean(fft_kaiser[band_mask]) / np.std(fft_kaiser[~band_mask & (freqs < 500)])
snr_cmst = np.mean(fft_cmst[band_mask]) / np.std(fft_cmst[~band_mask & (freqs < 500)])


