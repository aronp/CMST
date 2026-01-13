import numpy as np
import matplotlib.pyplot as plt

# Parameters for the CMST kernel
def cmst_kernel(t):
    # Standard CMST window: exp(1 / (t^2 - 1)) on (-1, 1)
    return np.where(np.abs(t) < 1, np.exp(t**4 / (t**2 - 1.0)), 0)

# 1. Generate Raw Data via FFT
N = 2**20  # High resolution for deep sidelobes
t = np.linspace(-1.5, 1.5, N)
dt = t[1] - t[0]
window = cmst_kernel(t)

# FFT
spec = np.fft.fftshift(np.fft.fft(np.fft.fftshift(window)))
freq = np.fft.fftshift(np.fft.fftfreq(N, d=dt))

# Normalize and convert to dB
mag_db = 20 * np.log10(np.abs(spec) / np.max(np.abs(spec)))

# We need distance in "bins". In a standard FFT, bins relate to the window width.
# The window width in time is T = 2 (from -1 to 1). 
# Bin spacing in frequency is df = 1/T = 0.5.
bins = freq / 0.5

# Extract envelope (positive frequencies)
mask = (bins >= 0) & (bins <= 100)
x_raw = bins[mask]
y_raw = mag_db[mask]

# 2. Theoretical Law Curve
# m = (ln(R))^2 / pi
# R = 10^(|dB| / 20)
# m = (ln(10^(|dB| / 20)))^2 / pi
# Rearranging to get dB(m):
# sqrt(pi * m) = ln(R)
# R = exp(sqrt(pi * m))
# dB = -20 * log10(exp(sqrt(pi * m)))
m_theory = np.linspace(0.1, 100, 400)
db_theory = -20 * (np.sqrt(np.pi * m_theory)) / np.log(10)

# 3. Plotting
plt.style.use('seaborn-v0_8-darkgrid')
plt.figure(figsize=(10, 6))

# Plot the Law as a dashed line
plt.plot(m_theory, db_theory, 'r--', linewidth=2, label=r'Theoretical Law: $m \approx \frac{(\ln R)^2}{\pi}$')

# Plot the Raw CMST Decay
plt.plot(x_raw, y_raw, color='navy', alpha=0.7, label='Raw CMST Spectral Decay (FFT)')

# Highlight benchmarks
benchmarks = [
    (15, -60, "10^3"),
    (27, -80, "10^4"),
    (43, -100, "10^5"),
    (72, -130, "10^6.5")
]

for m_val, db_val, r_text in benchmarks:
    plt.scatter(m_val, db_val, color='darkorange', zorder=5)
    plt.annotate(f"{db_val}dB", (m_val, db_val), textcoords="offset points", xytext=(5,5), ha='left', fontsize=9, fontweight='bold')

plt.title('CMST Spectral Decay: Raw Data vs. Theoretical Resolution Law', fontsize=14)
plt.xlabel('Distance from Main Lobe (Frequency Bins $m$)', fontsize=12)
plt.ylabel('Magnitude (dB)', fontsize=12)
plt.ylim(-150, 10)
plt.xlim(0, 100)
plt.legend(fontsize=11)
plt.grid(True, which='both', linestyle='--', alpha=0.5)

plt.savefig('cmst_decay_vs_law.png', dpi=300)

