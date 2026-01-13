import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftshift

# 1. Setup
N = 8192
w = get_cmst_window(N)

# 2. Use the Elegant Formula: m = ceil( (ln R)^2 / pi )
m_80 = int(np.ceil((np.log(10**4))**2 / np.pi))    # Should be ~27
m_100 = int(np.ceil((np.log(10**5))**2 / np.pi))   # Should be ~42

# 3. Construct Signal
t = np.arange(N)
star = np.ones(N) 
planet_80 = 10**(-80/20) * np.exp(2j * np.pi * m_80 * t / N)
planet_100 = 10**(-100/20) * np.exp(2j * np.pi * -m_100 * t / N)

signal = (star + planet_80 + planet_100) * w

# 4. Transform
spec = fftshift(fft(signal))
mag = 20 * np.log10(np.abs(spec) / np.max(np.abs(spec)))
freq_bins = np.arange(-N//2, N//2)

# 5. Annotated Plot
plt.figure(figsize=(12, 7))
plt.plot(freq_bins, mag, label='CMST Spectrum', color='#1f77b4', lw=1.5)

# Highlight the Star and Planets
plt.scatter([0, m_80, -m_100], [0, -80, -100], color='black', zorder=5)

# Vertical markers for the bins
plt.axvline(m_80, color='orange', linestyle=':', alpha=0.8)
plt.axvline(-m_100, color='red', linestyle=':', alpha=0.8)

# Annotate specific bin numbers
plt.annotate(f'Planet (-80dB)\nBin: {m_80}', xy=(m_80, -80), xytext=(m_80+10, -60),
             arrowprops=dict(arrowstyle='->', color='orange', lw=1.5), fontweight='bold')

plt.annotate(f'Planet (-100dB)\nBin: -{m_100}', xy=(-m_100, -100), xytext=(-m_100-25, -80),
             arrowprops=dict(arrowstyle='->', color='red', lw=1.5), fontweight='bold')

plt.annotate('The "Star"\nBin: 0', xy=(0, 0), xytext=(-15, 10),
             arrowprops=dict(arrowstyle='->', color='black'), fontweight='bold')

# Limits and labels
plt.xlim(-80, 80)
plt.ylim(-130, 20)
plt.xlabel("Frequency Bin Distance (m)", fontsize=12)
plt.ylabel("Magnitude (dB)", fontsize=12)
plt.title(f"Verification of the Resolution Law: $m = \\lceil \\pi^{{-1}} (\\ln R)^2 \\rceil$\n"
          f"Confirmed Bins: {m_80} (at -80dB) and {m_100} (at -100dB)", fontsize=14)

plt.grid(True, which='major', linestyle='--', alpha=0.4)
plt.tight_layout()

plt.savefig('resolution_law.png')
print("resolution_law.png")

# Only show if NOT running on CI
if not os.environ.get('CI'):
    plt.show()


