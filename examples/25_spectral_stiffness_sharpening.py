import numpy as np
import matplotlib.pyplot as plt
import os


def cmst_window(N):
    """Generates the mathematically strict [-1, 1] CMST window."""
    t = np.linspace(-1, 1, N)
    with np.errstate(divide='ignore', invalid='ignore'):
        # Note: Using the standard root-exponential form
        w = np.exp(t**4 / (t**2 - 1.0))
    w = np.where(np.abs(t) < 1.0, w, 0.0)
    return w

def calculate_alpha(window):
    """
    Calculates the sharpening parameter alpha (kappa^2) directly 
    from the spectral footprint of the discrete window.
    """
    # Take the FFT of the pure window (no signal)
    sum = np.sum(window)
    t = np.linspace(-1, 1, len(window))
    t2 = t*t
    t2sum = np.sum(t2*window)/sum
    return t2sum

fs = 1024.0        # Sample rate
N = 2048           # Buffer size
t = np.arange(N) / fs

# Two targets 2 bins away.
f1, f2 = 150.0, 151.0 
signal = 1.0 * np.sin(2 * np.pi * f1 * t) + 0.8 * np.sin(2 * np.pi * f2 * t)

# Inject a realistic, pristine -80dB noise floor
noise_amplitude = 10**(-80/20)
noise = np.random.normal(0, noise_amplitude, N)
signal += noise

window = cmst_window(N)
# Normalize to keep 0 dB reference
window /= np.sum(window) 

spec = np.abs(np.fft.rfft(signal * window))
freqs = np.fft.rfftfreq(N, 1/fs)

# Convert to dB for physical acoustics representation
spec_db = 20 * np.log10(spec + 1e-12)

# alpha represents (kappa^2). Set analytically based on the window.
alpha = calculate_alpha(window)
# Only sharpen actual signals, ignore the noise floor
threshold_db = -60.0 

spec_sharp = np.copy(spec)

# Boolean mask: only True where the signal beats the threshold
mask = spec_db[1:-1] > threshold_db

# Vectorized 3-tap Laplacian arrays
center = spec[1:-1]
above = spec[2:]
below = spec[:-2]

# Apply the energy-preserving Laplacian equation
sharpened = 1/(1.0 - alpha) * (center - (alpha / 2.0) * above - (alpha / 2.0) * below)

# Prevent the second-derivative from driving bins into negative power
sharpened = np.maximum(sharpened, 1e-12) 

# Inject the sharpened peaks back into the spectrum, leaving noise alone
spec_sharp[1:-1] = np.where(mask, sharpened, center)
spec_sharp_db = 20 * np.log10(spec_sharp + 1e-12)

# ---------------------------------------------------------
# 4. Visualization
# ---------------------------------------------------------
plt.figure(figsize=(10, 6))

# Zoom in on the active frequency band
zoom_mask = (freqs > 140) & (freqs < 160)

plt.plot(freqs[zoom_mask], spec_db[zoom_mask], 
         label='Original CMST FFT (Fused Peak)', 
         linewidth=4, color='#1f77b4', alpha=0.4)

plt.plot(freqs[zoom_mask], spec_sharp_db[zoom_mask], 
         label=f'Sharpened FFT ($\\alpha={alpha:.3f}$)', 
         linewidth=2, color='#d62728')

# Plot true signal locations
plt.axvline(f1, color='black', linestyle='--', alpha=0.5, label=f'True f1 ({f1} Hz)')
plt.axvline(f2, color='black', linestyle=':', alpha=0.5, label=f'True f2 ({f2} Hz)')

plt.title("Frequency-Domain Laplacian Sharpening vs Standard FFT")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.ylim(-200, 5)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

# Save and show
plt.savefig("Sharpening_Before_After.png", dpi=300)

# Only show if NOT running on CI
if not os.environ.get('CI'):
    plt.show()


