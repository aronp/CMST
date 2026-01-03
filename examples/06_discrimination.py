import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# 1. Define the Window Functions
# ---------------------------------------------------------

def cmst_window(N, p=2):
    """
    Generates the CMST window.
    p=2: Maximum spectral purity (The "Silencer").
    """
    # Normalized time grid (-1, 1)
    t = np.linspace(-1, 1, N)
    
    # Avoid singularity at exactly +/- 1
    mask = np.abs(t) < 1.0
    w = np.zeros(N)
    
    # CMST Formula
    val = t[mask]
    # The geometric mollifier: exp(1 + t^p - 1/(1-t^p))
    exponent = 1 + val**p - 1.0/(1.0 - val**p)
    w[mask] = np.exp(exponent)
    return w

def get_planck(N, epsilon=0.1):
    """
    Generates the standard Planck-taper window (LIGO standard).
    """
    w = np.zeros(N)
    for i in range(N):
        t = i / (N - 1)
        if t <= 0 or t >= 1: 
            w[i] = 0
        elif t < epsilon:
            val = (epsilon / t) - (epsilon / (epsilon - t))
            w[i] = 1.0 / (1.0 + np.exp(val))
        elif t >= 1 - epsilon:
            val = (epsilon / (1 - t)) - (epsilon / (epsilon - (1 - t)))
            w[i] = 1.0 / (1.0 + np.exp(val))
        else: 
            w[i] = 1.0
    return w

# ---------------------------------------------------------
# 2. Setup the "Lighthouse" Scenario
# ---------------------------------------------------------
N = 2048
fs = 10000            # Sample rate: 10 kHz
t = np.arange(N) / fs

# Target Parameters
f_strong = 1000       # Strong signal (Carrier)
f_weak   = 1300       # Weak signal (Target)
dB_target = -100      # The Challenge: Detect a -100 dB signal

# Calculate linear amplitude: 10^(-100/20) = 1e-5
amp_weak = 10**(dB_target / 20)

print(f"Generating signal: Strong at {f_strong}Hz, Weak at {f_weak}Hz ({dB_target} dB)")

# Create Signal
signal = 1.0 * np.sin(2*np.pi*f_strong*t) + \
         amp_weak * np.sin(2*np.pi*f_weak*t)

# ---------------------------------------------------------
# 3. Process and Plot
# ---------------------------------------------------------

# Apply Windows
w_cmst = cmst_window(N, p=2)
w_planck = get_planck(N, epsilon=0.1)

# Compute FFT (normalized to 0 dB peak)
fft_cmst = 20 * np.log10(np.abs(np.fft.fft(signal * w_cmst)) + 1e-15)
fft_planck = 20 * np.log10(np.abs(np.fft.fft(signal * w_planck)) + 1e-15)

fft_cmst -= np.max(fft_cmst)
fft_planck -= np.max(fft_planck)

# Frequency Axis
freqs = np.fft.fftfreq(N, 1/fs)
mask = freqs >= 0

# Plotting
plt.figure(figsize=(10, 6))

plt.plot(freqs[mask], fft_cmst[mask], label='CMST Window (p=2)', color='blue', linewidth=1.5)
plt.plot(freqs[mask], fft_planck[mask], label='Planck-taper (LIGO)', color='red', linestyle='--', linewidth=1.5)

# Formatting
plt.ylim(-160, 10)
plt.xlim(0, 2500)
plt.ylabel("Magnitude (dB)")
plt.xlabel("Frequency (Hz)")
plt.title(f"Weak Signal Detection Test ({dB_target} dB Target)\nCMST vs Planck-taper")
plt.legend(loc='upper right')
plt.grid(True, alpha=0.3)

# Highlight the weak signal for clarity
plt.annotate(f'Weak Signal ({dB_target} dB)', 
             xy=(f_weak, dB_target), 
             xytext=(f_weak+200, dB_target+10),
             arrowprops=dict(facecolor='black', shrink=0.05))

plt.tight_layout()
plt.savefig('weak_signal_comparison.png', dpi=300)
plt.savefig('weak_signal_comparison.pdf', dpi=300)
plt.show()
