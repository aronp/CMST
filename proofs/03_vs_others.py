import numpy as np
import matplotlib.pyplot as plt
import mpmath
from mpmath import mp
from scipy.signal.windows import kaiser, blackmanharris, hann
filename = "VsOtherFilters.pdf"

# --- 1. IEEE Plot Style Configuration ---
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 10,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'lines.linewidth': 1.0,
    'figure.figsize': (3.5, 2.8) # IEEE single column width
})

# Set 100-digit precision to see the "infinite" drop
mp.dps = 100  

# --- 2. Window Generators ---

def generate_cmst_high_precision(N, p):
    """Generates the Super-Flat CMST window: exp(t^p + 1/(t^p - 1))"""
    # Use mpmath.linspace to preserve precision
    t_vals = mpmath.linspace(-1, 1, N)
    w = []
    for t in t_vals:
        if abs(t) >= 1:
            w.append(mp.mpf(0))
        else:
            # The Super-Flat Formula
            val = mp.exp(t**p + 1/(t**p - 1))
            w.append(val)
    w_max = max(w)
    # Return as list of mpf
    return [x/w_max for x in w] 

def generate_planck(N, epsilon=0.1):
    """Generates a Planck-taper for comparison."""
    w = np.zeros(N)
    t = np.linspace(0, 1, N) 
    for i in range(N):
        if t[i] == 0 or t[i] == 1:
            w[i] = 0
        elif t[i] < epsilon:
            val = 1.0/t[i] + 1.0/(t[i]-epsilon)
            w[i] = 1.0 / (np.exp(val) + 1)
        elif t[i] > 1 - epsilon:
            val = 1.0/(1-t[i]) + 1.0/( (1-t[i])-epsilon )
            w[i] = 1.0 / (np.exp(val) + 1)
        else:
            w[i] = 1.0
    return w

def compute_spectrum_mpmath(signal):
    """High-precision DFT focusing on 0 to 0.5 normalized freq."""
    N = len(signal)
    j = mp.j
    spectrum = []
    # 500 points is enough for the curve visuals in the pdf
    freqs = np.linspace(0, 0.5, 500)
    
    for f in freqs:
        sum_val = mp.mpf(0)
        for n in range(N):
            # Centered t index
            t_idx = n - (N-1)/2 
            angle = -2 * mp.pi * f * t_idx
            sum_val += signal[n] * mp.exp(j * angle)
        
        # Floor at -1000 dB to avoid log(0) issues
        mag = mp.fabs(sum_val) + mp.mpf('1e-50') 
        spectrum.append(20 * float(mp.log10(mag)))
        
    return freqs, spectrum

# --- 3. Main Execution ---
N = 4096*2
p = 2 

print(f"Running Ultimate Comparison (N={N}, p={p})...")
print("This may take 1-2 minutes due to high precision calculations.")

# --- Generate Windows ---

# CMST (High Precision)
print("  1/5: Computing CMST (p=2)...")
w_cmst = generate_cmst_high_precision(N, p)
f_cmst, spec_cmst_raw = compute_spectrum_mpmath(w_cmst)
spec_cmst = np.array(spec_cmst_raw) - max(spec_cmst_raw)

# Kaiser (Beta 16 to match p=2 width)
print("  2/5: Computing Kaiser...")
w_kaiser = kaiser(N, 16.0)
spec_kaiser = 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(w_kaiser, 8192))))
spec_kaiser -= np.max(spec_kaiser)
f_std = np.linspace(-0.5, 0.5, 8192)
mask = (f_std >= 0) & (f_std <= 0.5)

# Blackman-Harris (4-term)
print("  3/5: Computing Blackman-Harris...")
w_bh = blackmanharris(N)
spec_bh = 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(w_bh, 8192))))
spec_bh -= np.max(spec_bh)

# Hann
print("  4/5: Computing Hann...")
w_hann = hann(N)
spec_hann = 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(w_hann, 8192))))
spec_hann -= np.max(spec_hann)

# Planck-taper
print("  5/5: Computing Planck...")
w_planck = generate_planck(N)
spec_planck = 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(w_planck, 8192))))
spec_planck -= np.max(spec_planck)


# --- 4. Plotting ---
fig, ax = plt.subplots()

# Plot Standard Windows (Thin lines)
ax.plot(f_std[mask], spec_hann[mask], 'k:', label='Hann', alpha=0.4, linewidth=0.8)
ax.plot(f_std[mask], spec_bh[mask], 'g-.', label='Blackman-Harris', linewidth=0.8, alpha=0.6)
ax.plot(f_std[mask], spec_planck[mask], 'm-.', label='Planck', linewidth=0.8, alpha=0.6)
ax.plot(f_std[mask], spec_kaiser[mask], 'r--', label='Kaiser', linewidth=1.0)

# Plot CMST (Thick, Distinct)
ax.plot(f_cmst, spec_cmst, color='#0072BD', label=f'CMST (p={p})', linewidth=1.5)

# Styling
ax.set_ylim(-500, 10) # The deep drop
ax.set_xlim(0, 0.5)
ax.set_xlabel('Normalized Frequency')
ax.set_ylabel('Magnitude (dB)')
ax.grid(True, which='major', linestyle='-', alpha=0.2)
ax.legend(loc='upper right', frameon=False, prop={'size': 6})

# Annotations
ax.annotate('Double Precision Floor\n(-320 dB)', xy=(0.2, -320), xytext=(0.15, -250),
            fontsize=7, arrowprops=dict(facecolor='red', arrowstyle='->', lw=0.5))

ax.annotate('Super-Algebraic\nDecay', xy=(0.08, -480), xytext=(0.02, -400),
            color='#0072BD', fontsize=7, 
            arrowprops=dict(color='#0072BD', arrowstyle='->', lw=0.5))
plt.savefig(filename, bbox_inches='tight', pad_inches=0.02)            
print("Success! " + filename + " Saved")
plt.show()




