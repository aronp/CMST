import os
import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import linregress
from scipy.stats import binned_statistic

def cmst_window(N):
    t = np.linspace(-1, 1, N)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = np.exp(t**4 / (t**2 - 1))
    w = np.where(np.abs(t) < 1, w, 0.0)
    return w



# Download ALL data (Quarters 1 through 17)
print("Downloading full mission data for Kepler-10...")
search = lk.search_lightcurve("Kepler-10", author="Kepler", cadence="long")
lc_collection = search.download_all()

# Stitch them together
lc_stitched = lc_collection.stitch()

# Fill Gaps & Center
lc_clean = lc_stitched.fill_gaps()
flux = lc_clean.flux.value
flux_centered = flux - np.nanmean(flux)
time = lc_clean.time.value

print(f"Total Duration: {time[-1] - time[0]:.1f} days")
print(f"New FFT Resolution: {1/(time[-1] - time[0]):.6f} cycles/day")

N = len(flux_centered)

window = cmst_window(N)
windowed_data = flux_centered * window

# RUN THE FFT
spectrum = np.fft.rfft(windowed_data)
freqs = np.fft.rfftfreq(N, d=(lc_clean.time.value[1] - lc_clean.time.value[0])) # 'd' is sample spacing in days

# FIND ALL CANDIDATE PEAKS
# We search in the range 0.5 to 16.0 cycles/day
search_mask = (freqs >= 0.5) & (freqs <= 16.0)
f_search = freqs[search_mask]
s_search = np.abs(spectrum[search_mask])

# Low threshold to catch even faint high harmonics
peaks_idx, _ = find_peaks(s_search, height=np.max(s_search)*0.05, distance=5)
candidate_freqs = f_search[peaks_idx]
candidate_amps = s_search[peaks_idx]

# THE "CHAMPION" BUCKET LOGIC
# We define the fundamental frequency (approximate)
approx_f0 = 1.1939 

# Calculate which harmonic "bucket" each peak belongs to
# e.g., 1.18 -> bucket 1,  2.41 -> bucket 2
candidate_n = np.round(candidate_freqs / approx_f0)

final_freqs = []
final_ns = []
final_amps = []

# Iterate through every unique bucket found (e.g., 1, 2, 3... 13)
unique_buckets = np.unique(candidate_n)

for n in unique_buckets:
    if n <= 0: continue # Skip DC or negative rubbish
    
    # Find all peaks that claimed to be in harmonic 'n'
    mask = (candidate_n == n)
    
    # Isolate their frequencies and amplitudes
    bucket_freqs = candidate_freqs[mask]
    bucket_amps = candidate_amps[mask]
    
    # --- THE FIX: Pick only the Highest Amplitude in this bucket ---
    champion_idx = np.argmax(bucket_amps)
    
    # Check if the champion is actually close to the predicted center (tolerance check)
    # This prevents picking a noise spike that rounded into the bucket by accident
    current_freq = bucket_freqs[champion_idx]
    if abs(current_freq - (n * approx_f0)) < 0.2: 
        final_freqs.append(current_freq)
        final_ns.append(n)
        final_amps.append(bucket_amps[champion_idx])

# Convert to arrays for regression
final_freqs = np.array(final_freqs)
final_ns = np.array(final_ns)

# 3. LEAST SQUARES REGRESSION THROUGH ORIGIN
x = final_ns
y = final_freqs

slope = np.sum(x * y) / np.sum(x**2)
intercept = 0.0

# CALCULATE STATISTICS MANUALLY
y_pred = slope * x

# R-squared Calculation
ss_res = np.sum((y - y_pred)**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

# Refined Period
refined_period = 1.0 / slope


# 4. PLOT
plt.figure(figsize=(10, 6))
plt.plot(final_ns, final_freqs, 'ro', markersize=8, label='Champion Peaks')
plt.plot(final_ns, slope*final_ns + intercept, 'k--', label=f'Fit ($R^2={r_squared:.6f}$)')
plt.title(f"Harmonic Regression (Best Peak Per Bucket)\nRefined Period: {refined_period:.6f} days")
plt.xlabel("Harmonic Index (n)")
plt.ylabel("Frequency (cycles/day)")
plt.legend()
plt.grid(True, alpha=0.3)
# Only show if NOT running on CI
if not os.environ.get('CI'):
    plt.show()

print(f"Refined Period: {refined_period:.7f} days")
print(f"R-squared: {r_squared:.6f}")



# PLOT

plt.figure(figsize=(12, 6))

plt.semilogy(freqs, np.abs(spectrum), label="CMST Spectrum")
plt.title("Frequency Domain: The Comb")
plt.xlabel("Frequency (cycles/day)")
plt.xlim(0, 14) # Zoom in on low frequencies first
plt.ylim(10**-3, 20) # Zoom in on low frequencies first


# 3. THE ANNOTATION LOOP
# We label the first peak clearly, and the others as harmonics
for i, (f, a) in enumerate(zip(final_freqs, final_amps)):
    harmonic_num = final_ns[i] # Get the integer index (1, 2, 3...)
    
    # Logic for the labels
    if harmonic_num == 1:
        text_label = "Kepler-10b\n(Fundamental)"
        y_offset_factor = 10 # Higher up for the main label
        color = 'blue'
        font_weight = 'bold'
    else:
        text_label = f"{int(harmonic_num)}f"
        y_offset_factor = 4  # Closer for the harmonics
        color = 'black'
        font_weight = 'normal'

    # The Annotation Function
    plt.annotate(
        text_label, 
        xy=(f, a),                # The point to label (The Peak)
        xytext=(f, a * y_offset_factor), # Where the text sits (Slightly above)
        textcoords='data',
        arrowprops=dict(facecolor=color, shrink=0.05, width=1, headwidth=6),
        horizontalalignment='center',
        verticalalignment='bottom',
        color=color,
        fontsize=9,
        weight=font_weight
    )

plt.legend()
plt.tight_layout()
plt.savefig('kepler_10.png', dpi=300)

# Only show if NOT running on CI
if not os.environ.get('CI'):
    plt.show()


# 1. SETUP PARAMETERS
# We assume 'time' and 'flux' are from your stitched lightkurve object
# Ensure we are using the raw values
t_all = lc_clean.time.value
f_all = lc_clean.flux.value

# 2. NORMALIZE TO PPM (Parts Per Million)
# This makes the Y-axis readable (e.g. -150 instead of 0.00015)
# We divide by the median to normalize to 1, then subtract 1, then * 1e6
f_ppm = 1e6 * (f_all / np.nanmedian(f_all) - 1)

# 3. LIMIT THE FOLD (The Critical Fix)
# We only fold the first 100 days to avoid "Drift Smearing"
# because our calculated period isn't precise enough for 4 years yet.
mask = (t_all < (t_all[0] + 100))
t_subset = t_all[mask]
f_subset = f_ppm[mask]

# 4. CALCULATE PHASE
# We add a 't0' (Epoch) estimation to center the dip
# A simple way to guess t0 for plotting is the time of the minimum flux
t0_guess = t_subset[np.argmin(f_subset)]
phase = ((t_subset - t0_guess) % refined_period) / refined_period

# Shift phase so 0.0 is in the middle (0.5)
phase = (phase + 0.5) % 1.0 - 0.5

# 5. BINNING
sort_idx = np.argsort(phase)
phase_sorted = phase[sort_idx]
f_sorted = f_subset[sort_idx]

bin_means, bin_edges, _ = binned_statistic(phase_sorted, f_sorted, statistic='mean', bins=50)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# 6. PLOT
plt.figure(figsize=(10, 6))

# Plot Raw Data (Grey)
plt.plot(phase_sorted, f_sorted, '.', color='lightgrey', alpha=0.5, markersize=2, label='Raw Data')

# Plot Binned Model (Red)
plt.plot(bin_centers, bin_means, 'r-', linewidth=3, label='Binned Signal')

plt.title(f"Kepler-10b Phase Fold\nUsing Calculated Period: {refined_period:.5f} days")
plt.xlabel("Orbital Phase (Centered)")
plt.ylabel("Flux (PPM)")
plt.ylim(-200, 100) # Zoom to see the dip
plt.xlim(-0.5, 0.5)
plt.legend()
plt.grid(True, alpha=0.3)


if not os.environ.get('CI'):
    plt.show()

