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



print("Downloading full mission data for Kepler-10...")
search = lk.search_lightcurve("Kepler-10", author="Kepler", cadence="long")
lc_collection = search.download_all()

lc_stitched = lc_collection.stitch()

lc_clean = lc_stitched.fill_gaps()
flux = lc_clean.flux.value
flux_centered = flux - np.nanmean(flux)
time = lc_clean.time.value

print(f"Total Duration: {time[-1] - time[0]:.1f} days")
print(f"New FFT Resolution: {1/(time[-1] - time[0]):.6f} cycles/day")

N = len(flux_centered)

window = cmst_window(N)
windowed_data = flux_centered * window

spectrum = np.fft.rfft(windowed_data)
freqs = np.fft.rfftfreq(N, d=(lc_clean.time.value[1] - lc_clean.time.value[0])) # 'd' is sample spacing in days

search_mask = (freqs >= 0.5) & (freqs <= 16.0)
f_search = freqs[search_mask]
s_search = np.abs(spectrum[search_mask])

peaks_idx, _ = find_peaks(s_search, height=np.max(s_search)*0.05, distance=5)
candidate_freqs = f_search[peaks_idx]
candidate_amps = s_search[peaks_idx]

approx_f0 = 1.1939 

candidate_n = np.round(candidate_freqs / approx_f0)

final_freqs = []
final_ns = []
final_amps = []

unique_buckets = np.unique(candidate_n)

for n in unique_buckets:
    if n <= 0: continue # Skip DC or negative rubbish
    
    mask = (candidate_n == n)
    
    bucket_freqs = candidate_freqs[mask]
    bucket_amps = candidate_amps[mask]
    
    champion_idx = np.argmax(bucket_amps)
    
    current_freq = bucket_freqs[champion_idx]
    if abs(current_freq - (n * approx_f0)) < 0.2: 
        final_freqs.append(current_freq)
        final_ns.append(n)
        final_amps.append(bucket_amps[champion_idx])

final_freqs = np.array(final_freqs)
final_ns = np.array(final_ns)

x = final_ns
y = final_freqs

slope = np.sum(x * y) / np.sum(x**2)
intercept = 0.0

y_pred = slope * x

ss_res = np.sum((y - y_pred)**2)
ss_tot = np.sum((y - np.mean(y))**2)
r_squared = 1 - (ss_res / ss_tot)

refined_period = 1.0 / slope





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


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic

# --- 1. SETUP PARAMETERS ---
# We assume 'lc_clean' and 'refined_period' are already defined from previous steps
t_all = lc_clean.time.value
f_all = lc_clean.flux.value

# Normalize to PPM
f_ppm = 1e6 * (f_all / np.nanmedian(f_all) - 1)

# --- 2. DRIFT CORRECTION (STACKING) ---
# We slice the data into 100-day chunks and re-center the dip locally
chunk_size = 100.0 # days
t_start = t_all[0]
t_end = t_all[-1]

corrected_phases = []
aligned_fluxes = []

print(f"Stacking 100-day chunks to correct drift...")

current_t = t_start
while current_t < t_end:
    mask = (t_all >= current_t) & (t_all < current_t + chunk_size)
    
    if np.sum(mask) < 100: 
        current_t += chunk_size
        continue
        
    t_chunk = t_all[mask]
    f_chunk = f_ppm[mask]
    
    local_phase = ((t_chunk - t_chunk[0]) % refined_period) / refined_period
    
    bin_means_temp, bin_edges_temp, _ = binned_statistic(local_phase, f_chunk, statistic='mean', bins=50)
    bin_centers_temp = 0.5 * (bin_edges_temp[1:] + bin_edges_temp[:-1])
    
    min_idx = np.argmin(bin_means_temp)
    phase_offset = bin_centers_temp[min_idx]
    
    shifted_phase = local_phase - phase_offset
    shifted_phase = (shifted_phase + 0.5) % 1.0 - 0.5
    
    corrected_phases.append(shifted_phase)
    aligned_fluxes.append(f_chunk)
    
    current_t += chunk_size

final_phase = np.concatenate(corrected_phases)
final_flux = np.concatenate(aligned_fluxes)

sort_idx = np.argsort(final_phase)
p_sorted = final_phase[sort_idx]
f_sorted = final_flux[sort_idx]

bin_means, bin_edges, _ = binned_statistic(p_sorted, f_sorted, statistic='mean', bins=100)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

window_size = 12
kernel = np.ones(window_size) / window_size

tiled_data = np.concatenate([bin_means, bin_means, bin_means])
smoothed_tiled = np.convolve(tiled_data, kernel, mode='same')

start = len(bin_means)
end = 2 * len(bin_means)
smoothed_means = smoothed_tiled[start:end]

plt.figure(figsize=(12, 6))

plt.plot(p_sorted, f_sorted, '.', color='lightgrey', markersize=1, alpha=0.05, label='Aligned Raw Data')

plt.plot(bin_centers, bin_means, 'r-', linewidth=1.5, alpha=0.4, label='Binned (Raw)')

plt.plot(bin_centers, smoothed_means, 'b-', linewidth=3, label=f'Moving Avg ({window_size} bins)')

plt.title(f"Kepler-10b: Drift-Corrected & Smoothed\n(Stacked 100-day segments)")
plt.xlabel("Orbital Phase")
plt.ylabel("Flux (PPM)")
plt.ylim(-100, 40)
plt.xlim(-0.5, 0.5)
plt.legend(loc='lower right')
plt.grid(True, alpha=0.3)

if not os.environ.get('CI'):
    plt.show()

