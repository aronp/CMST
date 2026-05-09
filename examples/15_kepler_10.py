import os
import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import linregress
from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
from scipy.optimize import minimize


def cmst_window(N):
    t = np.linspace(-1, 1, N)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = np.exp(t**8 / (t**4 - 1))
    w = np.where(np.abs(t) < 1, w, 0.0)
    return w

def sinc_envelope(n, delta, d):
    # n: harmonic index
    # delta: scaling factor (related to transit depth)
    # d: duty cycle (D/P)
    return delta * np.abs(np.sinc(n * d))


DATA_FILENAME = "kepler10_data.csv"

if os.path.exists(DATA_FILENAME):
    print(f"Loading local cached data: {DATA_FILENAME}")
    df = pd.read_csv(DATA_FILENAME)
    time = df['time'].values
    flux = df['flux'].values
else:
    print("Local data not found. Downloading from MAST...")
    search = lk.search_lightcurve("Kepler-10", author="Kepler", cadence="long")
    lc_collection = search.download_all()
    lc_stitched = lc_collection.stitch()
    lc_clean = lc_stitched.fill_gaps()
    
    time = lc_clean.time.value
    flux = lc_clean.flux.value
    
    # Save for future stability
    df = pd.DataFrame({'time': time, 'flux': flux})
    df.to_csv(DATA_FILENAME, index=False)
    print(f"Data saved to {DATA_FILENAME}. Future runs will be identical.")
    
print("Downloading full mission data for Kepler-10...")
flux_centered = flux - np.nanmean(flux)

print(f"Total Duration: {time[-1] - time[0]:.1f} days")
print(f"New FFT Resolution: {1/(time[-1] - time[0]):.6f} cycles/day")

N = len(flux_centered)

window = cmst_window(N)
windowed_data = flux_centered * window

spectrum = np.fft.rfft(windowed_data)
freqs = np.fft.rfftfreq(N, d=(time[1] - time[0])) # 'd' is sample spacing in days

search_mask = (freqs >= 0.5) & (freqs <= 24.0)
f_search = freqs[search_mask]
s_search = np.abs(spectrum[search_mask])

peaks_idx, _ = find_peaks(s_search, prominence=np.max(s_search)*0.05, distance=100)

top_indices = np.argsort(s_search[peaks_idx])[-45:]
candidate_freqs = f_search[peaks_idx][top_indices]
candidate_amps = s_search[peaks_idx][top_indices]

# Find the fundamental from these 'Clean' peaks
# The smallest gap between these top peaks is likely the true f0
gaps = np.diff(np.sort(candidate_freqs))
auto_f0 = np.median(gaps[gaps > 0.5])

f_candidates = np.linspace(0.5, 5.0, 2000)
comb_power = []

for f_try in f_candidates:
    match_score = 0
    hits = 0
    for n in range(1, 15):
        target = n * f_try
        dist = np.abs(candidate_freqs - target)
        idx = np.argmin(dist)
        if dist[idx] < 0.1:
            match_score += candidate_amps[idx]
            hits += 1 # Track how many rungs actually had a peak
    
    # We multiply by the "Occupancy Ratio" (hits / total harmonics checked)
    # This kills sub-harmonics where every other rung is empty
    comb_power.append(match_score * (hits / 14))

approx_f0 = f_candidates[np.argmax(comb_power)]
print(f"Auto-Detected Fundamental Frequency (Occupancy Weighted): {approx_f0:.6f} cycles/day")

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
plt.xlim(0, 24) # Zoom in on low frequencies first
plt.ylim(10**-3, 20) # Zoom in on low frequencies first


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

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('kepler_10.png', dpi=300)

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

if not os.environ.get('CI'):
    plt.show()

print(f"Refined Period: {refined_period:.7f} days")
print(f"R-squared: {r_squared:.6f}")



t_all =time
f_all =flux

# Normalize to PPM
f_ppm = 1e6 * (f_all / np.nanmedian(f_all) - 1)

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

window_size = 6
kernel = np.ones(window_size) / window_size

tiled_data = np.concatenate([bin_means, bin_means, bin_means])
smoothed_tiled = np.convolve(tiled_data, kernel, mode='same')

start = len(bin_means)
end = 2 * len(bin_means)
smoothed_means = smoothed_tiled[start:end]

# Find the index of the absolute minimum
min_idx = np.argmin(bin_means)

# Extract the minimum and its two neighbors
neighbors = bin_means[max(0, min_idx-1) : min_idx+2]

# Calculate the robust average depth
mean_trough_flux = np.mean(neighbors)
depth_ppm = abs(mean_trough_flux)

print(f"Raw Stacked Depth: {depth_ppm:.2f} PPM")


plt.figure(figsize=(12, 6))

plt.plot(p_sorted, f_sorted, '.', color='lightgrey', markersize=1, alpha=0.05, label='Aligned Raw Data')

plt.plot(bin_centers, bin_means, 'r-', linewidth=1.5, alpha=0.4, label='Binned (Raw)')

plt.plot(bin_centers, smoothed_means, 'b-', linewidth=3, label=f'Moving Avg ({window_size} bins)')

plt.title(f"Kepler-10b: Drift-Corrected & Smoothed\n(Stacked 100-day segments)")
plt.xlabel("Orbital Phase")
plt.ylabel("Flux (PPM)")
plt.ylim(-150, 40)
plt.xlim(-0.5, 0.5)
plt.legend(loc='lower right')
plt.grid(True, alpha=0.3)

if not os.environ.get('CI'):
    plt.show()

n_data = final_ns
a_data = final_amps


# Initial guess: 
# d = 1/12 (based on your observation of the null)
# delta = twice the fundamental amp (rough heuristic)

popt, pcov = curve_fit(sinc_envelope, n_data, a_data, p0=[max(a_data)*2, 1/12])

fit_delta, fit_d = popt

# 4. Calculate the Physical Duration
# D = P * d
calc_duration_hours = (refined_period * fit_d) * 24

print(f"--- Sinc Fit Results ---")
print(f"Deduced Duty Cycle (d): {fit_d:.6f}")
print(f"Deduced Transit Duration: {calc_duration_hours:.2f} hours")
print(f"--- ---")

# 5. Visualization of the Fit
plt.figure(figsize=(10, 5))
plt.plot(n_data, a_data, 'ro', label='Measured Harmonic Amps')
n_smooth = np.linspace(1, max(n_data), 500)
plt.plot(n_smooth, sinc_envelope(n_smooth, *popt), 'b--', alpha=0.7, label='Sinc Envelope Fit')
plt.title(f"Fourier Envelope Fit\nCalculated Duration: {calc_duration_hours:.2f} hrs")
plt.xlabel("Harmonic Index (n)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True, alpha=0.2)

if not os.environ.get('CI'):
    plt.show()

# --- PHYSICAL CALCULATIONS ---
# Constants for the Kepler-10 System
R_STAR = 1.056       # Solar Radii
M_STAR = 0.910       # Solar Masses
G = 6.67430e-11
R_EARTH_SUN = 0.00917 # Conversion ratio


# Planetary Radius from your 115.98 PPM depth
depth_ratio = depth_ppm / 1e6
rp_rs = np.sqrt(depth_ratio)
rp_earth = (rp_rs * R_STAR) / R_EARTH_SUN

# Semi-Major Axis
M_kg = M_STAR * 1.989e30
P_sec = refined_period * 86400
a_meters = ((G * M_kg * P_sec**2) / (4 * np.pi**2))**(1/3)
a_au = a_meters / 1.496e11


# Equilibrium Temperature (Albedo = 0.3)
T_star = 5627
t_eq_c = (T_star * np.sqrt(R_STAR * 0.00465 / (2 * a_au)) * (1 - 0.3)**0.25) - 273.15


print("\n" + "="*50)
print("          KEPLER-10b: SUMMARY")
print("="*50)
print(f"{'Orbital Period:':<25} {refined_period:.6f} days")
print(f"{'Transit Depth:':<25} {depth_ppm:.0f} PPM")
print(f"{'Planetary Radius:':<25} {rp_earth:.3f} R_earth")
print("-" * 50)
print(f"{'Effective Duration:':<25} {calc_duration_hours:.3f} hours (Sinc-derived)")
print("-" * 50)
print(f"{'Orbital Distance:':<25} {a_au:.4f} AU")
print(f"{'Equilibrium Temp:':<25} {t_eq_c:.0f} °C")
print("-" * 50)




