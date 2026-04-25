import numpy as np
import matplotlib.pyplot as plt
import h5py
import requests
import io
import json
import os
from scipy.signal import spectrogram
from scipy.interpolate import interp1d

def cmst_window(N):
    t = np.linspace(-1, 1, N)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = np.exp(t**4 / (t**2 - 1))
    w = np.where(np.abs(t) < 1, w, 0.0)
    return w

def get_gw150914_data(api_url):
    print("Fetching valid download URL from GWOSC API...")
    
    r_meta = requests.get(api_url)
    
    try:
        meta = r_meta.json()
        file_url = meta['url']
        print(f"Target Acquired: {file_url}")
    except json.JSONDecodeError:
        file_url = api_url

    print("Downloading raw HDF5 data...")
    r_file = requests.get(file_url)
    
    with h5py.File(io.BytesIO(r_file.content), 'r') as f:
        strain = f['strain/Strain'][:]
        dt = f['strain/Strain'].attrs['Xspacing']
    return strain, dt

def whiten_data(strain, dt, window):
    N = len(strain)
    spec = np.fft.rfft(strain * window(N))
    freqs = np.fft.rfftfreq(N, dt)
    
    mag = np.abs(spec)
    psd_smooth = interp1d(freqs, mag, kind='nearest', fill_value="extrapolate")(freqs)
    
    whitened_spec = spec / psd_smooth
    
    whitened_strain = np.fft.irfft(whitened_spec)
    return whitened_strain

def align_ligo_data(h1_strain, l1_strain, sample_rate=16384):
    """
    Aligns H1 and L1 data by shifting L1 by ~7ms and flipping sign.
    Ensures both arrays are returned with identical lengths.
    """
    # 1. Calculate the Shift (7ms is approximately 114 samples at 16384Hz)
    # Use 0.0069s for the precise NBI-style lag
    shift_samples = int(0.0069 * sample_rate) 
    
    # 2. Slice the arrays
    # L1 saw it first, so we skip the first 'shift' samples to 'delay' it
    l1_shifted = -1 * l1_strain[shift_samples:] 
    
    # H1 saw it last, so we drop the last 'shift' samples to align the end
    h1_shifted = h1_strain[:-shift_samples]
    
    # 3. Intersection Trim (The Robustness Guard)
    # If one file is 32.00s and another is 32.01s, this fixes it.
    min_len = min(len(l1_shifted), len(h1_shifted))
    
    h1_final = h1_shifted[:min_len]
    l1_final = l1_shifted[:min_len]
    
    return h1_final, l1_final

def plot_cmst_spectrogram(ax, data, title):
    f, t_spec, Sxx = spectrogram(data, fs, 
                                 window=cmst_window(NFFT), 
                                 nperseg=NFFT, 
                                 noverlap=noverlap)
    
    # Use log scale for Sxx to bring out the sub-threshold islands
    Sxx_log = np.log10(Sxx + 1e-250) 

    # Set the 'Brightest' point to the 99th percentile, not the Max.
    # This prevents the 60Hz line from washing out the plot.
    vmax = np.percentile(Sxx_log, 99.95) + 1

    # Set the 'Darkest' point to be 10-15 orders of magnitude below the signal.
    vmin = np.percentile(Sxx_log, 20) 
    
    mesh = ax.pcolormesh(t_spec, f, Sxx_log, shading='gouraud', cmap='magma',vmin=vmin, vmax=vmax)
    ax.contourf(t_spec, f, Sxx_log, levels=np.linspace(vmin, vmax, 60), cmap='magma', alpha=0.8)
    ax.contour(t_spec, f, Sxx_log, levels=30, colors='white', linewidths=0.3, alpha=0.4)
    
    ax.set_title(title)
    ax.set_ylabel('Frequency (Hz)')
    ax.set_ylim(0, 500)
    ax.set_xlim(zoom_center - zoom_width, zoom_center + zoom_width)
    return fig.colorbar(mesh, ax=ax, label='Log Power')



h1_api_url = "https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5"
l1_api_url = "https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3/L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"

strain_h1, dt_h1 = get_gw150914_data(h1_api_url)
strain_l1, dt_l1 = get_gw150914_data(l1_api_url)
    
    # to align the two chirps.
shift = 113
    
    # L1 was hit first, so we "delay" it to match H1
l1_aligned, h1_aligned = align_ligo_data(strain_l1,strain_h1)
    
fs = 1.0 / dt_h1
N = len(h1_aligned)
print(f"Data Loaded: {N} samples at {fs} Hz")

win_cmst = cmst_window(N)

dt = dt_h1
    # 4. FFT & Plot
print("Computing FFTs...")
spec_cmst_h1 = np.fft.rfft(h1_aligned * win_cmst)
spec_cmst_l1 = np.fft.rfft(l1_aligned * win_cmst)
freqs = np.fft.rfftfreq(N, dt)

fs = 1.0 / dt
time_vector = np.arange(len(h1_aligned)) * dt


print("Whitening the data (this extracts the chirp from the seismic noise)...")
    # We use the CMST window for the whitening step too
whitened_strain_h1 = whiten_data(h1_aligned, dt, cmst_window)
whitened_strain_l1 = whiten_data(l1_aligned, dt, cmst_window)

 # --- Spectrogram Parameters ---
NFFT = 256        
noverlap = 240    
zoom_width = 0.07  
zoom_center = 15.39



# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

# Plot Hanford
plot_cmst_spectrogram(ax1, whitened_strain_h1, 'Hanford (H1) Aligned')

# Plot Livingston
plot_cmst_spectrogram(ax2, whitened_strain_l1, 'Livingston (L1) Aligned (+7ms shift, Sign Inverted)')

ax2.set_xlabel('Time (s)')
plt.tight_layout()

plt.savefig("H1_L1_Comparison.png")
print("Saved H1_L1_Comparison.png")

if not os.environ.get('CI'):
    plt.show()
    

f_h1, t_spec_h1, Sxx_h1 = spectrogram(whitened_strain_h1, fs, 
                                 window=cmst_window(NFFT), 
                                 nperseg=NFFT, 
                                 noverlap=noverlap)

f_l1, t_spec_l1, Sxx_l1 = spectrogram(whitened_strain_l1, fs, 
                                 window=cmst_window(NFFT), 
                                 nperseg=NFFT, 
                                 noverlap=noverlap)
                                 
                                 
                                 # 1. Compute the Geometric Mean (Coherence)
# This acts as a 'logical AND' for the two detectors.
# If a feature only exists in one, it gets suppressed.
Sxx_cross = np.sqrt(Sxx_h1 * Sxx_l1)

# 2. Log-transform for visualization
# Use a deep floor to accommodate the CMST resolution
Sxx_cross_log = np.log10(Sxx_cross + 1e-50)

# 3. Set robust color limits based on percentiles
vmax = np.percentile(Sxx_cross_log, 99.9) 
vmin = vmax - 12  # Show 120dB of dynamic range

plt.figure(figsize=(12, 6))
    # Set the 'Brightest' point to the 99th percentile, not the Max.
    # This prevents the 60Hz line from washing out the plot.
vmax = np.percentile(Sxx_cross_log, 99.95) + 1

    # Set the 'Darkest' point to be 10-15 orders of magnitude below the signal.
vmin = np.percentile(Sxx_cross_log, 40) 
 
# Plot the Cross-Correlation Map
mesh = plt.pcolormesh(t_spec_h1, f_h1, Sxx_cross_log, 
                      shading='gouraud', cmap='magma', 
                      vmin=vmin, vmax=vmax)

# Add contours to highlight the 'islands'
plt.contour(t_spec_h1, f_h1, Sxx_cross_log, 
            levels=np.linspace(vmin, vmax, 20), 
            colors='white', linewidths=0.2, alpha=0.3)

plt.title('H1-L1 Coherence Map: GW150914 (CMST p=2 Kernel)')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.ylim(0, 500)
plt.xlim(15.32, 15.46) # Zoom on the event
plt.colorbar(mesh, label='Log Coherent Power')

plt.savefig("Coherence_Map.png")
plt.show()

# 1. The Coherent Average (Linear Scale)
# We sum the complex or real aligned strains to let the signals 
# add constructively and noise add destructively.
Sxx_consensus = (Sxx_h1 + Sxx_l1) / 2

# 2. Apply the Coherence Mask (The "Sober" Filter)
# We multiply by the normalized coherence to 'gate' the noise.
# This ensures that if Coherence is low, the output is pushed toward zero.
coherence_mask = np.sqrt(Sxx_h1 * Sxx_l1) / np.maximum(Sxx_h1, Sxx_l1)
Sxx_final = Sxx_consensus * coherence_mask

# 3. Log Transform for the Final Topography
Sxx_final_log = np.log10(Sxx_final + 1e-250)



plt.figure(figsize=(14, 8))

# 1. Plot the base topography (the "heat map")
mesh = plt.pcolormesh(t_spec_h1, f_h1, Sxx_final_log, 
                      shading='gouraud', cmap='inferno', # 'inferno' often contrasts well with contours
                      vmin=vmin, vmax=vmax)

# Add the Contours (The "Topographic Lines")
# Define the levels you want to draw.
# Let's draw 20 lines evenly spaced between your floor and your peak.
levels = np.linspace(vmin, vmax, 20)

contours = plt.contour(t_spec_h1, f_h1, Sxx_final_log, 
                       levels=levels, 
                       colors='white', # White contrasts well with inferno/magma
                       linewidths=0.2, # Keep them thin so they don't hide the data
                       alpha=0.4, linestyles='solid') # Slightly transparent for better integration

# Define specific 'islands' with thick contours
# This is how you prove the 15.33s precursor is a solid object.
# Let's pick a high level, say the top 10% of the signal power.
peak_contour_level = vmax - 1.6 # -15dB below the peak
plt.contour(t_spec_h1, f_h1, Sxx_final_log, 
            levels=[peak_contour_level], 
            colors='cyan', # A different color for emphasis
            linewidths=1.0, 
            alpha=0.7, linestyles='solid')

# Draw many thin background contours to show the "terrain"
# This will show the shape of everything, even the faint stuff.
bg_levels = np.linspace(vmin, vmax, 30)
plt.contour(t_spec_h1, f_h1, Sxx_final_log, levels=bg_levels, 
            colors='white', linewidths=0.1, alpha=0.2, linestyles='solid')

# Draw "Sub-threshold" highlight contours
# Instead of one peak level, we define a RANGE of interest for the islands.
# If vmax is -8.0, try looking between -9.2 and -8.5.
island_levels = [vmax - 1.2, vmax - 1.0, vmax - 0.8] 

plt.contour(t_spec_h1, f_h1, Sxx_final_log, levels=island_levels, 
            colors='cyan', linewidths=0.8, alpha=0.8,linestyles='solid')

# Final Polish
plt.title('GW150914: H1-L1 Consensus Topography with Sub-threshold Contours (CMST p=2)')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.ylim(0, 500)
plt.xlim(15.32, 15.46)
plt.colorbar(mesh, label='Log Consensus Power')

plt.savefig("Consensus_Topography_Contours.png", dpi=300)
plt.show()




