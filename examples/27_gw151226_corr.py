import os
import requests
import h5py
import io
import numpy as np
import matplotlib.pyplot as plt
import requests
import json
from scipy.signal import spectrogram
from scipy.interpolate import interp1d

h1_api_url = "https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5"
l1_api_url = "https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3/L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"
GW151226_h1_url = "https://gwosc.org/GW151226data/H-H1_LOSC_16_V2-1135136334-32.hdf5"
GW151226_l1_url = "https://gwosc.org/GW151226data/L-L1_LOSC_16_V2-1135136334-32.hdf5"

# Execution
h1_GW151226_local = "H1_GW151226.hdf5"
l1_GW151226_local = "L1_GW151226.hdf5"


def get_gw150914_data(api_url, local_filename):
    # Check if the file already exists locally
    if os.path.exists(local_filename):
        print(f"Loading {local_filename} from local storage...")
        with h5py.File(local_filename, 'r') as f:
            strain = f['strain/Strain'][:]
            dt = f['strain/Strain'].attrs['Xspacing']
        return strain, dt

    # If not local, download and save
    print(f"{local_filename} not found. Downloading from GWOSC...")
    r_meta = requests.get(api_url)
    try:
        meta = r_meta.json()
        file_url = meta['url']
    except (requests.JSONDecodeError, KeyError):
        file_url = api_url

    r_file = requests.get(file_url)
    
    # Save to local file
    with open(local_filename, 'wb') as f_out:
        f_out.write(r_file.content)
    print(f"Saved to {local_filename}")

    # Return the data
    with h5py.File(io.BytesIO(r_file.content), 'r') as f:
        strain = f['strain/Strain'][:]
        dt = f['strain/Strain'].attrs['Xspacing']
    return strain, dt

def cmst_window(N):
    t = np.linspace(-1, 1, N)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = np.exp(t**8 / (t**4 - 1))
    w = np.where(np.abs(t) < 1, w, 0.0)
    return w

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


def calculate_alpha(window):
    """
    Calculates the sharpening parameter alpha directly 
    from the spectral footprint of the discrete window.
    """
    sum_w = np.sum(window)
    t = np.linspace(-1, 1, len(window))
    t2 = t*t
    return np.sum(t2*window)/sum_w

def sharpen_spectrogram(Sxx, alpha):
    laplacian_f = np.zeros_like(Sxx)
    # Approximate 2nd derivative along frequency axis (axis 0)
    laplacian_f[1:-1, :] = Sxx[2:, :] - 2*Sxx[1:-1, :] + Sxx[0:-2, :]
    Sxx_sharp = Sxx - alpha * laplacian_f
    
    # Clip to prevent negative values (non-physical energy)
    Sxx_sharp = np.maximum(Sxx_sharp, 1e-20)
    return Sxx_sharp

def plot_grid(Sxx_sharp,low = 20):  
    plt.figure(figsize=(14, 8))

    # This prevents the 60Hz line from washing out the plot.
    vmax = np.percentile(Sxx_sharp, 99.9) + (np.abs(np.percentile(Sxx_sharp, 99.9)) * 0.01)
    # Set the 'Darkest' point to be 10-15 orders of magnitude below the signal.
    vmin = np.percentile(Sxx_sharp, low) 
    levels = np.linspace(vmin, vmax, 20)
    trange = np.max(t_spec) - np.min(t_spec)
    plt.pcolormesh(t_spec, f, Sxx_sharp, shading='gouraud', cmap='viridis',vmin=vmin, vmax=vmax)
    contour_filled = plt.contourf(t_spec, f, Sxx_sharp, levels=levels, cmap='inferno')
    contour_lines = plt.contour(t_spec, f, Sxx_sharp, levels=levels, colors='white', linewidths=0.5, alpha=0.5)

    # Zoom in on the Chirp
    plt.xlim(trange/2 - zoom_width, trange/2 + zoom_width)
    plt.ylim(0, 1000) # The "Audible" range of the black holes

    plt.title(f'LIGO GW151226: Spectrogram Analysis (CMST Window)')
    plt.ylabel('Frequency (Hz)')
    plt.xlabel('Time (s)')
    plt.colorbar(contour_filled, label='Signal-to-Noise Ratio')
    plt.grid(True)
    plt.show()


def generate_band_limited_map(h1, l1, freqs, freq_limit=500):
    """
    Computes correlation only for frequencies up to freq_limit.
    """
    # 1. Identify the index for 500Hz
    idx_limit = np.searchsorted(freqs, freq_limit)
    
    # 2. Slice the spectrograms to focus on the 0-500Hz band
    h1_slice = h1[:idx_limit, :]
    l1_slice = l1[:idx_limit, :]
    
    # 3. Standardize only the relevant data
    h1_n = (h1_slice - np.mean(h1_slice)) / np.std(h1_slice)
    l1_n = (l1_slice - np.mean(l1_slice)) / np.std(l1_slice)
    
    # 4. Compute the local correlation map (zero-lag product)
    # This results in values between -1 and 1 if normalized correctly
    corr_slice = h1_n * l1_n
    
    # 5. Reconstruct the full-size map with zeros for irrelevant frequencies
    # This prevents high-frequency noise from affecting the final product
    full_corr_map = np.zeros_like(h1)
    full_corr_map[:idx_limit, :] = corr_slice
    
    return full_corr_map
    

try:
    strain_full_h1, dt_h1 = get_gw150914_data(GW151226_h1_url,h1_GW151226_local)
    strain_full_l1, dt_l1 = get_gw150914_data(GW151226_l1_url,l1_GW151226_local)

    fs = 1.0 / dt_h1

    # Parameters for Spectrogram
    zoom_width, zoom_center = 2, 16.1
#    zoom_width, zoom_center = 1.2, 16.1-t_start


    t_start, t_end = zoom_center-2*zoom_width, zoom_center+2*zoom_width
    i_start, i_end = int(t_start * fs), int(t_end * fs)
    
    strain_h1 = strain_full_h1[i_start:i_end]
    strain_l1 = strain_full_l1[i_start:i_end]

    
    # to align the two chirps.
#    shift = 113
    
    # L1 was hit first, so we "delay" it to match H1
#    l1_aligned, h1_aligned = align_ligo_data(strain_l1,strain_h1)
    l1_aligned, h1_aligned = strain_l1, strain_h1

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
    # Parameters for Spectrogram
    NPERSEG = 2048  # The physical window size (your 'stiffness')
    NFFT = 2*NPERSEG    # The padded size (your 'interpolation')
    noverlap = int(NPERSEG*0.9)

    # 1. Prepare Window and calculate alpha
    win_seg = cmst_window(NPERSEG)
    alpha = calculate_alpha(win_seg)
    print(f"Calculated alpha: {alpha:.3f}")

    # 2. Compute Spectrogram
    f, t_spec, Sxx_h1_power = spectrogram(whitened_strain_h1, fs, 
                                       window=win_seg, 
                                   nperseg=NPERSEG, # Window is 256 long
                                   nfft=NFFT,       # Padded to 1024
                                   noverlap=noverlap,
                                   scaling='spectrum')

    f, t_spec, Sxx_l1_power = spectrogram(whitened_strain_l1, fs, 
                                   window=win_seg, 
                                   nperseg=NPERSEG, # Window is 256 long
                                   nfft=NFFT,       # Padded to 1024
                                   noverlap=noverlap,
                                   scaling='spectrum')


    # Convert Power to Amplitude (Linear Magnitude)
    Sxx_h1 = np.sqrt(Sxx_h1_power)
    Sxx_l1 = np.sqrt(Sxx_l1_power)


    Sxx_h1_sharp = sharpen_spectrogram(Sxx_h1,alpha*NPERSEG/NFFT)
    Sxx_l1_sharp = sharpen_spectrogram(Sxx_l1,alpha*NPERSEG/NFFT)
    eps = 1e-100
    Sxx_sharp = np.log(eps+np.sqrt(np.abs(Sxx_l1_sharp*Sxx_h1_sharp)))
    plot_grid(Sxx_sharp,70)
   

except Exception as e:
    print(f"Failed again: {e}")


