import numpy as np
import matplotlib.pyplot as plt
import h5py
import requests
import io
import json
import os
from scipy.signal import spectrogram
from scipy.interpolate import interp1d

h1_api_url = "https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5"
l1_api_url = "https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3/L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"

h1_local = "H1_data.hdf5"
l1_local = "L1_data.hdf5"


def cmst_window(N):
    t = np.linspace(-1, 1, N)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = np.exp(t**4 / (t**2 - 1))
    w = np.where(np.abs(t) < 1, w, 0.0)
    return w

def planck_taper(N, epsilon=0.1):
    t = np.linspace(0, 1, N)
    w = np.ones(N)
    # Left taper
    n_left = int(epsilon * (N - 1))
    for i in range(1, n_left):
        w[i] = 1.0 / (np.exp(epsilon * (N - 1) / i - epsilon * (N - 1) / (n_left - i)) + 1)
    # Right taper
    w[0] = 0
    w[-1] = 0
    w[N-n_left:-1] = w[1:n_left][::-1]
    return w

def calculate_alpha(window):
    """
    Calculates the sharpening parameter alpha directly 
    from the spectral footprint of the discrete window.
    """
    sum_w = np.sum(window)
    t = np.linspace(-1, 1, len(window))
    t2 = t*t
    t2sum = np.sum(t2*window)/sum_w
    return t2sum

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

def whiten_data(strain, dt, window_func):
    N = len(strain)
    spec = np.fft.rfft(strain * window_func(N))
    freqs = np.fft.rfftfreq(N, dt)
    mag = np.abs(spec)
    psd_smooth = interp1d(freqs, mag, kind='nearest', fill_value="extrapolate")(freqs)
    whitened_spec = spec / (psd_smooth + 1e-20)
    return np.fft.irfft(whitened_spec, n=N)

# --- Execution ---
try:
    strain_full, dt = get_gw150914_data(h1_api_url,h1_local)
    fs = 1.0 / dt
    
    # Slice around the event (approx 15.3 - 15.5s)
    t_start, t_end = 15.0, 16.0 
    i_start, i_end = int(t_start * fs), int(t_end * fs)
    
    # Chop it
    strain = strain_full[i_start:i_end]
    
    whitened_strain = whiten_data(strain, dt, cmst_window)

    # Parameters for Spectrogram
    zoom_width, zoom_center = 0.07, 15.39-t_start
    # Parameters for Spectrogram
    NPERSEG = 256  # The physical window size (your 'stiffness')
    NFFT = 1024    # The padded size (your 'interpolation')
    noverlap = 240

    # 1. Prepare Window and calculate alpha
    win_seg = cmst_window(NPERSEG)
    alpha = calculate_alpha(win_seg)
    print(f"Calculated alpha: {alpha:.3f}")

    # 2. Compute Spectrogram
#    f, t_spec, Sxx = spectrogram(whitened_strain, fs, 
#                                 window=win_seg, 
#                                 nperseg=NFFT, 
#                                 noverlap=noverlap)
    # Compute Spectrogram (defaults to Power)
#    f, t_spec, Sxx_power = spectrogram(whitened_strain, fs, 
#                                   window=win_seg, 
#                                   nperseg=NFFT, 
#                                   noverlap=noverlap,
#                                   scaling='spectrum') # 'spectrum' ensures units are V^2

    f, t_spec, Sxx_power = spectrogram(whitened_strain, fs, 
                                   window=win_seg, 
                                   nperseg=NPERSEG, # Window is 256 long
                                   nfft=NFFT,       # Padded to 1024
                                   noverlap=noverlap,
                                   scaling='spectrum')
                                   
    # Convert Power to Amplitude (Linear Magnitude)
    Sxx = np.sqrt(Sxx_power)


    # 3. Frequency-Domain Laplacian Sharpening
    # S_sharp = S - alpha * Laplacian_f(S)
    laplacian_f = np.zeros_like(Sxx)
    # Approximate 2nd derivative along frequency axis (axis 0)
    laplacian_f[1:-1, :] = Sxx[2:, :] - 2*Sxx[1:-1, :] + Sxx[0:-2, :]
    Sxx_sharp = Sxx - alpha * laplacian_f
    
    # Clip to prevent negative values (non-physical energy)
    Sxx_sharp = np.maximum(Sxx_sharp, 1e-200)
    Sxx_sharp = Sxx - (alpha / (NFFT/NPERSEG)) * laplacian_f
    # Plot
    plt.figure(figsize=(14, 8))
    
    plt.pcolormesh(t_spec, f, Sxx_sharp, shading='gouraud', cmap='viridis')
    contour_filled = plt.contourf(t_spec, f, Sxx_sharp, levels=100, cmap='inferno')
    contour_lines = plt.contour(t_spec, f, Sxx_sharp, levels=60, colors='white', linewidths=0.5, alpha=0.5)

    # Zoom in on the Chirp
    plt.xlim(zoom_center - 1.25 * zoom_width, zoom_center + 2/2 * zoom_width)
    plt.ylim(0, 500) # The "Audible" range of the black holes

    plt.title(f'LIGO GW150914: Spectrogram Analysis (CMST (p=2) Window)')
    plt.ylabel('Frequency (Hz)')
    plt.xlabel('Time (s)')
    plt.grid(False) # Turn off grid to see the track clearly
    plt.colorbar(contour_filled, label='Signal-to-Noise Ratio')
    plt.grid(True)
    
    plt.savefig("Chirp.png")
    print("Saved Chirp.png")
    
    # Only show if NOT running on CI
    if not os.environ.get('CI'):
        plt.show()

except Exception as e:
    print(f"Failed again: {e}")
