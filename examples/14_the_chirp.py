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

def get_gw150914_data():
    print("Fetching valid download URL from GWOSC API...")
    api_url = "https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5"
    
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


try:
    strain, dt = get_gw150914_data()
    fs = 1.0 / dt
    N = len(strain)
    print(f"Data Loaded: {N} samples at {fs} Hz")

    win_cmst = cmst_window(N)

    # 4. FFT & Plot
    print("Computing FFTs...")
    spec_cmst = np.fft.rfft(strain * win_cmst)
    freqs = np.fft.rfftfreq(N, dt)

    fs = 1.0 / dt
    time_vector = np.arange(len(strain)) * dt

    print("Whitening the data (this extracts the chirp from the seismic noise)...")
    # We use the CMST window for the whitening step too
    whitened_strain = whiten_data(strain, dt, cmst_window)

    print("Generating Waterfall...")

    # Parameters for the Spectrogram
    # We focus on the event time. In this file, the event is at ~15.4 seconds.
    NFFT = 256        
    noverlap = 240    # High overlap for smooth image

    zoom_width = 0.07  
    zoom_center = 15.39

    # Create the figure
    plt.figure(figsize=(10, 8))

    # Compute Spectrogram using YOUR window
    f, t_spec, Sxx = spectrogram(whitened_strain, fs, 
                             window=cmst_window(NFFT), 
                             nperseg=NFFT, 
                             noverlap=noverlap)

    # Plot
    plt.pcolormesh(t_spec, f, Sxx, shading='gouraud', cmap='viridis')

    # Zoom in on the Chirp
    plt.xlim(zoom_center - zoom_width, zoom_center + zoom_width)
    plt.ylim(30, 300) # The "Audible" range of the black holes

    plt.title(f'LIGO GW150914: Spectrogram Analysis (CMST (p=2) Window)')
    plt.ylabel('Frequency (Hz)')
    plt.xlabel('Time (s)')
    plt.colorbar(label='Signal-to-Noise Ratio')
    plt.grid(False) # Turn off grid to see the track clearly
    plt.savefig("Chirp.png")
    print("Saved Chirp.png")
    
    # Only show if NOT running on CI
    if not os.environ.get('CI'):
        plt.show()


except Exception as e:
    print(f"Failed again: {e}")
