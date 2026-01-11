import numpy as np
import matplotlib.pyplot as plt

# --- PARAMETERS ---
N = 2048               # Resolution
RADIUS = 50.0         # Distance of planets
SUN_SIZE = 15          # Radius of Sun
PLANET_DB = -80       # Brightness of planets
NOISE_DB = -90       # Noise floor
PLANET_SIZE = 3 #Planet size

# ------------------
def generate_clean_image(size):
    """
    Generates a CLEAN Spatial Image (Large Sun + Planets, No Noise).
    """
    # 1. Start with Empty Space (Pure Zeros)
    image = np.zeros((size, size), dtype=complex)
    
    # Grid coordinates
    y, x = np.ogrid[:size, :size]
    cy, cx = size // 2, size // 2

    # 2. Draw The Sun (0 dB)
    mask_sun = ((y - cy)**2 + (x - cx)**2) <= SUN_SIZE**2
    image[mask_sun] = 1.0

    # 3. Add The Planets (-80 dB)
    p_amp = 10**(PLANET_DB / 20.0)
    angles = np.linspace(0, 2*np.pi, 8, endpoint=False)
    
    for angle in angles:
        py = int(cy + RADIUS * np.sin(angle))
        px = int(cx + RADIUS * np.cos(angle))
        for ang in np.linspace(0, 2*np.pi, 8, endpoint=False):
            # Calculate center of planet
        
            # Draw Planet Disk
            dist_sq = (y - py)**2 + (x - px)**2
            mask_planet = dist_sq <= PLANET_SIZE**2
            image[mask_planet] = p_amp
        
            
    return image

def show_log_image(image, title="Log-Scale Image"):
    plt.figure(figsize=(8, 8))
    
    # Convert to Magnitude
    mag = np.abs(image)
    
    # Convert to dB
    # We add a tiny epsilon (1e-15) so the pure zero background 
    # becomes approx -300 dB instead of -Infinity.
    data_db = 20 * np.log10(mag + 1e-15)
    
    # Normalize Peak to 0 dB
    data_db -= np.max(data_db)
    
    # Display
    # vmin=-110 keeps the dynamic range consistent with previous tests
    plt.imshow(data_db, cmap='inferno', vmin=-110, vmax=0)
    plt.colorbar(label='Amplitude (dB)', fraction=0.046, pad=0.04)
    plt.title(title, fontsize=14)
    plt.axis('off')
    plt.show()

def get_planck_window(N, epsilon=0.1):
    """
    The Planck-Taper Window (Standard in GW Analysis).
    Transitions from 0 to 1 using a sigmoid function 1/(1+exp(z)).
    """
    w = np.ones(N)
    # Define transition region width
    n_trans = int(epsilon * N)
    
    for i in range(1, n_trans):
        # x goes from 0 to 1 over the transition
        x = i / n_trans
        # The Planck taper function
        z = (1.0/x) + (1.0/(x-1.0))
        val = 1.0 / (1.0 + np.exp(z))
        
        # Apply to rising edge (left) and falling edge (right)
        w[i] = val
        w[N - 1 - i] = val
        
    # Zero out the exact edges
    w[0] = 0.0
    w[-1] = 0.0
    return w

def get_cmst_window(N):
    """The Analytical CMST / Gevrey Window."""
    t = np.linspace(-1, 1, N)
    t = t[1:-1] # Avoid singularity
    w = np.exp(t**2 + 1 / (t**2 - 1))
    w = np.pad(w, (1, 1), 'constant')
    return w
    
def process_fft(image, window):
    # 1. Apply Window
    windowed_image = image * window
    
    shifted_fft = np.fft.ifft2(np.fft.ifftshift(windowed_image))
    
    spec_db = 20 * np.log10(np.abs(shifted_fft) + 1e-15)
    spec_db -= np.max(spec_db)
    
    return spec_db

def show_data(data, title="Data Visualization"):

    # 2. Fix Complex: If data is complex, take Magnitude
    if np.iscomplexobj(data):
        print("Converting Complex Data to Magnitude...")
        data = np.abs(data)
    # 1. Fix Shape: If it's (2, 512, 512), pick the first one
    if hasattr(data, 'shape') and len(data.shape) == 3 and data.shape[0] == 2:
        print("Warning: Input was shape (2, N, N). Visualizing index [0].")
        data = data[0]

    
    # 3. Fix Log Scale: If data is linear amplitude, convert to dB
    # We assume if max value is small (< 100), it's likely linear. 
    # If it's already dB, max would be 0 or negative.
    if np.max(data) > 0: # Avoid log(0) issues for all-zero arrays
        # Add tiny epsilon to avoid log(0) errors
        data_db = 20 * np.log10(data + 1e-15)
        
        # Normalize so the Peak is 0 dB (Standard for FFT)
        data_db -= np.max(data_db)
        data = data_db

    # 4. Plot
    plt.figure(figsize=(8, 8))
    start = N // 2 -100
    end =  N // 2 + 100   # e.g., for N=512, this is 128 to 384

    # Apply to your axis
    plt.xlim(start, end)
    plt.ylim(end, start) # Note: 'end' first flips Y to match image coordinates
    # We clip the bottom at -120dB so the noise floor is visible but not distracting
    plt.imshow(data, cmap='inferno', vmin=-120, vmax=0)
    
    plt.colorbar(label='Amplitude (dB)')
    plt.title(title, fontsize=14)
    plt.axis("off")
    plt.show()


def show_comparison_plots(data1, data2, data3,top_title="Comparison Chart", titles=["1. Input", "2. Planck", "3. CMST"]):
    """
    Takes 3 sets of data (complex or real) plots them side-by-side zoomed to Center +/- 100 pixels.
    """
    data_list = [data1, data2, data3]
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    fig.suptitle(top_title, fontsize=20, y=0.98)    
    # Loop through each subplot and data source
    for i, (ax, data) in enumerate(zip(axs, data_list)):
        
        data = np.abs(data)
            
        # 2. Convert to dB
        # data_db = -20 * np.log10(data + 1e-15)
        data_db = -data
        # 3. Normalize Peak to 0 dB
        data_db -= np.max(data_db)
        
        # 4. Determine Zoom Bounds (Center +/- 100)
        N = data.shape[0]
        center = N // 2
        start = center - 100
        end = center + 100
        
        # 5. Plot
        # vmin=-110 ensures the noise floor is visible
        im = ax.imshow(data_db, cmap='inferno', vmin=-120, vmax=0)
        
        # 6. Apply Zoom
        ax.set_xlim(start, end)
        ax.set_ylim(end, start) # 'end' first puts larger number at bottom (standard image view)
        
        # Styling
        ax.set_title(titles[i], fontsize=12)
        ax.axis('off')
    
    # Add a shared colorbar
    cbar = fig.colorbar(im, ax=axs, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('Amplitude (dB)')
    plt.savefig("planets.png")
    print("Plot saved as planets.png")
    plt.show()


# --- EXECUTION ---
img_clean = generate_clean_image(N)
img_clean_abs = np.abs(img_clean)
img_clean_log = -20 * np.log10(img_clean_abs + 1e-15)

fft_result = np.fft.fftshift(np.fft.fft2(img_clean))

size = N

noise_amp = 10**(NOISE_DB / 20.0)

np.random.seed(42) 
noise = (np.random.normal(scale=noise_amp, size=(size, size)) + 1j * np.random.normal(scale=noise_amp, size=(size, size)))

fft_result_n = fft_result + noise

w_planck_1d = get_planck_window(N, epsilon=0.1) # 10% Taper
w_planck_2d = np.outer(w_planck_1d, w_planck_1d)

w_cmst_1d = get_cmst_window(N)
w_cmst_2d = np.outer(w_cmst_1d, w_cmst_1d)

planck = process_fft(fft_result_n,w_planck_2d)

cmst = process_fft(fft_result_n,w_cmst_2d)

show_comparison_plots(img_clean_log,planck,cmst,"Image->FFT->+ noise -> * window -> FFT\n Amplitude shown in dB", titles = ["Input","Planck", "CMST (p=2)"])

