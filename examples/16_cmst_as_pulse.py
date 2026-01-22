import os
import numpy as np
import matplotlib.pyplot as plt

def cmst_window(N):
    t = np.linspace(-1, 1, N)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = np.exp(t**4 / (t**2 - 1))
        w = np.where(np.abs(t) < 1, w, 0.0)
    return w

def get_db_spec(s):
    w = s * cmst_window(len(s))
    spec = np.abs(np.fft.rfft(w))
    spec /= np.max(spec) # Normalize
    return 20 * np.log10(spec + 1e-15)

def generate_pulse_demo():
    print("Generating Long-Chain Pulse Shaping Comparison...")
    
    # Setup a long random digital stream
    np.random.seed(42)
    N_bits = 1024  # Increased for better spectral resolution
    bits = np.random.randint(0, 2, N_bits)
    
    pulse_width = 128
    gap = 12            
    samples_per_symbol = pulse_width + gap
    
    # Define the Shapes
    t = np.linspace(-1, 1, pulse_width)
    
    # Shape A: Rectangular (Standard Digital)
    kernel_rect = np.ones_like(t)
    
    # Shape B: CMST(2) 
    # Formula: exp(1 + t^2 + 1/(t^2-1))
    k = 1
    with np.errstate(divide='ignore', invalid='ignore'):
        val = 1 + t**(2*k) + 1/(t**(2*k) - 1)
        kernel_cmst = np.exp(val)
    kernel_cmst[np.abs(t) >= 1] = 0
    
    # Generate the Signal Chain (Convolution)
    # Create impulse train
    sig_len = len(bits) * samples_per_symbol
    impulses = np.zeros(sig_len)
    impulses[::samples_per_symbol] = bits
    
    # Convolve
    print(f"Convolving {N_bits} bits...")
    chain_rect = np.convolve(impulses, kernel_rect, mode='same')
    chain_cmst = np.convolve(impulses, kernel_cmst, mode='same')
    
        
    freqs = np.linspace(0, 0.5, len(get_db_spec(chain_rect)))
    spec_rect = get_db_spec(chain_rect)
    spec_cmst = get_db_spec(chain_cmst)
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Time Domain (Zoomed in to show pulse shape)
    zoom_subset = slice(0, 4000) # Show first few bits
    ax1.plot(chain_rect[zoom_subset], 'r', alpha=0.4, label='Rectangular Pulses (Standard)')
    ax1.plot(chain_cmst[zoom_subset], 'g', linewidth=2, label='CMST Pulses (Gevrey)')
    ax1.set_title(f"Time Domain: First few of {N_bits} bits")
    ax1.set_ylabel("Amplitude")
    ax1.legend(loc="upper right")
    ax1.grid(alpha=0.3)
    
    # Freq Domain (Full Chain Analysis)
    ax2.plot(freqs, spec_rect, 'r', alpha=0.3, linewidth=0.5, label='Rect Spectrum (Algebraic Decay)')
    ax2.plot(freqs, spec_cmst, 'g', alpha=0.8, linewidth=0.8, label='CMST Spectrum (Gevrey Decay)')
    
    # Add an annotation arrow
    ax2.annotate('Spectral Floor (-120dB)', xy=(0.25, -120), xytext=(0.3, -100),
                 arrowprops=dict(facecolor='black', shrink=0.05))
    
    ax2.set_title(f"Frequency Domain: Spectrum of {N_bits}-bit Chain")
    ax2.set_ylabel("Magnitude (dB)")
    ax2.set_xlabel("Normalized Frequency")
    ax2.set_ylim(-200, 0)
    ax2.legend(loc="upper right")
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('pulse_shaping_comparison.png')
    print("Done. Saved to 'pulse_shaping_comparison.png'")

    # Only show if NOT running on CI
    if not os.environ.get('CI'):
        plt.show()
    

if __name__ == "__main__":
    generate_pulse_demo()
