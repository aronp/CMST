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

    # 4. Helper: Calculate 99% Bandwidth
def calc_99_bw(signal):
    # Get Power Spectrum
    spec = np.abs(np.fft.rfft(signal))**2
        
    # Cumulative Sum
    cum_energy = np.cumsum(spec)
    total_energy = cum_energy[-1]
    norm_energy = cum_energy / total_energy
        
    # Find index crossing 0.99
    idx_99 = np.where(norm_energy >= 0.99)[0][0]
        
    return idx_99, norm_energy


def generate_pulse_demo():
    print("Generating Long-Chain Pulse Shaping Comparison...")
    
    # Setup a long random digital stream
    np.random.seed(42)
    N_bits = 1024  # Increased for better spectral resolution
    bits = np.random.randint(0, 2, N_bits)
    
    pulse_width = 256
    gap = 12            
    samples_per_symbol = pulse_width + gap
    
    # Define the Shapes
    t = np.linspace(-1, 1, pulse_width)
    
    # Shape A: Rectangular (Standard Digital)
    kernel_rect = np.ones_like(t)
    
    # Shape B: CMST(2) 
    # standard bump is best.
    k = 1
    with np.errstate(divide='ignore', invalid='ignore'):
        val = 1 + 1/(t**(2*k) - 1)
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

    idx_rect, _ = calc_99_bw(chain_rect)
    idx_cmst, _ = calc_99_bw(chain_cmst)
    
    print("Ratio of 99 BW ",idx_rect/idx_cmst)
    
    # Convert to Frequency Axis
    freqs = np.linspace(0, 0.5, len(np.fft.rfft(chain_rect)))
    bw_rect = freqs[idx_rect]
    bw_cmst = freqs[idx_cmst]

    
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Time Domain (Zoomed in to show pulse shape)
    zoom_subset = slice(0, 12000) # Show first few bits
    ax1.plot(chain_rect[zoom_subset], 'r', alpha=0.4, label='Rectangular Pulses (Standard)')
    ax1.plot(chain_cmst[zoom_subset], 'g', linewidth=2, label='CMST Pulses (Gevrey)')
    ax1.set_title(f"Time Domain: First few of {N_bits} bits")
    ax1.set_ylabel("Amplitude")
    ax1.legend(loc="upper right")
    ax1.grid(alpha=0.3)
    ax2.set_xlim(0,4000)
    
    # Freq Domain (Full Chain Analysis)
    ax2.plot(freqs, spec_rect, 'r', alpha=0.3, linewidth=0.5, label='Rect Spectrum (Algebraic Decay)')
    ax2.plot(freqs, spec_cmst, 'g', alpha=0.8, linewidth=0.8, label='CMST Spectrum (Gevrey Decay)')

    ax2.axvline(bw_cmst, color='green', linestyle='--', alpha=0.7)
    ax2.text(bw_cmst + 0.001, -155, f'CMST 99% BW\n{bw_cmst:.3f} Hz', color='darkgreen', fontweight='bold')
    
    # Mark the 99% BW for Rect
    ax2.axvline(bw_rect, color='red', linestyle='--', alpha=0.7)
    ax2.text(bw_rect + 0.001, -175, f'Rect 99% BW\n{bw_rect:.3f} Hz', color='darkred')
    
    ax2.set_title(f"Frequency Domain: Spectrum of {N_bits}-bit Chain")
    ax2.set_ylabel("Magnitude (dB)")
    ax2.set_xlabel("Normalized Frequency")
    ax2.set_ylim(-200, 0)
    ax2.legend(loc="upper right")
    ax2.grid(alpha=0.3)
    ax2.set_xlim(0,0.1)
    
    plt.tight_layout()
    plt.savefig('pulse_shaping_comparison.png')
    print("Done. Saved to 'pulse_shaping_comparison.png'")

    # Only show if NOT running on CI
    if not os.environ.get('CI'):
        plt.show()
    

if __name__ == "__main__":
    generate_pulse_demo()
