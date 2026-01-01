import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftshift

def standard_cmst(t, width=1.0, p=6):
    """
    The 'Old' Version.
    Formula: exp( -1 / (1-t^p) )
    Behaves like a standard 6th order Gaussian near zero.
    """
    t = np.abs(t) / width
    y = np.zeros_like(t)
    
    # Avoid singularity at boundary
    # We clip slightly for the comparison plot to avoid numerical warnings
    mask = t < 0.999999 
    
    t_m = t[mask]
    exponent = -1.0 / (1.0 - t_m**p)
    
    y[mask] = np.exp(exponent)
    
    # Normalize by peak (at t=0, val = exp(-1))
    return y / np.exp(-1.0)

def hyper_cmst(t, width=1.0, p=6):
    """
    The 'New' Version.
    Formula: exp( t^p - 1/(1-t^p) )
    The +t^p term cancels the 6th order curvature, creating 12th order flatness.
    """
    t = np.abs(t) / width
    y = np.zeros_like(t)
    mask = t < 0.999999
    
    t_m = t[mask]
    
    # The Compensated Formula
    exponent = t_m**p - (1.0 / (1.0 - t_m**p))
    
    y[mask] = np.exp(exponent)
    
    # Normalize by peak (at t=0, val = exp(-1))
    return y / np.exp(-1.0)

def run_civil_war():
    print("⚔️  Beginning Civil War: Standard vs Hyper...")
    
    # 1. Setup
    N = 4096
    t = np.linspace(-1.0, 1.0, N)
    
    # 2. Generate
    y_std = standard_cmst(t, p=6)
    y_hyper = hyper_cmst(t, p=6)
    
    # 3. Spectral Analysis
    # We use a simple FFT to check the side-lobe decay
    spec_std = 20 * np.log10(np.abs(fftshift(fft(y_std))) + 1e-20)
    spec_std -= np.max(spec_std)
    
    spec_hyper = 20 * np.log10(np.abs(fftshift(fft(y_hyper))) + 1e-20)
    spec_hyper -= np.max(spec_hyper)
    
    freqs = np.linspace(-0.5, 0.5, len(spec_std))

    # --- PLOTTING ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    
    # Plot 1: The Microscope (Zoomed Flatness)
    # We zoom in on the region where the window starts to roll off
    ax1.plot(t, y_std, 'b--', linewidth=2, label='Standard CMST (Dome-like)')
    ax1.plot(t, y_hyper, 'k-', linewidth=3, label='Hyper-CMST (Table-top)')
    
    ax1.set_title("1. The Flatness Upgrade (Zoomed)")
    ax1.set_ylabel("Amplitude")
    ax1.set_xlabel("Normalized Time")
    
    # CRITICAL ZOOM: 
    # This range highlights that Standard drops to 0.95 while Hyper is still 0.999
    ax1.set_xlim(0.0, 0.9) 
    ax1.set_ylim(0.8, 1.02)
    
    ax1.legend(loc='lower left')
    ax1.grid(True, alpha=0.3)
    
    # Annotation explaining the math
    ax1.text(0.2, 0.90, "Standard starts drooping\nat t=0.4", color='blue')
    ax1.text(0.5, 0.98, "Hyper stays flat\nuntil t=0.7", color='black', fontweight='bold')

    # Plot 2: Spectral Cost
    # Proving that we didn't ruin the frequency response to get the flatness
    ax2.plot(freqs, spec_std, 'b--', alpha=0.6, label='Standard Spectrum')
    ax2.plot(freqs, spec_hyper, 'k-', linewidth=1.5, label='Hyper Spectrum')
    
    ax2.set_title("2. The Spectral Cost (Frequency Domain)")
    ax2.set_ylabel("Magnitude (dB)")
    ax2.set_xlabel("Normalized Frequency")
    ax2.set_xlim(0, 0.1) # Zoom on main lobe and first few side lobes
    ax2.set_ylim(-150, 0)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("hyper_upgrade_proof.png", dpi=150)
    print("✅ Comparison saved to hyper_upgrade_proof.png")
    plt.show()

if __name__ == "__main__":
    run_civil_war()
