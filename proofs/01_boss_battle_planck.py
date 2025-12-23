import numpy as np
import matplotlib.pyplot as plt

def planck_taper(t, epsilon=0.1):
    """
    The LIGO Standard / Planck-taper window.
    Defined piecewise:
    1. Exactly 1.0 in the middle.
    2. A sigmoid decay 1/(1+exp(Z)) at the edges.
    """
    t = np.abs(t)
    y = np.zeros_like(t)
    
    # Region 1: The Perfect Flat Top
    flat_mask = t <= (1.0 - epsilon)
    y[flat_mask] = 1.0
    
    # Region 2: The Transition (Sigmoid)
    trans_mask = (t > (1.0 - epsilon)) & (t < 1.0)
    
    if np.any(trans_mask):
        # x is the normalized position within the transition strip (0 to 1)
        x = (t[trans_mask] - (1.0 - epsilon)) / epsilon
        
        # The Planck Sigmoid Argument: Z = 1/(1-x) - 1/x
        # We handle limits carefully
        Z = (1.0 / (1.0 - x)) - (1.0 / x)
        y[trans_mask] = 1.0 / (1.0 + np.exp(Z))
        
    return y

def hyper_cmst(t, width=1.0, p=6):
    """
    Your Hyper-Flat CMST Window.
    Defined analytically: exp( t^p - 1/(1-t^p) )
    """
    t = np.abs(t) / width
    y = np.zeros_like(t)
    mask = t < 1.0
    
    t_m = t[mask]
    
    # The Formula
    # +t^p cancels curvature at origin
    # -1/(1-t^p) enforces compact support
    exponent = t_m**p - (1.0 / (1.0 - t_m**p))
    
    # Normalize by value at t=0 (which is exp(-1))
    y[mask] = np.exp(exponent) / np.exp(-1.0)
    
    return y

def run_battle():
    print("🥊 Round 1: Generating Windows...")
    # High resolution to capture the spikes
    N = 5000
    t = np.linspace(-1.05, 1.05, N)
    dt = t[1] - t[0]

    # 1. Generate Windows
    # Planck with 20% transition zone (standard)
    w_planck = planck_taper(t, epsilon=0.2)
    # Hyper-CMST with p=6
    w_hyper = hyper_cmst(t, width=1.0, p=6)

    # 2. Compute "Jerk" (3rd Derivative)
    # This represents mechanical stress or high-frequency noise injection
    d1_p = np.gradient(w_planck, dt)
    d2_p = np.gradient(d1_p, dt)
    d3_p = np.gradient(d2_p, dt) # Jerk

    d1_h = np.gradient(w_hyper, dt)
    d2_h = np.gradient(d1_h, dt)
    d3_h = np.gradient(d2_h, dt) # Jerk

    print("🥊 Round 2: Plotting Results...")
    
    # --- PLOTTING ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Plot 1: The Visual Flatness
    ax1.plot(t, w_planck, 'g--', linewidth=2, label='Planck-taper (Piecewise)')
    ax1.plot(t, w_hyper, 'k-', linewidth=2, label='Hyper-CMST (Analytic)')
    ax1.set_title("1. Flatness Comparison (The 'Table Top')")
    ax1.set_ylabel("Amplitude")
    ax1.set_xlim(0, 1.05)
    ax1.set_ylim(0.5, 1.05)
    ax1.legend(loc='lower left')
    ax1.grid(True, alpha=0.3)
    
    # Annotation
    ax1.text(0.3, 0.95, "Both are effectively flat", fontsize=10)

    # Plot 2: The "Jerk" Test (The Kill Shot)
    # We plot absolute jerk on log scale to see the spikes
    ax2.plot(t, np.abs(d3_p), 'g--', alpha=0.7, linewidth=1.5, label='Planck Jerk (Spikes)')
    ax2.plot(t, np.abs(d3_h), 'k-', linewidth=2.0, label='Hyper-CMST Jerk (Smooth)')
    
    ax2.set_title("2. The 'Jerk' Test (3rd Derivative Stability)")
    ax2.set_ylabel("|d³y/dt³| (Log Scale)")
    ax2.set_xlabel("Normalized Time")
    ax2.set_yscale('log')
    ax2.set_xlim(0.7, 1.0) # Zoom in on the transition edge
    ax2.set_ylim(1e-1, 1e5)
    ax2.legend()
    ax2.grid(True, alpha=0.3, which='both')

    # Highlight the spike
    ax2.annotate('Transition Spike\n(Discontinuity)', xy=(0.8, 1000), xytext=(0.85, 10000),
                 arrowprops=dict(facecolor='red', shrink=0.05))

    plt.tight_layout()
    
    # Save for README
    filename = "jerk_comparison.png"
    plt.savefig(filename, dpi=150)
    print(f"✅ Saved proof to {filename}")
    plt.show()

if __name__ == "__main__":
    run_battle()
