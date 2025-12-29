import numpy as np
import matplotlib.pyplot as plt

# --- 1. Define CMST Implementation ---
def win_cmst(N, p):
    # Domain t from -1 to 1.
    # We use a tiny offset (1e-9) to avoid division by zero at the exact boundaries.
    t = np.linspace(-1 + 1e-9, 1 - 1e-9, N)
    w = np.exp(1 + t**p - 1/(1-t**p))
    return w

# --- 2. Manual Planck Window Implementation ---
def planck_window(N, epsilon=0.1):
    """
    Generates a Planck-taper window manually.
    Transitions occur over 'epsilon * N' samples at both ends.
    """
    w = np.ones(N)
    n_trans = int(epsilon * N)
    
    # Rising Edge
    k = np.arange(1, n_trans) 
    Z_rise = (epsilon * N) * (1/k + 1/(k - epsilon * N))
    w[1:n_trans] = 1 / (np.exp(Z_rise) + 1)
    w[0] = 0 
    
    # Falling Edge (Symmetric)
    w[N - n_trans : N - 1] = w[1:n_trans][::-1]
    w[N-1] = 0 
    
    return w

# --- 3. Configuration ---
N = 2048
cmst_orders = [2, 4, 6]
cmst_colors = ['#2196F3', '#FFC107', '#4CAF50'] # Blue, Amber, Green
planck_eps = 0.1
planck_color = 'black'

# --- 4. Plotting ---
plt.figure(figsize=(8, 6))

# Plot CMST variations
for p, col in zip(cmst_orders, cmst_colors):
    w = win_cmst(N, p)
    plt.plot(w, linewidth=2.5, label=f'CMST (p={p})', color=col, alpha=0.9)

# Plot Planck benchmark
w_pl = planck_window(N, epsilon=planck_eps)
plt.plot(w_pl, linewidth=2.5, linestyle='--', color=planck_color, label=f'Planck (ε={planck_eps})')

# Formatting
plt.title("Time Domain: Tunable 'Flat Top'", fontsize=14, fontweight='bold')
plt.ylabel("Amplitude (Normalized)", fontsize=12)
plt.xlabel("Sample Index", fontsize=12)
plt.grid(True, which='major', linestyle='-', alpha=0.4)
plt.legend(loc='lower center', frameon=True, fancybox=True, framealpha=0.9)
plt.ylim(0, 1.05)
plt.xlim(0, N)

# Add annotation
plt.annotate('Flatness increases with p', xy=(N/2, 0.95), xytext=(N/2, 0.5),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5),
            ha='center', fontsize=11)

plt.tight_layout()
plt.savefig('tunable_flatness_time_only.png', dpi=300)
plt.savefig('tunable_flatness_time_only.pdf', dpi=300)

plt.show()
