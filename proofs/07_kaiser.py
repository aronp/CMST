import numpy as np
import scipy.signal.windows as windows
import matplotlib.pyplot as plt

def cmst_window(N, p=2):
    t = np.linspace(-1, 1, N)
    mask = np.abs(t) < 1.0 - 1e-9
    w = np.zeros(N)
    val = t[mask]
    t_p = val**p
    exponent = -(t_p * t_p) / (1.0 - t_p)
    w[mask] = np.exp(exponent)
    return w

def run_200db_comparison():
    N = 4096 # Higher N to see the peaks clearly
    fs = 10000
    f_carrier = 1000
    f_target = 1500 # Move target slightly further to ensure separation if lobes are wide
    
    # Signal: Carrier (0dB) + Target (-200dB = 1e-10)
    t = np.arange(N) / fs
    sig = 1.0 * np.sin(2*np.pi*f_carrier*t) + \
          1e-10 * np.sin(2*np.pi*f_target*t)
    
    # 1. CMST p=2
    w_cmst = cmst_window(N, p=2)
    w_cmst /= np.sum(w_cmst)
    
    # 2. Kaiser "Standard High" (Beta=14) - The one we used before
    w_kaiser_std = windows.kaiser(N, beta=14)
    w_kaiser_std /= np.sum(w_kaiser_std)
    
    # 3. Kaiser "Extreme" (Beta=24) - Tuned for >200dB
    w_kaiser_tuned = windows.kaiser(N, beta=24)
    w_kaiser_tuned /= np.sum(w_kaiser_tuned)
    
    # FFTs
    fft_cmst = np.abs(np.fft.rfft(sig * w_cmst))
    fft_k_std = np.abs(np.fft.rfft(sig * w_kaiser_std))
    fft_k_tuned = np.abs(np.fft.rfft(sig * w_kaiser_tuned))
    
    # dB
    db_cmst = 20 * np.log10(fft_cmst + 1e-30)
    db_cmst -= np.max(db_cmst)
    
    db_k_std = 20 * np.log10(fft_k_std + 1e-30)
    db_k_std -= np.max(db_k_std)
    
    db_k_tuned = 20 * np.log10(fft_k_tuned + 1e-30)
    db_k_tuned -= np.max(db_k_tuned)
    
    freqs = np.fft.rfftfreq(N, 1/fs)
    
    plt.figure(figsize=(12, 7))
    
    # Plot Standard Kaiser
    plt.plot(freqs, db_k_std, 'r:', linewidth=1.5, alpha=0.7, label=r'Kaiser ($\beta=14$) - Fails')
    
    # Plot Tuned Kaiser
    plt.plot(freqs, db_k_tuned, 'g--', linewidth=2, label=r'Kaiser ($\beta=24$) - Works but Wide')
    
    # Plot CMST
    plt.plot(freqs, db_cmst, 'b-', linewidth=2, label=r'CMST ($p=2$) - Sharp & Deep')
    
    plt.title("Extreme Dynamic Range Challenge (-200 dB Target)", fontsize=16)
    plt.ylabel("Magnitude (dB)", fontsize=14)
    plt.xlabel("Frequency (Hz)", fontsize=14)
    plt.legend(fontsize=12, loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.ylim(-260, 10)
    plt.xlim(0, 3000)
    
    # Annotate Target
    idx_target = np.argmin(np.abs(freqs - f_target))
    plt.annotate('Target (-200 dB)', xy=(freqs[idx_target], -200), xytext=(freqs[idx_target]+400, -190),
                 arrowprops=dict(facecolor='black', shrink=0.05), fontsize=12)
    
    plt.tight_layout()
    plt.savefig("cmst_vs_kaiser_200db.png", dpi=150)
    plt.show()

run_200db_comparison()
