import numpy as np

# CMST Stability Constant (Psi) to 17 sig figs.
PSI = -0.079056818131899115

def cmst(N, p=2, sym=True):
    r"""
    Return a Cosh-Moment Sturm Transform (CMST) window.

    The CMST window is a "perfectly smooth" (C-infinity) test function generated 
    from a compensated log-concave mollifier. Unlike classical windows (Blackman, 
    Kaiser) which are analytic but finite-order, or piece-wise windows (Tukey, 
    Planck) which have discontinuous derivatives, the CMST window belongs to the 
    Gevrey class of regularity. It possesses infinite derivatives at the boundary 
    while maintaining a super-algebraic decay rate.

    This window is ideal for high-precision spectral analysis (e.g., gravitational 
    wave detection) where leakage from "distant" noise sources must be suppressed 
    below the machine precision floor.

    Parameters
    ----------
    N : int
        Number of points in the output window. If zero or less, an empty
        array is returned.
    p : int, optional
        The "flatness" order of the window. Must be an even integer. 
        Higher values of `p` create a flatter top (better amplitude accuracy) 
        and sharper transition, but increase the near-sidelobe levels.
        - p=2: Maximum dynamic range (deepest floor), widest main lobe.
        - p=6: Balanced resolution and dynamic range.
        Default is 2.
    sym : bool, optional
        When True (default), generates a symmetric window, for use in filter
        design.
        When False, generates a periodic window, for use in spectral analysis.

    Returns
    -------
    w : ndarray
        The window, with the maximum value normalized to 1 (though the value 1
        does not appear if `N` is even and `sym` is True).

    Notes
    -----
    The CMST window is defined mathematically is related to the standard "bump function" 
    from distribution theory, adapted for signal processing:
    
    .. math:: \Phi(t) = \exp\left(\frac{1}{t^p - 1}\right) \quad \text{for } |t| < 1

    References
    ----------
    .. [1] "On a Class of Zero-Preserving Bilateral Transforms and their Application to Precision Windowing", https://github.com/aronp/CMST/cmst.pdf

Examples
    --------
    Plot the window and its frequency response to demonstrate the super-algebraic 
    decay (no noise floor):

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.fft import fft, fftshift
    >>> # Assuming cmst is imported or defined
    
    >>> # 1. Generate Window
    >>> window = cmst(100, p=2)
    
    >>> # 2. Calculate Frequency Response (FFT)
    >>> # We zero-pad to 2048 to see the smooth sidelobe structure
    >>> A = fft(window, 2048)
    >>> mag = np.abs(fftshift(A))
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    
    >>> # Add a tiny epsilon (1e-35) to avoid log(0) warnings because
    >>> # the CMST decay can hit the machine precision floor.
    >>> response = 20 * np.log10(mag / np.max(mag) + 1e-35)
    
    >>> # 3. Plot
    >>> fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    
    >>> # Time Domain
    >>> ax1.plot(window)
    >>> ax1.set_title("CMST Window (Time Domain)")
    >>> ax1.set_xlabel("Sample")
    >>> ax1.set_ylabel("Amplitude")
    
    >>> # Frequency Domain
    >>> ax2.plot(freq, response)
    >>> ax2.set_title("Frequency Response (dB)")
    >>> ax2.set_ylabel("Magnitude [dB]")
    >>> ax2.set_xlabel("Normalized Frequency")
    >>> ax2.set_ylim(-200, 10)  # Show deep dynamic range
    >>> ax2.grid(True)
    
    >>> plt.tight_layout()
    >>> plt.show()           
    """
    
    if N <= 0:
        return np.array([])
    
    # 1. Generate the time grid [-1, 1]
    if sym:
        # Symmetric window (filter design: start=-1, end=1)
        # Note: If N is even, the exact center 0.0 is straddled, which is correct.
        t = np.linspace(-1, 1, N)
    else:
        # Periodic window (DFT/FFT: start=-1, end=1 - 1/N)
        # We exclude the last point because the DFT assumes the signal repeats.
        t = np.linspace(-1, 1, N, endpoint=False)
    
    return cmst_window(t, width=1.0, power=p)


def cmst_window(t, width=1.0, power=2):
    """
    Compute the Hyper-CMST (Compensated Mollifier) window function value.
    
    This is the continuous kernel function. It can be used directly for 
    non-uniform time grids or theoretical analysis.

    Parameters
    ----------
    t : array_like
        Input coordinates (time or frequency).
    width : float, optional
        The half-width of the support. The function is strictly zero for 
        abs(t) >= width. Default is 1.0.
    power : int, optional
        The order of the denominator polynomial. Must be even. Default is 2.

    Returns
    -------
    y : ndarray
        The window values.
    """
    # 1. Input Validation
    t = np.asarray(t, dtype=float)
    
    if width <= 0:
        raise ValueError(f"Window width must be positive. Received {width}.")
        
    if power % 2 != 0:
        raise ValueError(f"Power must be an even integer. Received {power}.")
    
    # 2. Initialization
    y = np.zeros_like(t)
    
    # 3. Normalization
    x = t / width
    
    # 4. Strict Compact Support Mask
    mask = np.abs(x) < 1.0
    
    if not np.any(mask):
        return y
        
    # 5. Robust Computation
    x_valid = x[mask]
    
    # Calculate power first
    x_p = x_valid ** power
    
    # Numerical Stability Check:
    # Floating point rounding can cause x^p to equal 1.0 even if x < 1.0.
    # We must exclude these points to prevent DivisionByZero.
    safe_indices = x_p < 1.0
    x_p_safe = x_p[safe_indices]
    
    # 6. The Hyper-CMST Formula
    exponent = 1.0 + x_p_safe - (1.0 / (1.0 - x_p_safe))
    
    # 7. Map back to full array
    # Create a temp buffer for the masked region
    y_masked = np.zeros_like(x_valid)
    # Fill the safe calculated values (underflow is handled natively by exp)
    y_masked[safe_indices] = np.exp(exponent)
    
    # Fill the original array
    y[mask] = y_masked
    
    return y
    
def cmst_pure(N, p=2, sym=True):
    if N <= 0:
        return np.array([])
    
    # 1. Generate the time grid [-1, 1]
    if sym:
        t = np.linspace(-1, 1, N)
    else:
        # Periodic window (DFT/FFT: start=-1, end=1 - 1/N)
        # We exclude the last point because the DFT assumes the signal repeats.
        t = np.linspace(-1, 1, N, endpoint=False)
    
    return cmst_pure_window(t, width=1.0)
    
def cmst_pure_window(t, width=1.0):
    """
    Compute the CMST window function value.
    
    This version is the pure standard bump.  It has the best side lobe decay properties, but cmst(p=2) does nearly as well and keeps more of the data.

    Parameters
    ----------
    t : array_like
        Input coordinates (time or frequency).
    width : float, optional
        The half-width of the support. The function is strictly zero for 
        abs(t) >= width. Default is 1.0.

    Returns
    -------
    y : ndarray
        The window values.
    """
    t = np.asarray(t, dtype=float)
    
    if width <= 0:
        raise ValueError(f"Window width must be positive. Received {width}.")
        
    y = np.zeros_like(t)
    
    x = t / width
    
    mask = np.abs(x) < 1.0
    
    if not np.any(mask):
        return y
        
    x_valid = x[mask]
    
    x_p = x_valid ** 2
    
    # Numerical Stability Check:
    # Floating point rounding can cause x^p to equal 1.0 even if x < 1.0.
    # We must exclude these points to prevent DivisionByZero.
    safe_indices = x_p < 1.0
    x_p_safe = x_p[safe_indices]
    
    exponent = 1.0 - (1.0 / (1.0 - x_p_safe))
    
    y_masked = np.zeros_like(x_valid)
    # Fill the safe calculated values (underflow is handled natively by exp)
    y_masked[safe_indices] = np.exp(exponent)
    
    # Fill the original array
    y[mask] = y_masked
    
    return y
    

def generate_cmst_fir(taps, fs, bandwidth, freq_power=4):
    """Generates linear FIR coefficients from the continuous CMST frequency mask."""
    if taps % 2 == 0:
        taps += 1

    freqs = np.fft.rfftfreq(taps, 1 / fs)
    u_50 = 0.5 * ( + np.sqrt(4.0*np.log(2) + np.log(2) ** 2) + -np.log(2))
    x_50 = u_50 ** (1.0 / freq_power)
    actual_width = bandwidth / x_50

    x = freqs / actual_width
    mask = np.zeros_like(x)
    valid = np.abs(x) < 1.0

    with np.errstate(divide='ignore', invalid='ignore'):
        x_valid = np.abs(x[valid]) ** freq_power
        exponent = 1 + x_valid - (1.0 / (1.0 - x_valid))
        mask[valid] = np.exp(exponent)

    impulse_response = np.fft.irfft(mask, n=taps)
    fir_coeffs = np.fft.fftshift(impulse_response)
    res_fir = fir_coeffs * cmst(taps, p=freq_power)

    return res_fir / np.sum(res_fir)


def cmst_lowpass(strain, fs, f_high, freq_power=4):
    """
    Ultra-fast, zero-phase CMST low-pass filter.
    Centers the window at 0 Hz to perfectly pass DC/low frequencies,
    with a mathematically exact 50% amplitude roll-off at f_high.
    """
    N = len(strain)
    if N < 8:
        return strain

    # 1. Global FFT (Direct, no padding)
    spec = np.fft.rfft(strain)
    freqs = np.fft.rfftfreq(N, 1 / fs)

    u_50 = 0.5 * ( + np.sqrt(4.0*np.log(2) + np.log(2) ** 2) + -np.log(2))

    # Cancel out the x^p exponent scaling
    x_50 = u_50 ** (1.0 / freq_power)

    # --- LOW-PASS GEOMETRY ---
    # Center the peak of the window exactly at 0 Hz
    f_center = 0.0

    # The mathematical distance from 0 Hz to the 50% mark is simply f_high
    actual_width = f_high / x_50
    # -------------------------

    # 3. Apply the Zero-Phase Frequency Mask
    freq_weight = cmst_window(freqs - f_center, width=actual_width, power=freq_power)
    filtered_spec = spec * freq_weight

    # 4. IFFT directly back to the time domain
    filtered = np.fft.irfft(filtered_spec, n=N)

    return filtered


def apply_spectral_inversion(lowpass_taps):
    taps = np.array(lowpass_taps, dtype=float)
    num_taps = len(taps)

    # Ensure an odd number of taps for a true center element
    if num_taps % 2 == 0:
        raise ValueError("Spectral inversion requires an odd number of taps.")

    # Find the central index
    nc = (num_taps - 1) // 2

    # Calculate the actual DC gain (sum of all unnormalized taps)
    dc_gain = np.sum(taps)

    # 1. Invert the sign of every tap
    highpass_taps = -taps

    # 2. Add the original DC gain to the central tap
    # This guarantees the new sum of all taps equals exactly 0.0
    highpass_taps[nc] += dc_gain

    return highpass_taps

def apply_spectral_modulation(lowpass_taps):
    taps = np.array(lowpass_taps, dtype=float)
    num_taps = len(taps)

    # Generate the [1, -1, 1, -1...] array
    alternating_sequence = np.power(-1, np.arange(num_taps))

    # Multiply element-wise
    highpass_taps = taps * alternating_sequence

    return highpass_taps

def generate_cmst_lp_fir(taps, fs, bandwidth, freq_power=4):
    return generate_cmst_fir(taps, fs, bandwidth, freq_power)

def generate_cmst_hp_fir(taps, fs, bandwidth, freq_power=4):
    return apply_spectral_inversion(generate_cmst_fir(taps, fs, bandwidth, freq_power))

def make_ola_50(window):
    """
    Converts a generic window array into a 50% overlap-add normalized window.
    Assumes the input array length is even.
    """
    window = np.asarray(window, dtype=float)
    N = len(window)
    
    if N % 2 != 0:
        raise ValueError("Window length N must be an even number for a clean 50% overlap.")
        
    hop = N // 2
    denom = np.zeros_like(window)
    
    # 1. Add the current window
    denom += window
    
    # 2. Add the overlap from the previous window (shifted right into the first half)
    denom[:hop] += window[hop:]
    
    # 3. Add the overlap from the next window (shifted left into the second half)
    denom[hop:] += window[:hop]
    
    # 4. Apply a deep safety floor to prevent division by zero in extreme decay tails
    denom[denom < 1e-100] = 1.0
    
    return window / denom




