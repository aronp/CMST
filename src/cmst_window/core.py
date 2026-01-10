import numpy as np


def cmst(N, p=2, sym=True):
    """
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
    

