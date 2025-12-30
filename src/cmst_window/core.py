import numpy as np

def cmst(N, p=2, sym=True):
    """
    Matches scipy.signal.get_window conventions.
    """
    # Generate the time grid [-1, 1]
    if sym:
        # Symmetric window (for filter design)
        t = np.linspace(-1, 1, N)
    else:
        # Periodic window (for spectral analysis/FFT)
        t = np.linspace(-1, 1, N, endpoint=False)
    
    return cmst_window(t, width=1.0, power=p)

def cmst_window(t, width=1.0, power=2):
    """
    Generates the Hyper-CMST (Compensated Mollifier) window function.
    
    This window uses a compensated log-concave mollifier to achieve maximal 
    flatness (equivalent to order 2*power) while maintaining C-infinity 
    analytic smoothness at the boundaries.
    
    Unlike piecewise functions (Planck-taper, Tukey), this function is 
    analytically defined everywhere, ensuring that high-order derivatives 
    (Jerk, Snap) remain continuous and bounded.
    
    Formula: 
        y(t) ~ exp( t^p - 1/(1-t^p) )  for |t| < 1
    
    Args:
        t (array-like): Input time or frequency coordinates.
        width (float): The half-width of the window. The function is strictly 
                       zero for |t| >= width.
        power (int): The base order of the mollifier. Must be an even integer. 
                     Default is 6 (yielding 12th-order flatness).
                     
    Returns:
        numpy.ndarray: The window values, normalized to a peak amplitude of 1.0.
    """
    # 1. Input Validation
    t = np.asarray(t, dtype=float)
    
    if width <= 0:
        raise ValueError(f"Window width must be positive. Received {width}.")
        
    if power % 2 != 0:
        raise ValueError(f"Power must be an even integer to ensure symmetry. Received {power}.")
    
    # 2. Initialization
    y = np.zeros_like(t)
    
    # 3. Normalization
    # Scale t so the boundaries are at -1 and +1
    x = t / width
    
    # 4. Strict Compact Support
    # We only compute the formula strictly inside the domain (-1, 1).
    # At exactly |x| = 1, the limit is 0, which is handled by the initialization.
    mask = np.abs(x) < 1.0
    
    if not np.any(mask):
        return y
        
    # 5. The Hyper-CMST Formula
    # Extract valid points to avoid DivisionByZero in the denominator
    x_valid = x[mask]
    
    # Compute x^p once
    x_p = x_valid ** power
    
    # The exponent has two terms:
    # Term A: +x^p (The Compensator) -> Cancels curvature at origin
    # Term B: -1/(1-x^p) (The Mollifier) -> Enforces zero at boundary
    # 1.0 Normalises the function.
    exponent = 1.0 + x_p - (1.0 / (1.0 - x_p))
    
    # 6. Apply
    y[mask] = np.exp(exponent) 
    
    return y
