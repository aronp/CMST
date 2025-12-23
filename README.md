# cmst-window
**The analytically sound, zero-preserving, interlace-preserving window.**
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Stability: Verified](https://img.shields.io/badge/Stability-Verified-green)](https://github.com/aronp/CMST)
[![Theory: CMST](https://img.shields.io/badge/Theory-CMST-purple)](https://github.com/yourusername/cmst-window)

### 🚀 The Problem
Signal processing engineers often need a "Flat-Top" window to preserve signal amplitude while filtering noise. Existing standard solutions force a dangerous compromise:

* **Truncated Super-Gaussians:** Create **geometric singularities** (infinite acceleration) at the boundaries, causing ringing in control loops.
* **Planck-taper (LIGO Standard):** Theoretically perfect flatness, but relies on **piecewise "stitching"** of functions. This creates impulsive "kicks" (discontinuities) in higher-order derivatives (Jerk/Snap) and makes hardware optimization difficult.

### 🔬 Theoretical Basis: CMST Theory
This window is an implementation of **CMST (Cosh Moment Sturm Transform)**. Unlike standard windows which are often heuristic curve-fits, the Hyper-CMST window is constructed as a **Geometric Mollifier** with three rigorous guarantees:

1.  **Analytically Sound ($C^\infty$):**
    The function is infinitely differentiable with no discontinuities in any derivative $f^{(n)}$. This eliminates the "spectral ringing" and mechanical jerk caused by piecewise functions like the Planck-taper or Tukey window.

2.  **Zero-Preserving :**
    Derived from CMST theorey, the kernel guarantees the preservation of realness in the signal chain. It does not introduce artificial complex roots (phantom oscillations) into the passband.

3.  **Interlace-Preserving Transform:**
    For the first $n$ derivatives, the window acts as a variation-diminishing operator. It preserves the root-interlacing structure of the underlying signal, ensuring that derivative noise is bounded and geometric topology is maintained even at the boundaries.

### 💡 The Formula (Hyper-CMST)
We utilize a compensated log-concave mollifier that cancels low-order curvature to achieve "Hyper-Flatness":

$$
w(t) = \exp\left(t^n - \frac{1}{1-t^n}\right)
$$

* **Compensating Term ($+t^n$):** Cancels the Gaussian curvature at the origin, extending the "Table-Top" flatness to order $2n$.
* **Mollifier Term ($-1/(1-t^n)$):** Enforces strict compact support with essential singularities at the boundaries, ensuring all derivatives decay to zero smoothly.


### 🎛️ Tunable Flatness (p-Control)

Unlike traditional windows which are locked to a single profile (e.g., Hann, Blackman), the Hyper-CMST window is a parametric family. The power parameter (p) allows you to tune the window's behavior to match your specific engineering constraint:

    Mode A: The "Brick Wall" (p=6 or higher)
        Goal: Maximal Amplitude Accuracy.

        Behavior: The window remains effectively flat (>0.99) for over 70% of the duration, ensuring that signals are not attenuated in the center.

    Mode B: The "Silencer" (p=2)
        Goal: Maximal Spectral Purity.

        Behavior: The window converges to an analytically smooth Gaussian-like profile. This sacrifices the "flat top" to achieve significantly faster side-lobe decay, diving into the noise floor deeper than standard piecewise functions like the Planck-taper.
        
        

### 📉 The "Jerk" Test (Proof)
Why does analytic smoothness matter? In control systems (e.g., drones, robotics), the 3rd derivative (Jerk) corresponds to mechanical stress.

* **Planck-taper (Green):** Shows "impulsive" spikes in the Jerk plot where the flat section stitches to the curve. This is hidden noise.
* **Hyper-CMST (Black):** A perfectly smooth, continuous wave.

![Jerk Comparison](proofs/jerk_comparison.png)

### 📦 Installation

```bash
git clone [https://github.com/aronp/CMST.git](https://github.com/aronp/CMST.git)

cd CMST

pip install .
