import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib.ticker import MaxNLocator

# 1. Define the exact integrals from your prompt
#    k=0 -> Base function
#    k=2 -> t^2 * exp(-t^4)
#    k=4 -> t^4 * exp(-t^4)
def integrand(t, s, k):
    return (t**k) * np.exp(-t**4) * np.cos(s * t)

def compute_integral(s, k):
    # Equivalent to NIntegrate[..., {t, 0, 10}]
    val, err = quad(integrand, 0, 10, args=(s, k))
    return val

# 2. Setup the Domain
s_values = np.linspace(0, 10, 500)

# 3. Compute and Normalize
#    We compute the "DC" component (s=0) first for normalization
norm_factors = {
    0: compute_integral(0, 0),
    2: compute_integral(0, 2),
    4: compute_integral(0, 4)
}

y0 = np.array([compute_integral(s, 0) for s in s_values]) / norm_factors[0]
y2 = np.array([compute_integral(s, 2) for s in s_values]) / norm_factors[2]
y4 = np.array([compute_integral(s, 4) for s in s_values]) / norm_factors[4]

# 4. Find Roots (for annotation)
def get_roots(x, y):
    spline = InterpolatedUnivariateSpline(x, y)
    roots = spline.roots()
    return roots[(roots >= 0) & (roots <= 10)]

r0 = get_roots(s_values, y0)
r2 = get_roots(s_values, y2)
r4 = get_roots(s_values, y4)

# 5. Plotting
fig, ax = plt.subplots(figsize=(10, 6))

# Plot the curves
ax.plot(s_values, y0, label=r'Base: $\int e^{-t^4} \cos(st) dt$', lw=2, color='navy')
ax.plot(s_values, y2, label=r'2nd Moment: $\int t^2 e^{-t^4} \cos(st) dt$', lw=2, linestyle='--', color='darkorange')
ax.plot(s_values, y4, label=r'4th Moment: $\int t^4 e^{-t^4} \cos(st) dt$', lw=2, linestyle='-.', color='green')

# Draw the zero line
ax.axhline(0, color='black', linewidth=0.8, alpha=0.3)

# Annotate the specific zeros to show Interlacing
# We just mark the first set of zeros to avoid clutter
ax.scatter(r4[0], 0, color='green', s=60, zorder=5, marker='s')
ax.scatter(r2[0], 0, color='darkorange', s=60, zorder=5, marker='o')
ax.scatter(r0[0], 0, color='navy', s=60, zorder=5, marker='D')

# Add text pointing out the shift
ax.annotate('Zeros shift inward\nwith higher moments', 
            xy=(r0[0], 0.1), xytext=(r0[0]+1, 0.4),
            arrowprops=dict(facecolor='black', arrowstyle='->'),
            fontsize=10)

# Styling
ax.set_title('Root Interlacing of CMST Transform Derivatives', fontsize=12)
ax.set_xlabel('Frequency $s$', fontsize=11)
ax.set_ylabel('Normalized Amplitude', fontsize=11)
ax.legend(loc='upper right')
ax.grid(True, alpha=0.2)
ax.set_xlim(0, 10)

# Halve the number of ticks (as requested previously)
ax.xaxis.set_major_locator(MaxNLocator(nbins=5))

plt.tight_layout()
plt.savefig('interlacing_chart.png', dpi=300)
plt.savefig('interlacing_chart.pdf', dpi=300)

plt.show()
