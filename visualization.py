import matplotlib.pyplot as plt
import numpy as np


def plot_individual_modes(fit, selected_weighted_data, ax=None):
    """
    Plot the data, the total fit, and each Gaussian mode in log-diameter space.

    Parameters
    ----------
    fit : dict                 result from fitting.fit_psd (has u, y_norm, y_fit, modes).
    selected_weighted_data : str  column name, used for the title.
    ax : matplotlib Axes|None  draw into this axes (no window) or create a figure.
    """
    u = fit["u"]
    d = np.exp(u)
    show = ax is None
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    ax.plot(d, fit["y_norm"], "k-", lw=1.5, label="Data")
    ax.plot(d, fit["y_fit"], "r--", lw=1.3, label="Total fit")

    colors = ["#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]
    for i, m in enumerate(fit["modes"]):
        comp = m["A"] * np.exp(-0.5 * ((u - m["mu"]) / m["sigma"]) ** 2)
        c = colors[i % len(colors)]
        ax.fill_between(d, comp, alpha=0.25, color=c)
        ax.plot(d, comp, color=c, lw=1.0,
                label=f"Mode {i+1}: Dg={m['Dg']:.1f} µm ({m['weight']*100:.0f}%)")

    ax.set_xscale("log")
    ax.set_xlim(0.5, 100)
    ax.set_xlabel("Particle diameter [µm]")
    ax.set_ylabel("Normalized PSD (per unit ln d)")
    ax.set_title(f"Original vs fitted PSD — {selected_weighted_data}  "
                 f"(R²={fit['r_squared']:.3f})")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", ls=":", alpha=0.5)
    if show:
        plt.show()
