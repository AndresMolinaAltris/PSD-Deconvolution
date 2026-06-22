"""
Distribution preparation and mode seeding.

The fit now runs directly on the raw instrument bins (no cubic interpolation:
splines overshoot in the near-zero tail and invent spurious modes that the fit
then chases to absurd diameters). Peaks are detected in log-diameter space,
which is where the modes are symmetric and evenly sampled, and serve only as
*seeds* for the fit — the number of modes is set explicitly by the operator.
"""
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

DIAMETER_COL = "Particle diameter  [µm]"


def psd_arrays(df, param):
    """Return clean, sorted (diameter, differential) arrays from a loaded frame."""
    d = pd.to_numeric(df[DIAMETER_COL], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[param], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(d) & np.isfinite(y) & (d > 0)
    d, y = d[mask], y[mask]
    order = np.argsort(d)
    return d[order], y[order]


def detect_peaks(d, y, prominence_rel=0.02, distance=2):
    """
    Detect peak diameters in log-diameter space.

    Parameters
    ----------
    d, y : array         diameter and differential weight (raw bins).
    prominence_rel : float  required peak prominence as a fraction of the peak
                            height (relative, so it adapts per file). Lower it
                            to catch smaller modes; raise it to ignore bumps.
    distance : int       minimum separation between peaks, in bins.

    Returns the diameters at detected peaks (always at least one — the global
    maximum — so a seed is guaranteed).
    """
    if y.max() <= 0:
        return np.array([d[len(d) // 2]])
    y_rel = y / y.max()
    idx, _ = find_peaks(y_rel, prominence=prominence_rel, distance=distance)
    if len(idx) == 0:
        idx = np.array([int(np.argmax(y_rel))])
    return d[idx]
