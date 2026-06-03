import numpy as np
from scipy.signal import find_peaks
import pandas as pd
from scipy.interpolate import interp1d

def detect_peaks(data_df, fitting_param,
                 prominence=0.001, distance=5,
                 height=None, width=None, threshold=None):
    """
    Detect peaks in the particle size distribution to determine modality.

    Parameters:
    - data_df: Preprocessed DataFrame.
    - fitting_param: Column name for the fitting parameter.
    - prominence: Required prominence of peaks (passed to scipy.signal.find_peaks).
      Lower it to capture small/low-amplitude peaks (default: 0.001).
    - distance: Required minimal horizontal distance (in samples) between peaks.
      Lower it to resolve peaks that sit close together (default: 5).
    - height: Optional required peak height (default: None = no height filter).
    - width: Optional required peak width in samples (default: None = no width filter).
    - threshold: Optional required vertical distance to neighbouring samples
      (default: None = no threshold filter).

    Returns:
    - diameters: Array of particle diameters.
    - diff_param: Normalized differential distribution.
    - num_modes: Number of detected peaks (modes).
    - peak_sizes: Diameters at peak locations.
    """
    data_df['Particle diameter  [µm]'] = pd.to_numeric(data_df['Particle diameter  [µm]'], errors='coerce')
    diameters = data_df['Particle diameter  [µm]'].values
    diff_param = data_df[fitting_param].values
    diff_param = diff_param / np.trapezoid(diff_param, diameters)
    peaks, _ = find_peaks(diff_param, prominence=prominence, distance=distance,
                          height=height, width=width, threshold=threshold)
    num_modes = len(peaks)
    peak_sizes = diameters[peaks]
    return diameters, diff_param, num_modes, peak_sizes

def prepare_distribution_data(data_df, fitting_param):
    """
    Prepare and normalize the particle size distribution data.
    
    Parameters:
    - data_df: Preprocessed DataFrame.
    - fitting_param: Column name for the fitting parameter.
    
    Returns:
    - diameters: Array of particle diameters.
    - diff_param: Normalized differential distribution.
    """
    data_df['Particle diameter  [µm]'] = pd.to_numeric(data_df['Particle diameter  [µm]'], errors='coerce')
    diameters = data_df['Particle diameter  [µm]'].values
    diff_param = data_df[fitting_param].values
    diff_param = diff_param / np.trapezoid(diff_param, diameters)
    return diameters, diff_param

def prepare_distribution_data_interpolation(data_df, fitting_param, interp_points=1000):
    """
    Prepare and normalize the particle size distribution data with interpolation.
    
    Parameters:
    - data_df: Preprocessed DataFrame.
    - fitting_param: Column name for the fitting parameter.
    - interp_points: Number of points for the interpolated grid (default: 1000).
    
    Returns:
    - diameters: Array of interpolated particle diameters.
    - diff_param: Normalized interpolated differential distribution.
    """
    data_df['Particle diameter  [µm]'] = pd.to_numeric(data_df['Particle diameter  [µm]'], errors='coerce')
    original_diameters = data_df['Particle diameter  [µm]'].values
    original_diff_param = data_df[fitting_param].values
    
    # Sort by diameters to ensure ascending order
    sort_idx = np.argsort(original_diameters)
    original_diameters = original_diameters[sort_idx]
    original_diff_param = original_diff_param[sort_idx]
    
    # Create a log-spaced interpolated grid to match typical particle size scaling
    if len(original_diameters) < 2:
        raise ValueError("Insufficient data points for interpolation.")
    min_d = np.min(original_diameters)
    max_d = np.max(original_diameters)
    if min_d <= 0:
        raise ValueError("Diameters must be positive for log-spacing.")
    diameters = np.logspace(np.log10(min_d), np.log10(max_d), interp_points)
    
    # Interpolate using linear interpolation to avoid oscillations/artifacts
    #interpolator = interp1d(original_diameters, original_diff_param, kind='linear', fill_value=0, bounds_error=False)
    #diff_param = interpolator(diameters)
    
    # Interpolate using cubic interpolation for smoother results (change to 'linear' if preferred)
    interpolator = interp1d(original_diameters, original_diff_param, kind='cubic', fill_value=0, bounds_error=False)
    diff_param = interpolator(diameters)

    # Ensure non-negative values
    diff_param = np.maximum(diff_param, 0)
    
    # Normalize the interpolated distribution
    integral = np.trapezoid(diff_param, diameters)
    if integral == 0:
        raise ValueError("Interpolated distribution integrates to zero.")
    diff_param = diff_param / integral


    # Normalize original for comparison
    original_integral = np.trapezoid(original_diff_param, original_diameters)
    if original_integral == 0:
        raise ValueError("Original distribution integrates to zero.")
    original_normalized = original_diff_param / original_integral
    
    return diameters, diff_param

def find_distribution_peaks(diameters, diff_param,
                            prominence=0.001, distance=20,
                            height=None, width=None, threshold=None):
    """
    Detect peaks in the normalized particle size distribution to determine modality.

    Parameters:
    - diameters: Array of particle diameters.
    - diff_param: Normalized differential distribution.
    - prominence: Required prominence of peaks (passed to scipy.signal.find_peaks).
      Lower it to capture small/low-amplitude peaks (default: 0.001).
    - distance: Required minimal horizontal distance (in samples) between peaks.
      Lower it to resolve peaks that sit close together (default: 20).
    - height: Optional required peak height (default: None = no height filter).
    - width: Optional required peak width in samples (default: None = no width filter).
    - threshold: Optional required vertical distance to neighbouring samples
      (default: None = no threshold filter).

    Returns:
    - num_modes: Number of detected peaks (modes).
    - peak_sizes: Diameters at peak locations.
    """
    peaks, _ = find_peaks(diff_param, prominence=prominence, distance=distance,
                          height=height, width=width, threshold=threshold)
    num_modes = len(peaks)
    peak_sizes = diameters[peaks]
    return num_modes, peak_sizes