from scipy.stats import lognorm
from scipy.interpolate import interp1d
import numpy as np
import os

def compute_percentiles(diameters, cum_weights, percentiles=[10, 50, 90]):
    # Linear interpolation on cumulative distribution
    interp_func = interp1d(cum_weights, diameters, bounds_error=False, fill_value=(diameters.min(), diameters.max()))
    return [float(interp_func(p)) for p in percentiles]

def monomodal_lognorm(x, mu1, sigma1):
    """
    Monomodal lognormal distribution (1 mode).
    Parameters: mu1 (log-mean), sigma1 (shape parameter).
    Weight is implicitly 1.0.
    """
    result = lognorm.pdf(x, s=sigma1, loc=0, scale=np.exp(mu1))
    return result

def bimodal_lognorm(x, w1, mu1, sigma1, mu2, sigma2):
    """
    Bimodal lognormal distribution (2 modes).
    Parameters: w1 (weight of first component), mu1, sigma1, mu2, sigma2.
    Weight w2 is computed as 1 - w1.
    """
    weights = [max(0, w1), max(0, 1.0 - w1)]
    weight_sum = sum(weights)
    if weight_sum == 0:
        return np.zeros_like(x, dtype=float)
    weights = [w / weight_sum for w in weights]
    result = (weights[0] * lognorm.pdf(x, s=sigma1, loc=0, scale=np.exp(mu1)) +
              weights[1] * lognorm.pdf(x, s=sigma2, loc=0, scale=np.exp(mu2)))
    return result

def trimodal_lognorm(x, w1, w2, mu1, sigma1, mu2, sigma2, mu3, sigma3):
    """
    Trimodal lognormal distribution (3 modes).
    Parameters: w1, w2 (weights of first two components), mu1, sigma1, mu2, sigma2, mu3, sigma3.
    Weight w3 is computed as 1 - (w1 + w2).
    """
    weights = [max(0, w1), max(0, w2), max(0, 1.0 - (w1 + w2))]
    weight_sum = sum(weights)
    if weight_sum == 0:
        return np.zeros_like(x, dtype=float)
    weights = [w / weight_sum for w in weights]
    result = (weights[0] * lognorm.pdf(x, s=sigma1, loc=0, scale=np.exp(mu1)) +
              weights[1] * lognorm.pdf(x, s=sigma2, loc=0, scale=np.exp(mu2)) +
              weights[2] * lognorm.pdf(x, s=sigma3, loc=0, scale=np.exp(mu3)))
    return result


def find_xlsx_files(parent_dir: str, 
                    subfolder_keyword: str = None, 
                    sub_subfolder_keyword: str = None, 
                    exact_subfolder: str = None, 
                    filename_keyword: str = None):
    """
    Finds .xlsx files in the given directory, with optional filters for subfolders and filename patterns.
    
    Parameters:
    - parent_dir: The root directory to search in.
    - subfolder_keyword: If provided, only search in subfolders where this keyword (case-insensitive) appears in the path components.
    - sub_subfolder_keyword: If provided, only search in subfolders where this additional keyword (case-insensitive) appears in the path components.
    - exact_subfolder: If provided, only search in subfolders with this exact name (case-insensitive).
    - filename_keyword: If provided, filter files where the second part after '-' in the filename (after stripping spaces) equals this keyword (case-insensitive).
    
    Returns:
    - A list of full paths to matching .xlsx files.
    """
    
    xlsx_files = []
    for root, dirs, files in os.walk(parent_dir):
        path_parts = root.lower().split(os.sep)
        
        if subfolder_keyword and subfolder_keyword.lower() not in path_parts:
            continue
        
        if sub_subfolder_keyword and sub_subfolder_keyword.lower() not in path_parts:
            continue
        
        if exact_subfolder and os.path.basename(root).lower() != exact_subfolder.lower():
            continue
        
        for file in files:
            if file.lower().endswith('.xlsx'):
                xlsx_files.append(os.path.join(root, file))
    
    if filename_keyword:
        filtered_files = []
        for file_path in xlsx_files:
            file_name = os.path.basename(file_path)
            file_name_no_ext = os.path.splitext(file_name)[0]
            parts = file_name_no_ext.split('-')
            if len(parts) >= 2:
                parts = [part.replace(' ', '') for part in parts]
                if parts[1].lower() == filename_keyword.lower():
                    filtered_files.append(file_path)
        return filtered_files
    else:
        return xlsx_files