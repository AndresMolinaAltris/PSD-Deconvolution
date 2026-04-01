import numpy as np
from scipy.stats import lognorm
from scipy.optimize import curve_fit
from scipy.integrate import cumulative_trapezoid
import pandas as pd


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

# def monomodal_lognorm(x, mu1, sigma1):
#     """
#     Monomodal lognormal distribution (1 mode).
#     Parameters: mu1 (log-mean), sigma1 (shape parameter).
#     Weight is implicitly 1.0.
#     """
#     result = lognorm.pdf(x, s=sigma1, loc=0, scale=np.exp(mu1))
#     return result

# def bimodal_lognorm(x, w1, mu1, sigma1, mu2, sigma2):
#     """
#     Bimodal lognormal distribution (2 modes).
#     Parameters: w1 (weight of first component), mu1, sigma1, mu2, sigma2.
#     Weight w2 is computed as 1 - w1.
#     """
#     weights = [max(0, w1), max(0, 1.0 - w1)]
#     weight_sum = sum(weights)
#     if weight_sum == 0:
#         return np.zeros_like(x, dtype=float)
#     weights = [w / weight_sum for w in weights]
#     result = (weights[0] * lognorm.pdf(x, s=sigma1, loc=0, scale=np.exp(mu1)) +
#               weights[1] * lognorm.pdf(x, s=sigma2, loc=0, scale=np.exp(mu2)))
#     return result

# def trimodal_lognorm(x, w1, w2, mu1, sigma1, mu2, sigma2, mu3, sigma3):
#     """
#     Trimodal lognormal distribution (3 modes).
#     Parameters: w1, w2 (weights of first two components), mu1, sigma1, mu2, sigma2, mu3, sigma3.
#     Weight w3 is computed as 1 - (w1 + w2).
#     """
#     weights = [max(0, w1), max(0, w2), max(0, 1.0 - (w1 + w2))]
#     weight_sum = sum(weights)
#     if weight_sum == 0:
#         return np.zeros_like(x, dtype=float)
#     weights = [w / weight_sum for w in weights]
#     result = (weights[0] * lognorm.pdf(x, s=sigma1, loc=0, scale=np.exp(mu1)) +
#               weights[1] * lognorm.pdf(x, s=sigma2, loc=0, scale=np.exp(mu2)) +
#               weights[2] * lognorm.pdf(x, s=sigma3, loc=0, scale=np.exp(mu3)))
#     return result

# def compute_sigma(diff_param, weighting_scheme=None):
#     """
#     Compute the sigma array for weighted least squares based on the specified scheme.
    
#     Parameters:
#     - diff_param: Normalized differential distribution (y-values).
#     - weighting_scheme: Optional string to specify weighting ('proportional_to_y' for weights ~ y + epsilon, None for unweighted).
    
#     Returns:
#     - sigma: Array of standard deviations for each data point (or None for unweighted).
#     """
#     sigma = None
#     if weighting_scheme == 'proportional_to_y':
#         # Choose epsilon based on data scale
#         max_val = np.max(diff_param[diff_param > 0]) if np.any(diff_param > 0) else 1.0
#         epsilon = 1e-4 * max_val  # Increased epsilon to reduce influence of near-zero values
#         # Threshold for very low values
#         threshold = 1e-3 * max_val  # Higher threshold to identify negligible data points
#         weights = np.where(diff_param > threshold, diff_param + epsilon, 1e-8)
#         print(f'Weights 1 are {weights}')
#         # Ensure weights are positive and clip to avoid numerical issues
#         weights = np.maximum(weights, 1e-8)
#         print(f'Weights 2 are {weights}')
#         sigma = 1 / np.sqrt(weights)
#         sigma = np.clip(sigma, None, 1e6)  # Prevent numerical overflow
#     # Add more schemes here in the future
#     return sigma

def compute_sigma(diff_param, weighting_scheme=None):
    """
    Compute the sigma array for weighted least squares based on the specified scheme.
    
    Parameters:
    - diff_param: Normalized differential distribution (y-values).
    - weighting_scheme: Optional string to specify weighting ('proportional_to_y' for weights ~ y + epsilon, 'sqrt_y' for weights ~ sqrt(y + epsilon), None for unweighted).
    
    Returns:
    - sigma: Array of standard deviations for each data point (or None for unweighted).
    """

    # Convert diff_param to a numpy array with float dtype to avoid dtype=object issues
    diff_param = np.asarray(diff_param, dtype=float)


    sigma = None
    max_val = np.max(diff_param[diff_param > 0]) if np.any(diff_param > 0) else 1.0
    if weighting_scheme == 'proportional_to_y':
        epsilon = 1e-4 * max_val
        threshold = 1e-3 * max_val  # Adjust as needed to downweight more small values
        weights = np.where(diff_param > threshold, diff_param + epsilon, 1e-8)
        weights = np.maximum(weights, 1e-8)  # Simplified safeguard
        sigma = 1 / np.sqrt(weights)
        sigma = np.clip(sigma, None, 1e6)
    elif weighting_scheme == 'sqrt_y':
        epsilon = 1e-4 * max_val
        threshold = 1e-3 * max_val  # Same threshold for consistency
        #print(f'Threshold is is {threshold}')
        #print(np.sqrt(diff_param + epsilon))
        weights = np.where(diff_param > threshold, np.sqrt(diff_param + epsilon), 1e-8)
        #print(f'Weights 1 is {weights}')
        weights = np.maximum(weights, 1e-8)  # Simplified safeguard
        #print(f'weights 2 is {weights}')
        sigma = 1 / np.sqrt(weights)
        sigma = np.clip(sigma, None, 1e6)
    return sigma

# def fit_multimodal_lognorm(diameters, diff_param, num_modes, peak_sizes, fitting_param):
#     """
#     Fit a multimodal lognormal distribution to the data and compute cumulative distribution.
    
#     Parameters:
#     - diameters: Array of particle diameters.
#     - diff_param: Normalized differential distribution.
#     - num_modes: Number of modes (1, 2, or 3).
#     - peak_sizes: Diameters at peak locations.
#     - fitting_param: Column name for the differential distribution (e.g., 'Size distribution Volume weighted [%]').
    
#     Returns:
#     - fit_params: Dictionary of fitted parameters (weights, means, sigmas, R-squared).
#     - fitted_df: DataFrame with 'Particle diameter [µm]', differential, and cumulative columns.
#     - r_squared: R-squared value of the fit.
#     """
#     num_modes = min(num_modes, 3)  # Cap at trimodal
#     modes = peak_sizes[:num_modes]
#     default_sigma = 0.5
#     p0 = []
#     total_height = sum(diff_param[np.searchsorted(diameters, modes)]) if len(modes) > 0 else 1.0
#     initial_weights = []
#     for i, mode in enumerate(modes):
#         if i < num_modes - 1:
#             peak_height = diff_param[np.searchsorted(diameters, mode)] if i < len(modes) else 0
#             weight = peak_height / total_height if total_height > 0 else 1.0 / num_modes
#             weight = max(0, min(weight, 0.99 / max(1, num_modes - 1)))
#             initial_weights.append(weight)
#         mu_init = np.log(mode + 1e-6) + default_sigma**2
#         p0.extend([weight] if i < num_modes - 1 else [])
#         p0.extend([mu_init, default_sigma])

#     if num_modes > 1 and initial_weights:
#         weight_sum = sum(initial_weights)
#         if weight_sum > 0.99:
#             initial_weights = [w * 0.99 / weight_sum for w in initial_weights]
#         p0[:num_modes-1] = initial_weights

#     if num_modes == 1:
#         fit_func = monomodal_lognorm
#         lower_bounds = [-np.inf, 0]
#         upper_bounds = [np.inf, np.inf]
#     elif num_modes == 2:
#         fit_func = bimodal_lognorm
#         lower_bounds = [0, -np.inf, 0, -np.inf, 0]
#         upper_bounds = [1, np.inf, np.inf, np.inf, np.inf]
#     else:
#         fit_func = trimodal_lognorm
#         lower_bounds = [0, 0, -np.inf, 0, -np.inf, 0, -np.inf, 0]
#         upper_bounds = [1, 1, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

#     try:
#         params, _ = curve_fit(fit_func, diameters, diff_param, p0=p0, 
#                               bounds=(lower_bounds, upper_bounds), maxfev=10000)
#     except RuntimeError:
#         print("Curve fitting failed, using initial guesses")
#         params = p0

#     if num_modes == 1:
#         weights = [1.0]
#     else:
#         weights = list(params[:(num_modes-1)]) + [1.0 - sum(params[:(num_modes-1)])]

#     fit_params = {}
#     for i in range(num_modes):
#         if num_modes == 1:
#             fit_params[f'Weight{i+1}'] = weights[i]
#             fit_params[f'Mean{i+1}'] = np.exp(params[i*2])
#             fit_params[f'Sigma{i+1}'] = params[i*2 + 1]
#         else:
#             fit_params[f'Weight{i+1}'] = max(0, weights[i])
#             fit_params[f'Mean{i+1}'] = np.exp(params[(num_modes-1) + i*2])
#             fit_params[f'Sigma{i+1}'] = params[(num_modes-1) + i*2 + 1]

#     y_data = diff_param
#     fitted_psd = fit_func(diameters, *params)
#     ss_res = np.sum((y_data - fitted_psd) ** 2)
#     ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
#     r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
#     fit_params['R_squared'] = r_squared

#     # Compute cumulative fitted distribution
#     cumulative_fitted = cumulative_trapezoid(fitted_psd, diameters, initial=0)
#     cumulative_fitted = 100 * cumulative_fitted / cumulative_fitted[-1] if cumulative_fitted[-1] > 0 else cumulative_fitted

#     # Create DataFrame for fitted data
#     cum_col = fitting_param.replace('Size distribution', 'Undersize')
#     fitted_df = pd.DataFrame({
#         'Particle diameter  [µm]': diameters,
#         fitting_param: fitted_psd,
#         cum_col: cumulative_fitted
#     })

#     # Print weights for each component
#     for i in range(num_modes):
#         print(f"Weight of Component {i+1}: {weights[i]*100:.1f}%")

#     return fit_params, fitted_df, r_squared


def fit_multimodal_lognorm(diameters_x_data, y_psd_data, num_modes, peak_sizes, fitting_param, weighting_scheme=None):
    """
    Fit a multimodal lognormal distribution to the data and compute cumulative distribution.
    
    Parameters:
    - diameters: Array of particle diameters.
    - diff_param: Normalized differential distribution.
    - num_modes: Number of modes (1, 2, or 3).
    - peak_sizes: Diameters at peak locations.
    - fitting_param: Column name for the differential distribution (e.g., 'Size distribution Volume weighted [%]').
    - weighting_scheme: Optional string to specify weighting (passed to compute_sigma).
    
    Returns:
    - fit_params: Dictionary of fitted parameters (weights, means, sigmas, R-squared).
    - fitted_df: DataFrame with 'Particle diameter [µm]', differential, and cumulative columns.
    - r_squared: R-squared value of the fit.
    """
    num_modes = min(num_modes, 3)  # Cap at trimodal
    modes = peak_sizes[:num_modes]
    default_sigma = 0.5
    p0 = []
    total_height = sum(y_psd_data[np.searchsorted(diameters_x_data, modes)]) if len(modes) > 0 else 1.0
    initial_weights = []
    for i, mode in enumerate(modes):
        if i < num_modes - 1:
            peak_height = y_psd_data[np.searchsorted(diameters_x_data, mode)] if i < len(modes) else 0
            weight = peak_height / total_height if total_height > 0 else 1.0 / num_modes
            weight = max(0, min(weight, 0.99 / max(1, num_modes - 1)))
            initial_weights.append(weight)
        mu_init = np.log(mode + 1e-6) + default_sigma**2
        p0.extend([weight] if i < num_modes - 1 else [])
        p0.extend([mu_init, default_sigma])

    if num_modes > 1 and initial_weights:
        weight_sum = sum(initial_weights)
        if weight_sum > 0.99:
            initial_weights = [w * 0.99 / weight_sum for w in initial_weights]
        p0[:num_modes-1] = initial_weights

    if num_modes == 1:
        fit_func = monomodal_lognorm
        lower_bounds = [-np.inf, 0]
        upper_bounds = [np.inf, np.inf]
    elif num_modes == 2:
        fit_func = bimodal_lognorm
        lower_bounds = [0, -np.inf, 0, -np.inf, 0]
        upper_bounds = [1, np.inf, np.inf, np.inf, np.inf]
    else:
        fit_func = trimodal_lognorm
        lower_bounds = [0, 0, -np.inf, 0, -np.inf, 0, -np.inf, 0]
        upper_bounds = [1, 1, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

    # Compute sigma using the separate weighting function
    sigma = compute_sigma(y_psd_data, weighting_scheme)

    try:
        params, _ = curve_fit(fit_func, 
                              diameters_x_data, 
                              y_psd_data, 
                              p0=p0, 
                              bounds=(lower_bounds, upper_bounds), 
                              maxfev=10000, 
                              sigma=sigma)
    except RuntimeError:
        print("Curve fitting failed, using initial guesses")
        params = p0

    if num_modes == 1:
        weights = [1.0]
    else:
        weights = list(params[:(num_modes-1)]) + [1.0 - sum(params[:(num_modes-1)])]

    fit_params = {}
    for i in range(num_modes):
        if num_modes == 1:
            fit_params[f'Weight{i+1}'] = weights[i]
            fit_params[f'Mean{i+1}'] = np.exp(params[i*2])
            fit_params[f'Sigma{i+1}'] = params[i*2 + 1]
        else:
            fit_params[f'Weight{i+1}'] = max(0, weights[i])
            fit_params[f'Mean{i+1}'] = np.exp(params[(num_modes-1) + i*2])
            fit_params[f'Sigma{i+1}'] = params[(num_modes-1) + i*2 + 1]

    y_data = y_psd_data
    fitted_psd = fit_func(diameters_x_data, *params)

    # import matplotlib.pyplot as plt
    # plt.plot(y_data)
    # plt.xscale('log')
    # plt.show




    ss_res = np.sum((y_data - fitted_psd) ** 2)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    fit_params['R_squared'] = r_squared

    # Compute cumulative fitted distribution
    cumulative_fitted = cumulative_trapezoid(fitted_psd, diameters_x_data, initial=0)
    cumulative_fitted = 100 * cumulative_fitted / cumulative_fitted[-1] if cumulative_fitted[-1] > 0 else cumulative_fitted

    # Create DataFrame for fitted data
    cum_col = fitting_param.replace('Size distribution', 'Undersize')
    fitted_df = pd.DataFrame({
        'Particle diameter  [µm]': diameters_x_data,
        fitting_param: fitted_psd,
        cum_col: cumulative_fitted
    })
    
    # Print weights for each component
    for i in range(num_modes):
        print(f"Weight of Component {i+1}: {weights[i]*100:.1f}%")

    return fit_params, fitted_df, r_squared

# Monomodal Cumulative distribution function
def monomodal_lognorm_cdf(x, mu1, sigma1):
    return lognorm.cdf(x, s=sigma1, loc=0, scale=np.exp(mu1))

# Bimodal Cumulative Distribution Function
def bimodal_lognorm_cdf(x, w1, mu1, sigma1, mu2, sigma2):
    weights = [max(0, w1), max(0, 1.0 - w1)]
    weight_sum = sum(weights)
    if weight_sum == 0:
        return np.zeros_like(x, dtype=float)
    weights = [w / weight_sum for w in weights]
    result = (weights[0] * lognorm.cdf(x, s=sigma1, loc=0, scale=np.exp(mu1)) +
              weights[1] * lognorm.cdf(x, s=sigma2, loc=0, scale=np.exp(mu2)))
    return result * 100  # Scale to [%] to match data

# Trimodal Cumulative Distribution Function
def trimodal_lognorm_cdf(x, w1, w2, mu1, sigma1, mu2, sigma2, mu3, sigma3):
    weights = [max(0, w1), max(0, w2), max(0, 1.0 - (w1 + w2))]
    weight_sum = sum(weights)
    if weight_sum == 0:
        return np.zeros_like(x, dtype=float)
    weights = [w / weight_sum for w in weights]
    result = (weights[0] * lognorm.cdf(x, s=sigma1, loc=0, scale=np.exp(mu1)) +
              weights[1] * lognorm.cdf(x, s=sigma2, loc=0, scale=np.exp(mu2)) +
              weights[2] * lognorm.cdf(x, s=sigma3, loc=0, scale=np.exp(mu3)))
    return result * 100  # Scale to [%]


# Cumulative Distribution Function (CDF) fitting function 
def fit_multimodal_lognorm_cdf(diameters, data_df, num_modes, peak_sizes, fitting_param, weighting_scheme=None):
    """
    Fit a multimodal lognormal distribution to the cumulative data and derive differential.
    
    Parameters:
    - diameters: Array of particle diameters.
    - data_df: Preprocessed DataFrame containing both differential and cumulative columns.
    - num_modes: Number of modes (1, 2, or 3).
    - peak_sizes: Diameters at peak locations.
    - fitting_param: Column name for the differential distribution (e.g., 'Size distribution Volume weighted [%]').
    - weighting_scheme: Optional string to specify weighting (passed to compute_sigma).
    
    Returns:
    - fit_params: Dictionary of fitted parameters (weights, means, sigmas, R-squared on cumulative).
    - fitted_df: DataFrame with 'Particle diameter [µm]', differential (derived), and cumulative columns.
    - r_squared: R-squared value of the fit (on cumulative).
    """
    num_modes = min(num_modes, 3)
    cum_col = fitting_param.replace('Size distribution', 'Undersize')
    cum_param = data_df[cum_col].values  # Target: cumulative in [%]

    # Use differential for initial guesses (peaks still useful)
    diff_param = data_df[fitting_param].values
    modes = peak_sizes[:num_modes]
    default_sigma = 0.5
    p0 = []
    total_height = sum(diff_param[np.searchsorted(diameters, modes)]) if len(modes) > 0 else 1.0
    initial_weights = []
    for i, mode in enumerate(modes):
        if i < num_modes - 1:
            peak_height = diff_param[np.searchsorted(diameters, mode)] if i < len(modes) else 0
            weight = peak_height / total_height if total_height > 0 else 1.0 / num_modes
            weight = max(0, min(weight, 0.99 / max(1, num_modes - 1)))
            initial_weights.append(weight)
        mu_init = np.log(mode + 1e-6) + default_sigma**2
        p0.extend([weight] if i < num_modes - 1 else [])
        p0.extend([mu_init, default_sigma])

    if num_modes > 1 and initial_weights:
        weight_sum = sum(initial_weights)
        if weight_sum > 0.99:
            initial_weights = [w * 0.99 / weight_sum for w in initial_weights]
        p0[:num_modes-1] = initial_weights

    # Select CDF function
    if num_modes == 1:
        fit_func = monomodal_lognorm_cdf
        lower_bounds = [-np.inf, 0]
        upper_bounds = [np.inf, np.inf]
    elif num_modes == 2:
        fit_func = bimodal_lognorm_cdf
        lower_bounds = [0, -np.inf, 0, -np.inf, 0]
        upper_bounds = [1, np.inf, np.inf, np.inf, np.inf]
    else:
        fit_func = trimodal_lognorm_cdf
        lower_bounds = [0, 0, -np.inf, 0, -np.inf, 0, -np.inf, 0]
        upper_bounds = [1, 1, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

    # Sigma for weighted fit (normalize cum_param to [0,1] for weighting)
    sigma = compute_sigma(cum_param / 100, weighting_scheme)

    try:
        params, _ = curve_fit(fit_func, diameters, cum_param, p0=p0, 
                              bounds=(lower_bounds, upper_bounds), maxfev=10000, sigma=sigma)
    except RuntimeError:
        print("Curve fitting failed, using initial guesses")
        params = p0

    # Extract params
    if num_modes == 1:
        weights = [1.0]
    else:
        weights = list(params[:(num_modes-1)]) + [1.0 - sum(params[:(num_modes-1)])]

    fit_params = {}
    for i in range(num_modes):
        if num_modes == 1:
            fit_params[f'Weight{i+1}'] = weights[i]
            fit_params[f'Mean{i+1}'] = np.exp(params[i*2])
            fit_params[f'Sigma{i+1}'] = params[i*2 + 1]
        else:
            fit_params[f'Weight{i+1}'] = max(0, weights[i])
            fit_params[f'Mean{i+1}'] = np.exp(params[(num_modes-1) + i*2])
            fit_params[f'Sigma{i+1}'] = params[(num_modes-1) + i*2 + 1]

    # Compute fitted cumulative and derive differential
    cumulative_fitted = fit_func(diameters, *params)
    fitted_psd = np.gradient(cumulative_fitted, diameters) / 100  # Derive PDF (scale back from %)

    # R-squared on cumulative (target)
    ss_res = np.sum((cum_param - cumulative_fitted) ** 2)
    ss_tot = np.sum((cum_param - np.mean(cum_param)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    fit_params['R_squared'] = r_squared

    # Fitted DF
    fitted_df = pd.DataFrame({
        'Particle diameter  [µm]': diameters,
        fitting_param: fitted_psd,  # Derived PSD (not yet normalized)
        cum_col: cumulative_fitted
    })

    # Normalize derived PSD
    integral = np.trapezoid(fitted_df[fitting_param], diameters)
    if integral > 0:
        fitted_df[fitting_param] /= integral

    # Print weights for each component
    for i in range(num_modes):
        print(f"Weight of Component {i+1}: {weights[i]*100:.1f}%")

    return fit_params, fitted_df, r_squared