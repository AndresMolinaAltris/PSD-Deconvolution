import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import lognorm

def plot_distribution(diameters, diff_param, fitted_psd, parameter_fitted):
    """
    Plot the original and fitted particle size distributions.
    
    Parameters:
    - diameters: Array of particle diameters.
    - diff_param: Original differential distribution.
    - fitted_psd: Fitted distribution values.
    - num_modes: Number of modes in the fit.
    """
    plt.figure(figsize=(8, 5))
    plt.plot(diameters, diff_param, 'b-', label='Original Data')
    plt.plot(diameters, fitted_psd, 'r--', label=f'Fitted')
    plt.xscale('log')
    plt.xlim(0, 100)
    plt.xlabel('Particle Diameter [µm]')
    plt.ylabel('Density')
    plt.title(f'Fit vs. Original Data ({parameter_fitted})')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()

def plot_cumulative(data_df, diameters, diff_param, fitting_param, num_modes):
    """
    Plot the original and fitted cumulative distributions.
    
    Parameters:
    - data_df: Preprocessed DataFrame.
    - diameters: Array of particle diameters.
    - diff_param: Original differential distribution.
    - fitting_param: Column name for the fitting parameter.
    - num_modes: Number of modes in the fit.
    """
    cum_col = fitting_param.replace('Size distribution', 'Undersize')
    fitted_psd = diff_param
    cumulative_fitted = cumulative_trapezoid(fitted_psd, diameters, initial=0)
    cumulative_fitted = 100 * cumulative_fitted / cumulative_fitted[-1] if cumulative_fitted[-1] > 0 else cumulative_fitted
    
    plt.figure(figsize=(8, 5))
    plt.plot(data_df['Particle diameter  [µm]'], data_df[cum_col], 'b-', label='Original Cumulative')
    plt.plot(diameters, cumulative_fitted, 'r--', label='Fitted Cumulative')
    #plt.xscale('log')
    plt.xlim(0.1,200)
    plt.xlabel('Particle Diameter [µm]')
    plt.ylabel('Cumulative Distribution [%]')
    plt.title(f'Cumulative vs. Fitted Distribution ({fitting_param})')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()


def plot_individual_modes(diameters, original_psd, fitted_psd, fit_params, num_modes, selected_weighted_data):
    """
    Plot the original PSD, fitted total PSD, and individual mode contributions.
    
    Parameters:
    - diameters: Array of particle diameters (x-values).
    - original_psd: Original normalized PSD (frequency) data.
    - fitted_psd: Fitted total normalized PSD (frequency) data.
    - fit_params: Dictionary of fitted parameters (weights, means, sigmas).
    - num_modes: Number of modes (1, 2, or 3).
    - selected_weighted_data: String describing the data type (e.g., 'Size distribution Volume weighted [%]') for the title.
    """
    plt.figure(figsize=(10, 6))
    
    # Plot original data
    plt.plot(diameters, original_psd, 'k-', label='Original Data')
    
    # Plot fitted total
    #plt.plot(diameters, fitted_psd, 'r--', label='Fitted Total')
    
    # Colors for modes (up to 3)
    colors = ['b', 'g', 'm']
    
    # Plot each individual mode
    for i in range(num_modes):
        weight = fit_params[f'Weight{i+1}']
        mu = np.log(fit_params[f'Mean{i+1}'])  # Convert median back to log-mean if needed, but since Mean is exp(mu), log(Mean) = mu
        sigma = fit_params[f'Sigma{i+1}']
        
        # Compute individual mode: weight * lognorm.pdf
        individual_mode = weight * lognorm.pdf(diameters, s=sigma, loc=0, scale=fit_params[f'Mean{i+1}'])
        
        plt.plot(diameters, individual_mode, color=colors[i], linestyle='--', label=f'Fitted Mode {i+1}')
    
    # Formatting
    plt.xscale('log')
    plt.xlim(0.5,100)
    plt.xlabel('Particle Diameter [µm]')
    plt.ylabel('Normalized PSD')
    plt.title(f'Original vs Fitted PSD {selected_weighted_data}')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()