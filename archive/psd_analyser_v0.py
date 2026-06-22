import pandas as pd
import os
from scipy.stats import lognorm
from scipy.signal import find_peaks, find_peaks_cwt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import re
from pathlib import Path
from scipy.integrate import cumulative_trapezoid
from functions import compute_percentiles, monomodal_lognorm, bimodal_lognorm, trimodal_lognorm, find_xlsx_files



######################################
### IMPORT DATA ######################
######################################
# Parent directory
parent_dir = 'C:\\Users\\AndresMolina\\Documents\\03_Projects\\15_PSD_Analysis\\Data'
#parent_dir = 'C:\\Users\\AndresMolina\\OneDrive - Altris AB\MD - Kungsgatan - 600L 2025'

xlsx_files = find_xlsx_files(parent_dir)


# Here we select which mode we want to fit
fitting_param_list = ['Size distribution Surface weighted [%]', 
                      'Size distribution Volume weighted [%]',
                      'Size distribution Number weighted [%]']
param_fit = 1

# Main analysis loop
#for filename in psd_qc_files:
for filename in xlsx_files: 
    path = Path(filename)
    file_ID = path.name
    print(f'File analysed:{file_ID}')

    headers = pd.read_excel(filename, usecols="F:V", nrows=3, skiprows=4)
    new_cols = headers.astype(str).apply(lambda x: ' '.join(x), axis=0)
    data_df = pd.read_excel(filename, usecols="F:V", skiprows=5, decimal=".")
    data_df.columns = new_cols
    data_df.columns = data_df.columns.str.replace("nan", "", regex=False).str.strip()
    data_df = data_df.drop(index=[0,1]).reset_index(drop=True)

    data_df = data_df.sort_values('Particle diameter  [µm]')
    data_df = data_df.dropna(subset=['Particle diameter  [µm]', 'Undersize Volume weighted [%]', 'Size distribution Volume weighted [%]'])

    ################################################
    ### Step 1: Basic Statistics and Percentiles ###
    ################################################
    weightings = ['Volume weighted'] # Volume weighted is the current value being analysed in our team
    results = {}
    for weight in weightings:
        cum_col = f'Undersize {weight} [%]'
        diff_col = f'Size distribution {weight} [%]'
        d10, d50, d90 = compute_percentiles(data_df['Particle diameter  [µm]'], data_df[cum_col])
        mean = np.average(data_df['Particle diameter  [µm]'], weights=data_df[diff_col])
        median = d50
        std = np.sqrt(np.average((data_df['Particle diameter  [µm]'] - mean)**2, weights=data_df[diff_col]))
        results[weight] = {'D10': d10, 'D50': d50, 'D90': d90, 'Mean': mean, 'Std': std}
    print(f'Results before deconvolution (volume weighted): D10: {d10}, D50: {d50}, D90: {d90}')

    ############################################################################
    ### Step 2: Multimodality Detection (using Volume weighted as example) #####
    ############################################################################
    data_df['Particle diameter  [µm]'] = pd.to_numeric(data_df['Particle diameter  [µm]'], errors='coerce')
    diameters = data_df['Particle diameter  [µm]'].values
    diameters = np.array(diameters)
    log_diameters = np.log(diameters + 1e-6)

    diff_param = data_df[fitting_param_list[param_fit]].values
    diff_param = diff_param / np.trapezoid(diff_param, diameters)
    peaks, properties = find_peaks(diff_param, prominence=0.001, distance=5)
    num_modes = len(peaks)
    print(f'Number of peaks {num_modes}')
    peak_sizes = diameters[peaks]

    ############################################################################
    #### Step 3: Fit Multimodal Lognormal ######################################
    ############################################################################
    if num_modes >= 1:  # Fit for 1, 2, or 3 modes
        peaks = sorted(peaks)
        modes = [diameters[p] for p in peaks]
        num_modes = min(num_modes, 3)  # Cap at trimodal

        # Initialize parameters based on num_modes
        default_sigma = 0.5
        p0 = []
        total_height = sum(diff_param[peaks]) if len(peaks) > 0 else 1.0
        initial_weights = []
        for i, mode in enumerate(modes[:num_modes]):
            if i < num_modes - 1:  # Only add weights for n-1 modes
                peak_height = diff_param[peaks[i]] if i < len(peaks) else 0
                weight = peak_height / total_height if total_height > 0 else 1.0 / num_modes
                weight = max(0, min(weight, 0.99 / max(1, num_modes - 1)))
                initial_weights.append(weight)
            mu_init = np.log(mode + 1e-6) + default_sigma**2
            p0.extend([weight] if i < num_modes - 1 else [])
            p0.extend([mu_init, default_sigma])

        # Normalize initial weights to sum to < 1
        if num_modes > 1 and initial_weights:
            weight_sum = sum(initial_weights)
            if weight_sum > 0.99:
                initial_weights = [w * 0.99 / weight_sum for w in initial_weights]
            p0[:num_modes-1] = initial_weights

        # Define bounds based on num_modes and define fit_func
        if num_modes == 1:
            fit_func = monomodal_lognorm
            lower_bounds = [-np.inf, 0]  # mu1, sigma1
            upper_bounds = [np.inf, np.inf]
        elif num_modes == 2:
            fit_func = bimodal_lognorm
            lower_bounds = [0, -np.inf, 0, -np.inf, 0]  # w1, mu1, sigma1, mu2, sigma2
            upper_bounds = [1, np.inf, np.inf, np.inf, np.inf]
        else:  # num_modes == 3
            fit_func = trimodal_lognorm
            lower_bounds = [0, 0, -np.inf, 0, -np.inf, 0, -np.inf, 0]  # w1, w2, mu1, sigma1, mu2, sigma2, mu3, sigma3
            upper_bounds = [1, 1, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

        # Perform fit
        try:
            params, _ = curve_fit(fit_func, diameters, diff_param, p0=p0, 
                                bounds=(lower_bounds, upper_bounds), maxfev=10000)
        except RuntimeError:
            print("Curve fitting failed, using initial guesses")
            params = p0

        # Reconstruct weights for output (sum to 1)
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

        # Calculate R-squared
        y_data = diff_param
        y_pred = fit_func(diameters, *params)
        ss_res = np.sum((y_data - y_pred) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        fit_params['R_squared'] = r_squared
        print(f'R-squared is {r_squared}')

        # Compute statistics for each component
        for i in range(num_modes):
            if num_modes == 1:
                dist = lognorm(s=params[i*2 + 1], loc=0, scale=np.exp(params[i*2]))
            else:
                dist = lognorm(s=params[(num_modes-1) + i*2 + 1], 
                            loc=0, scale=np.exp(params[(num_modes-1) + i*2]))
            d10 = dist.ppf(0.1)
            d50 = dist.ppf(0.5)
            d90 = dist.ppf(0.9)
            print(f"Component {i+1}: D10={d10:.2f} µm, D50={d50:.2f} µm, D90={d90:.2f} µm")
            print(f"Weight of Component {i+1}: {weights[i]*100:.1f}%")

        # Plotting
        plt.figure(figsize=(8, 5))
        plt.plot(diameters, diff_param, 'b-', label='Original Data')
        plt.plot(diameters, fit_func(diameters, *params), 'r--', label=f'Fitted {num_modes}-Modal')
        plt.xscale('log')
        plt.xlim(0,100)
        plt.xlabel('Particle Diameter [µm]')
        plt.ylabel('Density')
        plt.title(f'Log-Normal Fit vs. Original Data (Log Scale, {num_modes} modes)')
        plt.legend()
        plt.grid(True, which="both", ls="--")
        plt.show()

        #plt.xscale('log')
        #plt.yscale('log')
    
        # Plot individual components
        plt.figure(figsize=(8, 5))
        #plt.plot(diameters, diff_param, 'b-', label='Original Data')
        #plt.plot(diameters, fit_func(diameters, *params), 'r--', label=f'Fitted {num_modes}-Modal')

        # Compute total fitted PSD
        fitted_psd = fit_func(diameters, *params)
        # Integrate fitted PSD to get cumulative distribution
        cumulative_fitted = cumulative_trapezoid(fitted_psd, diameters, initial=0)
        # Normalize to match percentage scale (0-100)
        cumulative_fitted = 100 * cumulative_fitted / cumulative_fitted[-1] if cumulative_fitted[-1] > 0 else cumulative_fitted
        
        
        # This column is originally in the data and represent the integrated version of the Weighted values
        # Dynamic column selection
        cum_col = fitting_param_list[param_fit].replace('Size distribution', 'Undersize')  
        
        # Plot original cumulative data and fitted cumulative using selected weighting
        # Integral plot for total PSD curve
        plt.figure(figsize=(8, 5))
        plt.plot(data_df['Particle diameter  [µm]'], data_df[cum_col], 'b-', label='Original Cumulative')
        plt.plot(diameters, cumulative_fitted, 'r--', label='Fitted Cumulative')
        plt.xscale('log')
        plt.xlabel('Particle Diameter [µm]')
        plt.ylabel('Cumulative Distribution [%]')
        plt.title(f'Cumulative Distribution vs. Fitted Integral (Log Scale, {num_modes} modes, {fitting_param_list[param_fit]})')
        plt.legend()
        plt.grid(True, which="both", ls="--")
        plt.show()



        