import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from data_processing import load_data, compute_basic_statistics, find_xlsx_files
from peak_detection import detect_peaks, prepare_distribution_data, find_distribution_peaks, prepare_distribution_data_interpolation
from fitting import fit_multimodal_lognorm, fit_multimodal_lognorm_cdf, monomodal_lognorm
from visualization import plot_distribution, plot_cumulative, plot_individual_modes
from scipy.stats import lognorm
from scipy.integrate import cumulative_trapezoid, trapezoid
import re
import config

def main():
    parser = argparse.ArgumentParser(description='PSD Analysis')
    parser.add_argument('--param_fit', type=int, default=config.PARAM_FIT,
                        help='Weighting index: 0=Surface, 1=Volume, 2=Number (default: config.PARAM_FIT)')
    parser.add_argument('--files', nargs='+', default=None,
                        help='One or more filenames to analyse (e.g. sample1.xlsx sample2.xlsx). '
                             'If omitted, all .xlsx files in data_dir are processed.')
    args = parser.parse_args()

    parent_dir = config.DATA_DIR
    param_fit = args.param_fit
    fitting_param_list = [
        'Size distribution Surface weighted [%]',
        'Size distribution Volume weighted [%]',
        'Size distribution Number weighted [%]'
    ]
    
    # Dynamically extract the weighting based on param_fit (e.g., 'Number weighted' for index 2)
    selected_weighted_data = fitting_param_list[param_fit] 
    selected_parameter = re.search(r'Size distribution (.*?) \[%\]', selected_weighted_data).group(1)
     
    

    # Step 1: Import Data
    if args.files:
        xlsx_files = [str(Path(parent_dir) / f) for f in args.files]
    else:
        xlsx_files = find_xlsx_files(parent_dir)
    #xlsx_files = find_xlsx_files(parent_dir, 'PSD', None, None, 'QC')


    all_results = []
    # Main analysis loop
    for filename in xlsx_files:
        #try:
        path = Path(filename)
        file_ID = path.name
        current_results = {}
        current_results['Filename'] = file_ID
        current_results['Parameter'] = selected_parameter
        print(f'File analysed: {file_ID}')

        # Step 2: Load and preprocess data        
        data_df = load_data(filename, preprocess=True, param_fit=selected_weighted_data, min_diameter=0.01)

        # Step 3: Basic Statistics on the original weighted data
        results = compute_basic_statistics(data_df, selected_parameter)        
        
        print(f'Results before deconvolution ({selected_parameter}): '
            f'D10: {results[selected_parameter]["D10"]:.2f}, '
            f'D50: {results[selected_parameter]["D50"]:.2f}, '
            f'D90: {results[selected_parameter]["D90"]:.2f}, '
            f'Mean: {results[selected_parameter]["Mean"]:.2f}, '
            f'Std: {results[selected_parameter]["Std"]:.2f}')
        
        current_results['D10_raw'] = round(results[selected_parameter]["D10"], 3)
        current_results['D50_raw'] = round(results[selected_parameter]["D50"],3)
        current_results['D90_raw'] = round(results[selected_parameter]["D90"],3)

        # Step 4: Peak Detection. This works on the frequency and not the cumulative data
        # Here I prepare the data by using interpolation and normalizing it to 1
        diameters_x_data, y_psd_data_normalized = prepare_distribution_data_interpolation(data_df, selected_weighted_data)

        # Here I use the find peaks algorithm. This could be changed. CUrrently the finding peaks only works for frequency data.
        # Test with second differential
        num_modes, peak_sizes = find_distribution_peaks (diameters_x_data, y_psd_data_normalized)
        print(f'number of peaks: {num_modes}')
        print(f'peak_sizes are {peak_sizes}')
        current_results['num_modes'] = num_modes
        current_results['peak_sizes'] = peak_sizes

        # Step 5: Multimodal Fitting
        fit_params, fitted_df, r_squared = fit_multimodal_lognorm(diameters_x_data, 
                                                                y_psd_data_normalized, 
                                                                num_modes, 
                                                                peak_sizes, 
                                                                selected_weighted_data,
                                                                weighting_scheme='proportional_to_y')
        print(f'R-squared is {r_squared}')
        current_results['R_squared'] = round(r_squared, 3)
        

        # Compute the statistics of each mode
        diameters = fitted_df['Particle diameter  [µm]'].values

        for i in range(1, num_modes + 1):
            # Extract params for this mode (note: Mean is the scale, mu = log(Mean))
            weight = fit_params[f'Weight{i}']
            mu = np.log(fit_params[f'Mean{i}'])
            sigma = fit_params[f'Sigma{i}']
            
            # Compute the individual mode's contribution to PSD (weighted)
            #mode_psd = weight * lognorm.pdf(diameters, s=sigma, loc=0, scale=np.exp(mu))
            calculated_psd = weight * monomodal_lognorm(diameters, mu, sigma)
            
            # Normalize to treat as standalone distribution (integrates to 1)
            integral = trapezoid(calculated_psd, diameters)
            if integral > 0:
                normalized_mode_psd = calculated_psd / integral
            else:
                normalized_mode_psd = np.zeros_like(diameters)
                print(f"Warning: Mode {i} has zero integral; skipping statistics.")
                continue
            
            # Compute cumulative for the mode (in %)
            mode_cumulative = cumulative_trapezoid(normalized_mode_psd, diameters, initial=0) * 100
            
            # Create a temp DataFrame matching the structure expected by compute_basic_statistics
            diff_col = f'Size distribution {selected_parameter} [%]'
            cum_col = f'Undersize {selected_parameter} [%]'
            temp_df = pd.DataFrame({
                'Particle diameter  [µm]': diameters,
                diff_col: normalized_mode_psd * 100,  # Scale to % for consistency
                cum_col: mode_cumulative
            })
            
            # Compute statistics on this mode
            mode_results = compute_basic_statistics(temp_df, selected_parameter)
            
            # Print results
            print(f'Results for Mode {i} ({selected_parameter}): '
                f'D10: {mode_results[selected_parameter]["D10"]:.2f}, '
                f'D50: {mode_results[selected_parameter]["D50"]:.2f}, '
                f'D90: {mode_results[selected_parameter]["D90"]:.2f}, '
                f'Mean: {mode_results[selected_parameter]["Mean"]:.2f}, '
                f'Std: {mode_results[selected_parameter]["Std"]:.2f}')
                       
            # Store mode-specific results
            current_results[f'Mode{i}_Weight'] = round(weight, 3)
            current_results[f'Mode{i}_D10'] = round(mode_results[selected_parameter]["D10"], 3)
            current_results[f'Mode{i}_D50'] = round(mode_results[selected_parameter]["D50"], 3)
            current_results[f'Mode{i}_D90'] = round(mode_results[selected_parameter]["D90"], 3)
            current_results[f'Mode{i}_Mean'] = round(mode_results[selected_parameter]["Mean"], 3)
            current_results[f'Mode{i}_Std'] = round(mode_results[selected_parameter]["Std"], 3)
        
        # Step 6: Visualization
        # Call the new plot_individual_modes function
        plot_individual_modes(
            diameters_x_data, 
            y_psd_data_normalized, 
            fitted_df[selected_weighted_data], 
            fit_params, 
            num_modes, 
            selected_weighted_data
        )

        # Step 6: Visualization
        #plot_distribution(diameters, y_psd_data_normalized, fitted_df[fitting_param_list[param_fit]], selected_parameter)
        #plot_cumulative(data_df, diameters, fitted_df[fitting_param_list[param_fit]], fitting_param_list[param_fit], num_modes)
        all_results.append(current_results)
        #except: 
            #print(f'Error with file {file_ID}')
        
    all_results_df = pd.DataFrame(all_results)
    output_path = Path(config.OUTPUT_DIR) / config.OUTPUT_FILENAME
    all_results_df.to_excel(output_path, index=False)
    print(f'\nResults saved to: {output_path}')


if __name__ == "__main__":
    main()
