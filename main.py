import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from data_processing import load_data, compute_basic_statistics, find_xlsx_files
from peak_detection import find_distribution_peaks, prepare_distribution_data_interpolation
from fitting import fit_multimodal_lognorm, monomodal_lognorm
from visualization import plot_individual_modes
from scipy.integrate import cumulative_trapezoid, trapezoid
import re
import config

def main():
    parser = argparse.ArgumentParser(
        description='PSD Analysis. Optional fitting arguments default to the values in '
                    'config.py; omit them to reproduce the standard fit.')
    parser.add_argument('--param_fit', type=int, default=config.PARAM_FIT,
                        help='Weighting index: 0=Surface, 1=Volume, 2=Number '
                             f'(default: {config.PARAM_FIT})')
    parser.add_argument('--files', nargs='+', default=None,
                        help='One or more filenames to analyse (e.g. sample1.xlsx sample2.xlsx). '
                             'If omitted, all .xlsx files in data_dir are processed.')
    parser.add_argument('--min_diameter', type=float, default=config.MIN_DIAMETER,
                        help='Lower particle-diameter cutoff in µm; rows below this are '
                             f'dropped before fitting (default: {config.MIN_DIAMETER})')
    parser.add_argument('--prominence', type=float, default=config.PEAK_REL_PROMINENCE,
                        help='Peak prominence for mode detection, as a fraction of the '
                             'distribution peak height (relative, adapts per file). Lower '
                             f'it to catch smaller modes (default: {config.PEAK_REL_PROMINENCE})')
    parser.add_argument('--distance', type=int, default=config.PEAK_DISTANCE,
                        help='Minimum separation (in samples) between detected peaks; lower '
                             f'it to resolve close peaks (default: {config.PEAK_DISTANCE})')
    parser.add_argument('--weighting', choices=['proportional_to_y', 'sqrt_y', 'none'],
                        default=config.WEIGHTING_SCHEME,
                        help='Least-squares weighting scheme for the fit. "none" (unweighted) '
                             'best resolves small low-amplitude modes; "proportional_to_y" '
                             f'follows the dominant peak (default: {config.WEIGHTING_SCHEME})')
    parser.add_argument('--num_modes', type=int, default=None,
                        help='Force the number of modes to fit (1-3), overriding automatic '
                             'peak detection. If omitted, the mode count is detected.')
    args = parser.parse_args()

    # 'none' on the command line maps to Python None (unweighted fit).
    weighting_scheme = None if args.weighting == 'none' else args.weighting

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
        xlsx_files = []
        for f in args.files:
            path = Path(parent_dir) / f
            if path.exists():
                xlsx_files.append(str(path))
            else:
                print(f'Warning: file not found, skipping: {path}')
    else:
        xlsx_files = find_xlsx_files(parent_dir)

    if not xlsx_files:
        print(f'No .xlsx files to analyse in {parent_dir}. Nothing to do.')
        return

    all_results = []
    # Main analysis loop. Each file is isolated so one failure does not abort the batch.
    for filename in xlsx_files:
        path = Path(filename)
        file_ID = path.name
        try:
            current_results = {}
            current_results['Filename'] = file_ID
            current_results['Parameter'] = selected_parameter
            print(f'File analysed: {file_ID}')

            # Step 2: Load and preprocess data
            data_df = load_data(filename, preprocess=True, param_fit=selected_weighted_data,
                                min_diameter=args.min_diameter)

            # Step 3: Basic Statistics on the original weighted data
            results = compute_basic_statistics(data_df, selected_parameter)

            print(f'Results before deconvolution ({selected_parameter}): '
                f'D10: {results[selected_parameter]["D10"]:.2f}, '
                f'D50: {results[selected_parameter]["D50"]:.2f}, '
                f'D90: {results[selected_parameter]["D90"]:.2f}, '
                f'Mean: {results[selected_parameter]["Mean"]:.2f}, '
                f'Std: {results[selected_parameter]["Std"]:.2f}')

            current_results['D10_raw'] = round(results[selected_parameter]["D10"], 3)
            current_results['D50_raw'] = round(results[selected_parameter]["D50"], 3)
            current_results['D90_raw'] = round(results[selected_parameter]["D90"], 3)

            # Step 4: Peak Detection. Works on the frequency data, not the cumulative.
            # Prepare the data via interpolation, normalized so it integrates to 1.
            diameters_x_data, y_psd_data_normalized = prepare_distribution_data_interpolation(
                data_df, selected_weighted_data, interp_points=config.INTERP_POINTS)

            # Detect the modes (peaks) of the distribution. --prominence is relative
            # (a fraction of the peak height), so convert it to the absolute value
            # scipy.signal.find_peaks expects for this particular distribution.
            abs_prominence = args.prominence * float(np.max(y_psd_data_normalized))
            num_modes, peak_sizes = find_distribution_peaks(
                diameters_x_data, y_psd_data_normalized,
                prominence=abs_prominence, distance=args.distance)
            print(f'number of peaks: {num_modes}')
            print(f'peak_sizes are {peak_sizes}')

            # Optional override of the detected mode count (clamped to 1-3).
            if args.num_modes is not None:
                forced = max(1, min(args.num_modes, 3))
                if forced != num_modes:
                    print(f'Overriding detected modes ({num_modes}) with --num_modes {forced}')
                # Seed any extra modes in the small-particle region below the smallest
                # detected peak. In this data the missing mode is almost always the
                # small (~2-5 µm) population, so biasing the seeds low resolves it
                # far more reliably than seeding across the full diameter range.
                if len(peak_sizes) < forced:
                    extra_needed = forced - len(peak_sizes)
                    upper = float(np.min(peak_sizes)) if len(peak_sizes) else float(diameters_x_data[-1])
                    seeds = np.logspace(np.log10(diameters_x_data[0]),
                                        np.log10(upper),
                                        extra_needed + 2)[1:-1]
                    peak_sizes = np.sort(np.concatenate(
                        [np.asarray(peak_sizes, dtype=float), seeds]))
                peak_sizes = peak_sizes[:forced]
                num_modes = forced
            current_results['num_modes'] = num_modes
            current_results['peak_sizes'] = peak_sizes

            # Step 5: Multimodal Fitting
            fit_params, fitted_df, r_squared = fit_multimodal_lognorm(
                diameters_x_data,
                y_psd_data_normalized,
                num_modes,
                peak_sizes,
                selected_weighted_data,
                weighting_scheme=weighting_scheme)
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
            plot_individual_modes(
                diameters_x_data,
                y_psd_data_normalized,
                fitted_df[selected_weighted_data],
                fit_params,
                num_modes,
                selected_weighted_data
            )

            all_results.append(current_results)
        except Exception as exc:
            print(f'Error analysing {file_ID}: {exc}')

    if not all_results:
        print('\nNo files were analysed successfully; no results written.')
        return

    all_results_df = pd.DataFrame(all_results)
    output_path = Path(config.OUTPUT_DIR) / config.OUTPUT_FILENAME
    all_results_df.to_excel(output_path, index=False)
    print(f'\nResults saved to: {output_path}')


if __name__ == "__main__":
    main()
