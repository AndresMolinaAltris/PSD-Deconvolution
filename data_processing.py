import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import os

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

def compute_percentiles(diameters, cum_weights, percentiles=[10, 50, 90]):
    """
    Compute percentiles (e.g., D10, D50, D90) from cumulative distribution.
    
    Parameters:
    - diameters: Array of particle diameters.
    - cum_weights: Array of cumulative weights.
    - percentiles: List of percentiles to compute (default: [10, 50, 90]).
    
    Returns:
    - List of diameters corresponding to the specified percentiles.
    """
    interp_func = interp1d(cum_weights, diameters, bounds_error=False, fill_value=(diameters.min(), diameters.max()))
    return [float(interp_func(p)) for p in percentiles]

# def load_data(filename, preprocess=False, param_fit=None):
#     """
#     Load and optionally preprocess Excel data.
    
#     Parameters:
#     - filename: Path to the Excel file or directory.
#     - preprocess: If True, preprocess the data (sort, clean, etc.).
#     - param_fit: Column name for the fitting parameter (used in preprocessing).
    
#     Returns:
#     - If filename is a directory, returns a list of Excel file paths.
#     - If filename is a file and preprocess=False, returns raw DataFrame.
#     - If preprocess=True, returns preprocessed DataFrame.
#     """
#     if isinstance(filename, str) and os.path.isdir(filename):
#         return find_xlsx_files(filename)
    
#     headers = pd.read_excel(filename, usecols="F:V", nrows=3, skiprows=4)
#     new_cols = headers.astype(str).apply(lambda x: ' '.join(x), axis=0)
#     data_df = pd.read_excel(filename, usecols="F:V", skiprows=5, decimal=".")
#     data_df.columns = new_cols
#     data_df.columns = data_df.columns.str.replace("nan", "", regex=False).str.strip()
    
#     if preprocess:
#         data_df = data_df.drop(index=[0,1]).reset_index(drop=True)
#         data_df = data_df.sort_values('Particle diameter  [µm]')
#         data_df = data_df.dropna(subset=['Particle diameter  [µm]', param_fit.replace('Size distribution', 'Undersize'), param_fit])
#     return data_df

def load_data(filename, preprocess=False, param_fit=None, min_diameter=None):
    """
    Load and optionally preprocess Excel data.
    
    Parameters:
    - filename: Path to the Excel file or directory.
    - preprocess: If True, preprocess the data (sort, clean, etc.).
    - param_fit: Column name for the fitting parameter (used in preprocessing).
    - min_diameter: Optional minimum particle diameter threshold. If provided, rows with 'Particle diameter  [µm]' below this value are removed during preprocessing.
    
    Returns:
    - If filename is a directory, returns a list of Excel file paths.
    - If filename is a file and preprocess=False, returns raw DataFrame.
    - If preprocess=True, returns preprocessed DataFrame.
    """
    if isinstance(filename, str) and os.path.isdir(filename):
        return find_xlsx_files(filename)
    
    headers = pd.read_excel(filename, usecols="F:V", nrows=3, skiprows=4)
    new_cols = headers.astype(str).apply(lambda x: ' '.join(x), axis=0)
    data_df = pd.read_excel(filename, usecols="F:V", skiprows=5, decimal=".")
    data_df.columns = new_cols
    data_df.columns = data_df.columns.str.replace("nan", "", regex=False).str.strip()
    
    if preprocess:
        data_df = data_df.drop(index=[0,1]).reset_index(drop=True)
        data_df = data_df.sort_values('Particle diameter  [µm]')
        data_df = data_df.dropna(subset=['Particle diameter  [µm]', param_fit.replace('Size distribution', 'Undersize'), param_fit])
        if min_diameter is not None:
            data_df = data_df[data_df['Particle diameter  [µm]'] >= min_diameter]
    return data_df

def compute_basic_statistics(data_df, weightings):
    """
    Compute basic statistics (D10, D50, D90, mean, std) for the data.
    
    Parameters:
    - data_df: Preprocessed DataFrame.
    - weightings: List of weightings (e.g., ['Volume weighted']).
    - fitting_param: Column name for the fitting parameter.
    
    Returns:
    - Dictionary with statistics for each weighting.
    """
    results = {}
    #for weight in weightings:
    cum_col = f'Undersize {weightings} [%]'
    diff_col = f'Size distribution {weightings} [%]'
    d10, d50, d90 = compute_percentiles(data_df['Particle diameter  [µm]'], data_df[cum_col])
    mean = np.average(data_df['Particle diameter  [µm]'], weights=data_df[diff_col])
    std = np.sqrt(np.average((data_df['Particle diameter  [µm]'] - mean)**2, weights=data_df[diff_col]))
    results[weightings] = {'D10': d10, 'D50': d50, 'D90': d90, 'Mean': mean, 'Std': std}
    
    return results
    #return results