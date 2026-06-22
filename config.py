import json
from pathlib import Path

_config_path = Path(__file__).parent / 'config.json'
with open(_config_path) as f:
    _cfg = json.load(f)

# Path to the folder containing your .xlsx data files (set in config.json)
DATA_DIR = _cfg['data_dir']

# Directory where results.xlsx is saved (set in config.json)
OUTPUT_DIR = _cfg['output_dir']

# Name of the output Excel file (set in config.json)
OUTPUT_FILENAME = _cfg['output_filename']

# Default weighting for --param_fit argument:
#   0 = Surface weighted
#   1 = Volume weighted
#   2 = Number weighted
PARAM_FIT = 1

# --- Fitting defaults (values that work well across most distributions) ---
# These back the optional command-line arguments in main.py. Omitting a CLI
# argument falls back to the value here.
#
# The defaults are tuned to resolve the small (~2-5 µm) secondary mode that is
# present in the volume-weighted data in a robust way:
#   * PEAK_REL_PROMINENCE is RELATIVE (a fraction of the distribution's peak
#     height), so the threshold adapts to each file. ~0.1% catches the small
#     real modes while rejecting the tiny large-diameter interpolation artefacts.
#   * WEIGHTING_SCHEME = 'none' (unweighted). 'proportional_to_y' weights points
#     by their height and therefore suppresses the low-amplitude small mode; an
#     unweighted fit gives every point equal say and resolves it (and yields
#     equal-or-better R² on the bundled samples).
MIN_DIAMETER        = 0.5     # lower particle-diameter cutoff [µm]
PEAK_REL_PROMINENCE = 0.001   # find_peaks prominence as a fraction of peak height
PEAK_DISTANCE       = 20      # find_peaks min separation (samples)
WEIGHTING_SCHEME    = 'none'  # 'proportional_to_y' | 'sqrt_y' | 'none' (unweighted)
INTERP_POINTS       = 1000    # interpolation grid size
