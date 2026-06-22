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
# The fit is a sum of Gaussians in log-diameter (see fitting.py), run on the
# raw instrument bins. Peak detection only *seeds* the fit; the number of modes
# is set explicitly with --num_modes when the default detection is not wanted.
#   * PEAK_REL_PROMINENCE is RELATIVE (a fraction of the distribution's peak
#     height), so the threshold adapts per file. ~2% catches genuine secondary
#     modes on the raw bins while ignoring bin-to-bin noise. A small mode that
#     is only a shoulder has no local maximum: force it with --num_modes.
#   * WEIGHTING_SCHEME = 'none' (unweighted). In log-space an unweighted fit
#     gives every decade equal say; 'sqrt_y' can emphasise low-amplitude modes.
MIN_DIAMETER        = 0.5     # lower particle-diameter cutoff [µm]
PEAK_REL_PROMINENCE = 0.02    # peak prominence as a fraction of peak height (raw bins)
PEAK_DISTANCE       = 2       # peak min separation (bins)
WEIGHTING_SCHEME    = 'none'  # 'proportional_to_y' | 'sqrt_y' | 'none' (unweighted)
DEFAULT_SIGMA       = 0.35    # initial Gaussian width (log-space) per mode
MC_SAMPLES          = 300     # Monte-Carlo draws for D-value uncertainty (0 disables)
