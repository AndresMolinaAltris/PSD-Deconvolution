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
