# PSD Analysis

Automated particle size distribution (PSD) analysis tool. Loads `.xlsx` measurement files, fits multimodal lognormal distributions, computes D10/D50/D90 statistics per mode, and exports results to Excel.

## Setup

**Requirements:** Python 3.8+

```bash
pip install -r requirements.txt
```

**First-time configuration:**

```bash
cp config.template.json config.json
```

Then edit `config.json` with your local paths (this file is git-ignored and never committed).

## Configuration

All paths are set in **`config.json`** — edit this file before running:

```json
{
  "data_dir": "./Data_Test",
  "output_dir": ".",
  "output_filename": "results.xlsx"
}
```

| Key | Description |
|-----|-------------|
| `data_dir` | Folder containing your `.xlsx` data files |
| `output_dir` | Folder where `results.xlsx` is saved |
| `output_filename` | Name of the output Excel file |

## Usage

1. **Put your data files** in the folder specified by `data_dir` in `config.json`.

2. **Edit `config.json`** to point `data_dir` at your data folder.

3. **Run the analysis:**

   ```bash
   python main.py
   ```

   Optionally select a weighting mode with `--param_fit`:

   ```bash
   python main.py --param_fit 1   # 0 = Surface, 1 = Volume (default), 2 = Number
   ```

   Optionally limit analysis to specific files with `--files`:

   ```bash
   python main.py --files "sample1.xlsx"
   python main.py --files "sample1.xlsx" "sample2.xlsx"
   ```

   If `--files` is omitted, all `.xlsx` files in `data_dir` are processed.

## Examples

Run with default settings (Volume weighted, paths from `config.json`):
```bash
python main.py
```

Run with Surface weighted distribution:
```bash
python main.py --param_fit 0
```

Run with Volume weighted distribution (explicit):
```bash
python main.py --param_fit 1
```

Run with Number weighted distribution:
```bash
python main.py --param_fit 2
```

Analyse a single file:
```bash
python main.py --files "UKP240701 QC Rep1.xlsx"
```

Analyse multiple specific files:
```bash
python main.py --files "sample1.xlsx" "sample2.xlsx"
```

Analyse specific files with a chosen weighting:
```bash
python main.py --files "sample1.xlsx" --param_fit 0
```

Check available options:
```bash
python main.py --help
```

4. **Outputs:**
   - Console: per-file statistics printed as files are processed
   - `results.xlsx` (in `output_dir`): summary table with D10/D50/D90 and per-mode statistics for all files
   - Plot windows: fitted distribution and individual modes for each file

## Expected Excel format

The tool reads files exported from Malvern Mastersizer or equivalent instruments. Each `.xlsx` file must have:
- Data in columns **F through V**
- Header rows in rows **5–7** (rows 1–4 are skipped)
- Columns: `Particle diameter  [µm]`, `Size distribution <X> weighted [%]`, `Undersize <X> weighted [%]`

## Output columns

| Column | Description |
|--------|-------------|
| `Filename` | Source file name |
| `Parameter` | Weighting used (e.g. `Volume weighted`) |
| `D10_raw`, `D50_raw`, `D90_raw` | Percentile diameters before deconvolution |
| `num_modes` | Number of detected peaks |
| `R_squared` | Goodness of fit |
| `Mode1_Weight`, `Mode1_D10/D50/D90/Mean/Std` | Per-mode statistics (repeated for Mode2, Mode3 if present) |

## Project structure

```
main.py                  # Entry point
config.template.json     # Template config committed to git (copy to config.json)
config.json              # Your local config with real paths (git-ignored)
config.py                # Loads config.json and exposes settings to the code
data_processing.py       # File loading and basic statistics
peak_detection.py        # Peak detection and data interpolation
fitting.py               # Lognormal distribution fitting
visualization.py         # Plotting
requirements.txt         # Python dependencies
```
