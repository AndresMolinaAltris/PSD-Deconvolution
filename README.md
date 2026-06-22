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

`config.json` holds only **paths**. The default **fitting** values (diameter cutoff, peak
thresholds, weighting scheme) live as named constants in `config.py` — see
[Controlling the fit](#controlling-the-fit).

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

   If `--files` is omitted, all `.xlsx` files in `data_dir` are processed. Files that
   don't exist are skipped with a warning, and a failure on one file no longer aborts the
   rest of the batch.

   The fit itself can be tuned with the optional arguments described in
   [Controlling the fit](#controlling-the-fit) below. **Omitting them reproduces the
   standard fit**, so you only pass them when a particular sample needs different
   treatment.

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

Force a fit with a set number of modes (useful when a small mode is only a shoulder):
```bash
python main.py --files "sample1.xlsx" --num_modes 2
```

Reproduce the legacy peak-following fit (weight points by their height):
```bash
python main.py --files "sample1.xlsx" --weighting proportional_to_y
```

Check available options:
```bash
python main.py --help
```

4. **Outputs:**
   - Console: per-file statistics printed as files are processed
   - `results.xlsx` (in `output_dir`): summary table with D10/D50/D90 and per-mode statistics for all files
   - Plot windows: fitted distribution and individual modes for each file

## Controlling the fit

The values that shape the fit have sensible defaults that work well across most
distributions, defined as constants in **`config.py`**. Each can be overridden per run
with a command-line argument — **omit the argument and the default from `config.py` is
used**, so the standard fit is reproduced exactly.

| Argument | Default (`config.py`) | What it does |
|----------|-----------------------|--------------|
| `--min_diameter` | `MIN_DIAMETER` = `0.5` | Lower particle-diameter cutoff (µm). Rows below this are dropped before fitting. Affects raw and per-mode statistics. |
| `--prominence` | `PEAK_REL_PROMINENCE` = `0.001` | Peak prominence for detection, **relative** — a fraction of each file's peak height (so the threshold adapts per file). Lower it to catch smaller modes; raise it to ignore minor bumps. |
| `--distance` | `PEAK_DISTANCE` = `20` | Minimum separation (in samples) between detected peaks. Lower it to resolve peaks that sit close together. |
| `--weighting` | `WEIGHTING_SCHEME` = `none` | Least-squares weighting: `none` (unweighted), `proportional_to_y`, or `sqrt_y`. |
| `--num_modes` | _(auto-detected)_ | Force the number of lognormal modes to fit (`1`–`3`), bypassing automatic peak detection. Useful when a small mode is a shoulder rather than a distinct peak. |

`--prominence` and `--distance` decide **how many modes** are detected; `--num_modes`
lets you override that count outright. To make a new value the default for every run,
change the corresponding constant in `config.py` instead of passing the flag each time.

### Resolving the small (~2–5 µm) volume-weighted mode

The volume-weighted distributions in this dataset contain a small secondary population
around 2–5 µm. Because it carries very little *volume*, it is easy to miss. Two default
choices are tuned to catch it robustly:

- **Unweighted fit (`--weighting none`).** `proportional_to_y` weights each point by its
  height, which suppresses the low small mode — with it, the fit puts both modes on the
  dominant peak and the small mode is never recovered. An unweighted fit gives every
  point equal say and resolves the small mode (and gives equal-or-better R² on the
  bundled samples).
- **Relative prominence.** A threshold of ~0.1 % of the peak height catches the genuine
  small peaks while rejecting the tiny large-diameter interpolation artefacts.

With these defaults, the small mode is resolved automatically on most files:

```bash
python main.py --files "UKP240701 QC Rep1.xlsx"
# -> 2 modes: small ~5 µm + main ~29 µm
```

When the small mode is only a **shoulder** (no distinct local maximum, e.g.
`PWA-Altris-2501-Rep2`), peak detection cannot see it — force the mode count instead.
The extra mode is seeded in the small-particle region, so it lands on the small
population:

```bash
python main.py --files "PWA-Altris-2501-Rep2_....xlsx" --num_modes 2
# -> small ~3.4 µm + main ~16 µm

# A sample with small + shoulder + main needs three modes:
python main.py --files "PSD_PC-4-1-1 0.8 mol Pressed-15-t.xlsx" --num_modes 3
# -> small ~4 µm, R² ~0.999
```

You can sanity-check peak detection against the bundled samples with:

```bash
python test_peak_detection.py
```

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
config.py                # Loads config.json paths + holds default fitting constants
data_processing.py       # File loading and basic statistics
peak_detection.py        # Peak detection and data interpolation
fitting.py               # Lognormal distribution fitting
visualization.py         # Plotting
test_peak_detection.py   # Verifies peak-detection thresholds against Data_Test samples
requirements.txt         # Python dependencies
```
