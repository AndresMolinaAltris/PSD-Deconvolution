# PSD Analysis

Automated particle size distribution (PSD) analysis tool. Loads `.xlsx` measurement files, fits a multimodal lognormal model (a sum of Gaussians in log-diameter) directly to the raw bins, computes D10/D50/D90 statistics per mode **with propagated uncertainties**, and exports results to Excel.

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
| `--min_diameter` | `MIN_DIAMETER` = `0.5` | Lower particle-diameter cutoff (µm). Rows below this are dropped before fitting. Affects measured and per-mode statistics. |
| `--prominence` | `PEAK_REL_PROMINENCE` = `0.02` | Seed-peak prominence, **relative** — a fraction of each file's peak height (so the threshold adapts per file). Lower it to catch smaller modes; raise it to ignore bin-to-bin noise. |
| `--distance` | `PEAK_DISTANCE` = `2` | Minimum separation between detected peaks, in bins. Lower it to resolve peaks that sit close together. |
| `--weighting` | `WEIGHTING_SCHEME` = `none` | Least-squares weighting: `none` (unweighted), `proportional_to_y`, or `sqrt_y`. |
| `--num_modes` | _(auto-detected)_ | Force the number of modes to fit, bypassing seed detection. **This is the main control for deconvolution depth** — use it when a small mode is a shoulder rather than a distinct peak, or to stop over-fitting a single population into several modes. |

`--prominence` and `--distance` only choose the **seeds** for the fit; `--num_modes`
sets the mode count outright and is the parameter to reach for. To make a new value the
default for every run, change the corresponding constant in `config.py`.

### How the fit works (and why it is robust)

Each mode is lognormal in diameter `d`, i.e. a **Gaussian in `u = ln(d)`**. The fit is a
sum of `N` independent, non-negative Gaussians in `u`, fitted directly to the **raw
instrument bins** (no interpolation). This matters:

- **Right coordinate.** Fitting in `u` keeps the model well-conditioned and the residual
  weighting uniform across decades. Every statistic is then analytic: `D50 = exp(µ)`,
  `D_p = exp(µ + σ·z_p)`.
- **No interpolation artefacts.** Cubic interpolation overshoots in the near-zero tail
  and invents spurious large-diameter modes (the fit used to chase one to ~10¹² µm on
  some surface-weighted files). Fitting raw bins removes that entire failure class, and a
  guard drops any mode whose centre lands outside the measured range.
- **Honest uncertainty.** Because `dD/D = dµ + z·dσ`, a small log-space error is a fixed
  *relative* diameter error. The fit covariance is captured and every D-value is reported
  with its relative uncertainty (`Mode*_D50_relerr`) plus a Monte-Carlo CI on the overall
  D50. A poorly-constrained small mode shows up as a large `relerr`, not a false-precision
  number.

### Resolving the small (~2–5 µm) secondary mode

Many distributions contain a small secondary population around 2–5 µm. When it forms a
distinct peak it is seeded and fitted automatically:

```bash
python main.py --files "Recalculation of PWA-Altris-2503-Rep1_....xlsx"
# -> 2 modes: small ~2 µm (high relerr) + main ~28 µm
```

When the small mode is only a **shoulder** (no local maximum, e.g. `2501`, `2613`), seed
detection cannot see it — force the count. Extra modes are seeded in the small-particle
region, so they land on the small population:

```bash
python main.py --files "PWA-Altris-2501-Rep1-Mie_....xlsx" --num_modes 2
# -> small ~1.4 µm + main ~12 µm, R² 0.994 -> 0.997
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
| `D10_raw`, `D50_raw`, `D90_raw` | Measured percentile diameters from the data (independent of the fit — a QC reference) |
| `D10_fit`, `D50_fit`, `D90_fit` | Percentiles of the fitted mixture |
| `num_modes` | Number of fitted (in-range) modes |
| `R_squared` | Goodness of fit (log-diameter space) |
| `Mode1_Weight` | Mode mass fraction |
| `Mode1_Dg`, `Mode1_sigma_g` | Mode geometric mean (`exp µ`) and geometric SD (`exp σ`) |
| `Mode1_D10/D50/D90` | Per-mode percentile diameters (analytic) |
| `Mode1_D50_relerr` | Relative uncertainty of the mode's D50 (`dD/D`) |
| `Mode1_Mean`, `Mode1_Std` | Arithmetic mean/SD of the mode (linear diameter) |

(per-mode columns repeat for Mode2, Mode3, … up to the fitted count.)

## Project structure

```
main.py                  # Entry point
batch_report.py          # Batch-fit a folder -> figure grid + table + HTML report
config.template.json     # Template config committed to git (copy to config.json)
config.json              # Your local config with real paths (git-ignored)
config.py                # Loads config.json paths + holds default fitting constants
data_processing.py       # File loading and measured statistics
peak_detection.py        # Raw-bin extraction and log-space mode seeding
fitting.py               # Gaussian-in-log-diameter fit, analytic stats, uncertainties
visualization.py         # Plotting
test_peak_detection.py   # Verifies seed detection + no out-of-range modes on Data_Test
requirements.txt         # Python dependencies
```
