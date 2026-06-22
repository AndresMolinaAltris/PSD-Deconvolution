"""
Verification script for mode-seed detection (log-diameter space, raw bins).

For every .xlsx in Data_Test it detects seed peaks twice:
  1. with the default prominence (config.PEAK_REL_PROMINENCE)
  2. with a lowered prominence tuned to catch small / low-amplitude modes

It prints the detected peak diameters for each case and checks that lowering
the threshold never finds *fewer* peaks than the default. It also runs a fit
with the detected seeds and asserts no mode runs outside the measured range.

Run with the project venv:
    .venv/Scripts/python.exe test_peak_detection.py
"""
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import config
from data_processing import load_data
from peak_detection import psd_arrays, detect_peaks
from fitting import fit_psd

DATA_DIR = Path(config.DATA_DIR)
PARAM = "Size distribution Volume weighted [%]"
LOW_PROMINENCE = 0.002  # lowered threshold to surface small modes


def run_file(path: Path):
    print(f"\n=== {path.name} ===")
    try:
        df = load_data(str(path), preprocess=True, param_fit=PARAM, min_diameter=0.5)
    except Exception as exc:  # report and skip non-conforming files
        print(f"  SKIPPED (could not load/preprocess): {exc}")
        return None

    d, y = psd_arrays(df, PARAM)
    seeds_def = detect_peaks(d, y, prominence_rel=config.PEAK_REL_PROMINENCE,
                             distance=config.PEAK_DISTANCE)
    seeds_low = detect_peaks(d, y, prominence_rel=LOW_PROMINENCE,
                             distance=config.PEAK_DISTANCE)
    print(f"  default prominence: {len(seeds_def)} seed(s) {[round(s, 2) for s in seeds_def]}")
    print(f"  lowered prominence: {len(seeds_low)} seed(s) {[round(s, 2) for s in seeds_low]}")

    fit = fit_psd(d, y, len(seeds_def), seeds_def, mc_samples=0)
    print(f"  fit: {fit['n_modes']} mode(s), R²={fit['r_squared']:.4f}, "
          f"dropped(out-of-range)={len(fit['dropped'])}")

    assert len(seeds_low) >= len(seeds_def), "Lowered prominence found fewer peaks"
    assert not fit["dropped"], f"Fit produced out-of-range mode(s): {fit['dropped']}"
    return len(seeds_def), len(seeds_low)


def main():
    files = sorted(DATA_DIR.glob("*.xlsx"))
    if not files:
        print(f"No .xlsx files found in {DATA_DIR}")
        return
    print(f"Testing seed detection on {len(files)} file(s) in {DATA_DIR}")
    extra = 0
    for path in files:
        res = run_file(path)
        if res is not None:
            n_def, n_low = res
            extra += n_low - n_def
    print(f"\nAll assertions passed. Lowered prominence surfaced {extra} extra "
          f"seed(s) across the dataset that the default missed.")


if __name__ == "__main__":
    main()
