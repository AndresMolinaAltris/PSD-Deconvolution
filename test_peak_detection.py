"""
Verification script for the configurable peak-detection thresholds.

For every .xlsx file in Data_Test it runs peak detection twice:
  1. with the default thresholds (current behaviour, no extra kwargs)
  2. with lowered thresholds tuned to catch small / low-amplitude peaks

It prints the detected number of modes and peak diameters for each case and
checks that lowering the thresholds never finds *fewer* peaks than the default.

Run with the project venv:
    .venv/Scripts/python.exe test_peak_detection.py
"""
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from data_processing import load_data
from peak_detection import (
    detect_peaks,
    prepare_distribution_data_interpolation,
    find_distribution_peaks,
)

DATA_DIR = Path(__file__).parent / "Data_Test"
PARAM = "Size distribution Volume weighted [%]"

# Lowered thresholds that make small peaks easier to catch.
SMALL_PEAK_KWARGS = dict(prominence=1e-4, distance=2)


def run_file(path: Path):
    print(f"\n=== {path.name} ===")
    try:
        df = load_data(str(path), preprocess=True, param_fit=PARAM, min_diameter=0.5)
    except Exception as exc:  # noqa: BLE001 - report and skip non-conforming files
        print(f"  SKIPPED (could not load/preprocess): {exc}")
        return None

    # Raw distribution path (detect_peaks): default vs lowered thresholds.
    _, _, n_def, sizes_def = detect_peaks(df.copy(), PARAM)
    _, _, n_low, sizes_low = detect_peaks(df.copy(), PARAM, **SMALL_PEAK_KWARGS)
    print(f"  detect_peaks            default: {n_def} modes {sizes_def}")
    print(f"  detect_peaks            lowered: {n_low} modes {sizes_low}")

    # Interpolated distribution path (find_distribution_peaks), as used by main.py.
    diameters, y = prepare_distribution_data_interpolation(df, PARAM)
    n_def_i, sizes_def_i = find_distribution_peaks(diameters, y)
    n_low_i, sizes_low_i = find_distribution_peaks(diameters, y, **SMALL_PEAK_KWARGS)
    print(f"  find_distribution_peaks default: {n_def_i} modes {sizes_def_i}")
    print(f"  find_distribution_peaks lowered: {n_low_i} modes {sizes_low_i}")

    assert n_low >= n_def, "Lowered thresholds found fewer peaks (detect_peaks)"
    assert n_low_i >= n_def_i, "Lowered thresholds found fewer peaks (find_distribution_peaks)"
    return n_def, n_low, n_def_i, n_low_i


def main():
    files = sorted(DATA_DIR.glob("*.xlsx"))
    if not files:
        print(f"No .xlsx files found in {DATA_DIR}")
        return
    print(f"Testing peak detection on {len(files)} file(s) in {DATA_DIR}")
    extra_found = 0
    for path in files:
        result = run_file(path)
        if result is not None:
            n_def, n_low, n_def_i, n_low_i = result
            extra_found += (n_low - n_def) + (n_low_i - n_def_i)
    print(
        f"\nAll assertions passed. Lowered thresholds surfaced {extra_found} "
        f"extra peak(s) across the dataset that the defaults missed."
    )


if __name__ == "__main__":
    main()
