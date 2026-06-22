import argparse
import re
from pathlib import Path

import pandas as pd

import config
from data_processing import load_data, compute_basic_statistics, find_xlsx_files
from peak_detection import psd_arrays, detect_peaks
from fitting import fit_psd
from visualization import plot_individual_modes

FITTING_PARAM_LIST = [
    "Size distribution Surface weighted [%]",
    "Size distribution Volume weighted [%]",
    "Size distribution Number weighted [%]",
]


def build_parser():
    p = argparse.ArgumentParser(
        description="PSD analysis. Fits a sum of Gaussians in log-diameter to the "
                    "raw bins and reports per-mode D-values with uncertainties. "
                    "Optional arguments default to config.py; omit them for the "
                    "standard fit.")
    p.add_argument("--param_fit", type=int, default=config.PARAM_FIT,
                   help=f"Weighting: 0=Surface, 1=Volume, 2=Number (default: {config.PARAM_FIT})")
    p.add_argument("--files", nargs="+", default=None,
                   help="Filenames to analyse. If omitted, all .xlsx in data_dir.")
    p.add_argument("--min_diameter", type=float, default=config.MIN_DIAMETER,
                   help=f"Lower diameter cutoff in µm (default: {config.MIN_DIAMETER})")
    p.add_argument("--prominence", type=float, default=config.PEAK_REL_PROMINENCE,
                   help="Peak prominence as a fraction of peak height, for seed "
                        f"detection (default: {config.PEAK_REL_PROMINENCE})")
    p.add_argument("--distance", type=int, default=config.PEAK_DISTANCE,
                   help=f"Min separation between peaks, in bins (default: {config.PEAK_DISTANCE})")
    p.add_argument("--weighting", choices=["proportional_to_y", "sqrt_y", "none"],
                   default=config.WEIGHTING_SCHEME,
                   help=f"Least-squares weighting (default: {config.WEIGHTING_SCHEME})")
    p.add_argument("--num_modes", type=int, default=None,
                   help="Force the number of modes to fit, overriding seed detection. "
                        "Use it when a small mode is only a shoulder (no local peak).")
    p.add_argument("--no_plots", action="store_true", help="Skip plot windows.")
    return p


def analyse_file(path, swd, parameter, args, weighting):
    """Fit one file; return a flat result dict for the summary table."""
    file_id = Path(path).name
    res = {"Filename": file_id, "Parameter": parameter}

    data_df = load_data(path, preprocess=True, param_fit=swd,
                        min_diameter=args.min_diameter)

    # Measured (pre-fit) percentiles, independent of the fit — a QC reference.
    raw = compute_basic_statistics(data_df, parameter)[parameter]
    res["D10_raw"] = round(raw["D10"], 3)
    res["D50_raw"] = round(raw["D50"], 3)
    res["D90_raw"] = round(raw["D90"], 3)

    d, y = psd_arrays(data_df, swd)
    seeds = detect_peaks(d, y, prominence_rel=args.prominence, distance=args.distance)
    num_modes = args.num_modes if args.num_modes is not None else len(seeds)
    num_modes = max(1, num_modes)

    fit = fit_psd(d, y, num_modes, seeds, weighting=weighting,
                  default_sigma=config.DEFAULT_SIGMA, mc_samples=config.MC_SAMPLES)

    res["num_modes"] = fit["n_modes"]
    res["R_squared"] = round(fit["r_squared"], 4)
    ov = fit["overall"]
    res["D10_fit"] = round(ov[10], 3)
    res["D50_fit"] = round(ov[50], 3)
    res["D90_fit"] = round(ov[90], 3)

    print(f"\nFile analysed: {file_id}  ({parameter})")
    print(f"  measured  D10/D50/D90 = {raw['D10']:.2f}/{raw['D50']:.2f}/{raw['D90']:.2f} µm")
    ci = ov.get("D50_ci")
    ci_txt = f"  [68% CI {ci[0]:.2f}-{ci[1]:.2f}]" if ci else ""
    print(f"  fitted    D10/D50/D90 = {ov[10]:.2f}/{ov[50]:.2f}/{ov[90]:.2f} µm{ci_txt}"
          f"   R²={fit['r_squared']:.4f}, {fit['n_modes']} mode(s)")
    if fit["dropped"]:
        print(f"  (dropped {len(fit['dropped'])} out-of-range mode(s) at "
              f"{[round(x, 1) for x in fit['dropped']]} µm)")

    for i, m in enumerate(fit["modes"], start=1):
        res[f"Mode{i}_Weight"] = round(m["weight"], 3)
        res[f"Mode{i}_Dg"] = round(m["Dg"], 3)
        res[f"Mode{i}_sigma_g"] = round(m["sigma_g"], 3)
        res[f"Mode{i}_D10"] = round(m["D10"], 3)
        res[f"Mode{i}_D50"] = round(m["D50"], 3)
        res[f"Mode{i}_D90"] = round(m["D90"], 3)
        res[f"Mode{i}_D50_relerr"] = round(m["D50_rel_err"], 4)
        res[f"Mode{i}_Mean"] = round(m["Mean"], 3)
        res[f"Mode{i}_Std"] = round(m["Std"], 3)
        print(f"    Mode {i}: weight {m['weight']*100:5.1f}%  Dg={m['Dg']:7.2f} µm  "
              f"σg={m['sigma_g']:.2f}  D50={m['D50']:.2f}±{m['D50_rel_err']*100:.1f}%")

    if not args.no_plots:
        plot_individual_modes(fit, swd)

    return res


def main():
    args = build_parser().parse_args()
    weighting = None if args.weighting == "none" else args.weighting
    swd = FITTING_PARAM_LIST[args.param_fit]
    parameter = re.search(r"Size distribution (.*?) \[%\]", swd).group(1)
    parent_dir = config.DATA_DIR

    if args.files:
        xlsx_files = []
        for f in args.files:
            path = Path(parent_dir) / f
            if path.exists():
                xlsx_files.append(str(path))
            else:
                print(f"Warning: file not found, skipping: {path}")
    else:
        xlsx_files = find_xlsx_files(parent_dir)

    if not xlsx_files:
        print(f"No .xlsx files to analyse in {parent_dir}. Nothing to do.")
        return

    all_results = []
    for path in xlsx_files:
        try:
            all_results.append(analyse_file(path, swd, parameter, args, weighting))
        except Exception as exc:
            print(f"Error analysing {Path(path).name}: {exc}")

    if not all_results:
        print("\nNo files were analysed successfully; no results written.")
        return

    out_path = Path(config.OUTPUT_DIR) / config.OUTPUT_FILENAME
    pd.DataFrame(all_results).to_excel(out_path, index=False)
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
