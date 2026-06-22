"""
Batch-fit every .xlsx in the configured data_dir and build a single, clear
report (figure grid + summary table + HTML) of the results.

Reuses the project's own pipeline (load -> interpolate -> peak-detect -> fit)
with the standard defaults from config.py. Problematic files are skipped and
reported, never aborting the batch. Each sample is labelled with the 4-digit
ID parsed from its filename (2501, 2502, 2503, 2601-2613).
"""
import argparse
import re
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive: never blocks on plt.show()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import lognorm
from scipy.integrate import cumulative_trapezoid, trapezoid

import config
from data_processing import load_data, compute_basic_statistics, find_xlsx_files
from peak_detection import find_distribution_peaks, prepare_distribution_data_interpolation
from fitting import fit_multimodal_lognorm, monomodal_lognorm

warnings.filterwarnings("ignore")

FITTING_PARAM_LIST = [
    "Size distribution Surface weighted [%]",
    "Size distribution Volume weighted [%]",
    "Size distribution Number weighted [%]",
]

_parser = argparse.ArgumentParser(description="Batch PSD fit report.")
_parser.add_argument("--param_fit", type=int, default=config.PARAM_FIT,
                     help="Weighting index: 0=Surface, 1=Volume, 2=Number "
                          f"(default: {config.PARAM_FIT})")
_args = _parser.parse_args()

PARAM_FIT = _args.param_fit
selected_weighted_data = FITTING_PARAM_LIST[PARAM_FIT]
selected_parameter = re.search(r"Size distribution (.*?) \[%\]", selected_weighted_data).group(1)

# Per-weighting output folder so Surface/Volume/Number reports don't clobber each other.
OUT_DIR = Path(config.OUTPUT_DIR) / f"report_{selected_parameter.split()[0].lower()}"
OUT_DIR.mkdir(exist_ok=True)


def sample_label(filename: str) -> str:
    """Pull the 4-digit sample ID (2501/2502/2503/2601-2613) from the name."""
    m = re.search(r"Altris-(\d{4})", filename)
    if m:
        return m.group(1)
    m = re.search(r"\b(25\d{2}|26\d{2})\b", filename)
    return m.group(1) if m else Path(filename).stem


def fit_one(path: str):
    """Run the standard pipeline on one file; return a result dict or raise."""
    data_df = load_data(path, preprocess=True, param_fit=selected_weighted_data,
                        min_diameter=config.MIN_DIAMETER)
    raw = compute_basic_statistics(data_df, selected_parameter)[selected_parameter]

    x, y = prepare_distribution_data_interpolation(
        data_df, selected_weighted_data, interp_points=config.INTERP_POINTS)

    abs_prom = config.PEAK_REL_PROMINENCE * float(np.max(y))
    num_modes, peak_sizes = find_distribution_peaks(
        x, y, prominence=abs_prom, distance=config.PEAK_DISTANCE)
    if num_modes == 0:
        # No peak cleared the threshold; fall back to the global maximum.
        num_modes, peak_sizes = 1, np.array([x[int(np.argmax(y))]])
    # The fitter supports at most trimodal; extra detected peaks are tiny
    # artefacts. Clamp to 3 so the per-mode loop stays in sync with the fit.
    num_modes = min(num_modes, 3)

    weighting = None if config.WEIGHTING_SCHEME == "none" else config.WEIGHTING_SCHEME
    fit_params, fitted_df, r2 = fit_multimodal_lognorm(
        x, y, num_modes, peak_sizes, selected_weighted_data, weighting_scheme=weighting)

    diameters = fitted_df["Particle diameter  [µm]"].values
    modes = []
    for i in range(1, num_modes + 1):
        weight = fit_params[f"Weight{i}"]
        mu = np.log(fit_params[f"Mean{i}"])
        sigma = fit_params[f"Sigma{i}"]
        psd = weight * monomodal_lognorm(diameters, mu, sigma)
        integral = trapezoid(psd, diameters)
        if integral <= 0:
            continue
        norm = psd / integral
        cum = cumulative_trapezoid(norm, diameters, initial=0) * 100
        tmp = pd.DataFrame({
            "Particle diameter  [µm]": diameters,
            f"Size distribution {selected_parameter} [%]": norm * 100,
            f"Undersize {selected_parameter} [%]": cum,
        })
        ms = compute_basic_statistics(tmp, selected_parameter)[selected_parameter]
        modes.append({"weight": weight, "D10": ms["D10"], "D50": ms["D50"],
                      "D90": ms["D90"], "mean": ms["Mean"], "std": ms["Std"],
                      "sigma": sigma})

    modes.sort(key=lambda m: m["D50"])  # smallest mode first for readability
    return {
        "x": x, "y": y, "fit_total": fitted_df[selected_weighted_data].values,
        "diameters": diameters, "fit_params": fit_params, "num_modes": num_modes,
        "r2": r2, "raw": raw, "modes": modes,
    }


def main():
    files = sorted(find_xlsx_files(config.DATA_DIR), key=lambda f: sample_label(f))
    results, skipped = {}, []

    for path in files:
        label = sample_label(path)
        try:
            res = fit_one(path)
            res["filename"] = Path(path).name
            results[label] = res
            print(f"[ok]   {label}: {res['num_modes']} mode(s), "
                  f"D50={res['raw']['D50']:.2f} um, R2={res['r2']:.3f}")
        except Exception as exc:
            skipped.append((label, Path(path).name, str(exc)))
            print(f"[skip] {label}: {exc}")

    if not results:
        print("No files fitted successfully.")
        return

    labels = sorted(results.keys())

    # ---- Figure grid: one fit panel per sample -------------------------------
    n = len(labels)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows))
    axes = np.atleast_1d(axes).ravel()
    mode_colors = ["#1f77b4", "#2ca02c", "#9467bd"]

    for ax, label in zip(axes, labels):
        r = results[label]
        x, diam = r["x"], r["diameters"]
        ax.plot(x, r["y"], color="black", lw=1.6, label="Data")
        ax.plot(diam, r["fit_total"], "r--", lw=1.4, label="Total fit")
        fp = r["fit_params"]
        for i in range(r["num_modes"]):
            comp = fp[f"Weight{i+1}"] * lognorm.pdf(
                diam, s=fp[f"Sigma{i+1}"], loc=0, scale=fp[f"Mean{i+1}"])
            ax.fill_between(diam, comp, alpha=0.25, color=mode_colors[i % 3])
            ax.plot(diam, comp, color=mode_colors[i % 3], lw=1.0,
                    label=f"Mode {i+1} (~{fp[f'Mean{i+1}']:.1f} um)")
        ax.set_xscale("log")
        ax.set_xlim(0.5, 100)
        ax.set_title(f"Sample {label}", fontweight="bold")
        ax.set_xlabel("Particle diameter [um]")
        ax.set_ylabel("Norm. PSD")
        ax.grid(True, which="both", ls=":", alpha=0.5)
        ax.text(0.97, 0.95,
                f"D50={r['raw']['D50']:.1f} um\n$R^2$={r['r2']:.3f}\n{r['num_modes']} mode(s)",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
        ax.legend(fontsize=6.5, loc="upper left", framealpha=0.7)

    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle(f"PSD multimodal lognormal fits - {selected_parameter}",
                 fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    grid_png = OUT_DIR / "psd_fits_grid.png"
    fig.savefig(grid_png, dpi=130)
    plt.close(fig)

    # ---- Summary table -------------------------------------------------------
    rows = []
    for label in labels:
        r = results[label]
        row = {
            "Sample": label,
            "D10": round(r["raw"]["D10"], 2),
            "D50": round(r["raw"]["D50"], 2),
            "D90": round(r["raw"]["D90"], 2),
            "Modes": r["num_modes"],
            "R2": round(r["r2"], 3),
        }
        for i, m in enumerate(r["modes"], start=1):
            row[f"M{i} center [um]"] = round(m["mean"], 2)
            row[f"M{i} weight [%]"] = round(m["weight"] * 100, 1)
        rows.append(row)
    table_df = pd.DataFrame(rows).fillna("")

    # Save the standard Excel output too.
    table_df.to_excel(OUT_DIR / "psd_summary.xlsx", index=False)

    # Render the table as an image for the report.
    fig2, ax2 = plt.subplots(figsize=(min(2 + 1.1 * len(table_df.columns), 20),
                                       0.5 + 0.4 * len(table_df)))
    ax2.axis("off")
    tbl = ax2.table(cellText=table_df.values, colLabels=table_df.columns,
                    cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.4)
    for j, col in enumerate(table_df.columns):
        tbl[0, j].set_facecolor("#34495e")
        tbl[0, j].set_text_props(color="white", fontweight="bold")
    for i in range(1, len(table_df) + 1):
        if i % 2 == 0:
            for j in range(len(table_df.columns)):
                tbl[i, j].set_facecolor("#f0f3f7")
    ax2.set_title(f"PSD fit summary - {selected_parameter}",
                  fontweight="bold", pad=12)
    fig2.tight_layout()
    table_png = OUT_DIR / "psd_summary_table.png"
    fig2.savefig(table_png, dpi=150, bbox_inches="tight")
    plt.close(fig2)

    # ---- HTML report ---------------------------------------------------------
    skipped_html = ""
    if skipped:
        items = "".join(f"<li><b>{lab}</b> ({fn}): {err}</li>"
                        for lab, fn, err in skipped)
        skipped_html = f"<h2>Skipped files</h2><ul>{items}</ul>"
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>PSD Analysis Report</title>
<style>
body{{font-family:Segoe UI,Arial,sans-serif;margin:2rem;color:#222;background:#fafafa}}
h1{{color:#2c3e50}} h2{{color:#34495e;border-bottom:2px solid #ddd;padding-bottom:.3rem}}
img{{max-width:100%;height:auto;border:1px solid #ddd;border-radius:6px;background:#fff}}
.meta{{color:#666;font-size:.9rem}}
</style></head><body>
<h1>Particle Size Distribution - Batch Fit Report</h1>
<p class="meta">Weighting: <b>{selected_parameter}</b> &middot;
Samples fitted: <b>{len(results)}</b> &middot; Skipped: <b>{len(skipped)}</b> &middot;
Source: <code>{config.DATA_DIR}</code></p>
<h2>Summary table</h2>
<img src="psd_summary_table.png" alt="summary table">
<h2>Fitted distributions</h2>
<img src="psd_fits_grid.png" alt="fit grid">
{skipped_html}
</body></html>"""
    (OUT_DIR / "PSD_Report.html").write_text(html, encoding="utf-8")

    print(f"\nFitted {len(results)} samples, skipped {len(skipped)}.")
    print(f"Report written to: {OUT_DIR.resolve()}")
    print("  - PSD_Report.html (open this)")
    print("  - psd_fits_grid.png")
    print("  - psd_summary_table.png")
    print("  - psd_summary.xlsx")


if __name__ == "__main__":
    main()
