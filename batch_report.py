"""
Batch-fit every .xlsx in the configured data_dir and build a single, clear
report (figure grid + summary table + HTML) of the results.

Reuses the project pipeline (raw bins -> seed detection -> Gaussian-in-log-d
fit) with the standard defaults from config.py. Problematic files are skipped
and reported, never aborting the batch. Each sample is labelled with the
4-digit ID parsed from its filename.

Usage:  python batch_report.py [--param_fit 0|1|2]   (0=Surface, 1=Volume, 2=Number)
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

import config
from data_processing import load_data, compute_basic_statistics, find_xlsx_files
from peak_detection import psd_arrays, detect_peaks
from fitting import fit_psd
from visualization import plot_individual_modes

warnings.filterwarnings("ignore")

FITTING_PARAM_LIST = [
    "Size distribution Surface weighted [%]",
    "Size distribution Volume weighted [%]",
    "Size distribution Number weighted [%]",
]

_parser = argparse.ArgumentParser(description="Batch PSD fit report.")
_parser.add_argument("--param_fit", type=int, default=config.PARAM_FIT,
                     help=f"0=Surface, 1=Volume, 2=Number (default: {config.PARAM_FIT})")
_args = _parser.parse_args()

PARAM_FIT = _args.param_fit
selected_weighted_data = FITTING_PARAM_LIST[PARAM_FIT]
selected_parameter = re.search(r"Size distribution (.*?) \[%\]", selected_weighted_data).group(1)

OUT_DIR = Path(config.OUTPUT_DIR) / f"report_{selected_parameter.split()[0].lower()}"
OUT_DIR.mkdir(exist_ok=True)


def sample_label(filename: str) -> str:
    m = re.search(r"Altris-(\d{4})", filename)
    if m:
        return m.group(1)
    m = re.search(r"\b(25\d{2}|26\d{2})\b", filename)
    return m.group(1) if m else Path(filename).stem


def fit_one(path: str):
    data_df = load_data(path, preprocess=True, param_fit=selected_weighted_data,
                        min_diameter=config.MIN_DIAMETER)
    raw = compute_basic_statistics(data_df, selected_parameter)[selected_parameter]
    d, y = psd_arrays(data_df, selected_weighted_data)
    seeds = detect_peaks(d, y, prominence_rel=config.PEAK_REL_PROMINENCE,
                         distance=config.PEAK_DISTANCE)
    weighting = None if config.WEIGHTING_SCHEME == "none" else config.WEIGHTING_SCHEME
    fit = fit_psd(d, y, len(seeds), seeds, weighting=weighting,
                  default_sigma=config.DEFAULT_SIGMA, mc_samples=config.MC_SAMPLES)
    fit["raw"] = raw
    return fit


def main():
    files = sorted(find_xlsx_files(config.DATA_DIR), key=sample_label)
    results, skipped = {}, []
    for path in files:
        label = sample_label(path)
        try:
            results[label] = fit_one(path)
            r = results[label]
            print(f"[ok]   {label}: {r['n_modes']} mode(s), "
                  f"D50_fit={r['overall'][50]:.2f} µm, R²={r['r_squared']:.3f}")
        except Exception as exc:
            skipped.append((label, Path(path).name, str(exc)))
            print(f"[skip] {label}: {exc}")

    if not results:
        print("No files fitted successfully.")
        return

    labels = sorted(results.keys())

    # ---- Figure grid: one fit panel per sample (reuses plot_individual_modes) -
    n, ncols = len(labels), 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 3.4 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, label in zip(axes, labels):
        r = results[label]
        plot_individual_modes(r, selected_weighted_data, ax=ax)
        ax.set_title(f"Sample {label}", fontweight="bold")
        ov = r["overall"]
        ax.text(0.97, 0.95,
                f"D50={ov[50]:.1f} µm\n$R^2$={r['r_squared']:.3f}\n{r['n_modes']} mode(s)",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
        ax.legend(fontsize=6.5, loc="upper left", framealpha=0.7)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(f"PSD Gaussian-in-log-d fits — {selected_parameter}",
                 fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(OUT_DIR / "psd_fits_grid.png", dpi=130)
    plt.close(fig)

    # ---- Summary table -------------------------------------------------------
    rows = []
    for label in labels:
        r = results[label]
        ov = r["overall"]
        row = {
            "Sample": label,
            "D10_meas": round(r["raw"]["D10"], 2),
            "D50_meas": round(r["raw"]["D50"], 2),
            "D90_meas": round(r["raw"]["D90"], 2),
            "D50_fit": round(ov[50], 2),
            "Modes": r["n_modes"],
            "R2": round(r["r_squared"], 3),
        }
        for i, m in enumerate(r["modes"], start=1):
            row[f"M{i} Dg [µm]"] = round(m["Dg"], 2)
            row[f"M{i} wt [%]"] = round(m["weight"] * 100, 1)
            row[f"M{i} D50±%"] = round(m["D50_rel_err"] * 100, 1)
        rows.append(row)
    table_df = pd.DataFrame(rows).fillna("")
    table_df.to_excel(OUT_DIR / "psd_summary.xlsx", index=False)

    fig2, ax2 = plt.subplots(figsize=(min(2 + 1.05 * len(table_df.columns), 22),
                                       0.5 + 0.4 * len(table_df)))
    ax2.axis("off")
    tbl = ax2.table(cellText=table_df.values, colLabels=table_df.columns,
                    cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1, 1.4)
    for j in range(len(table_df.columns)):
        tbl[0, j].set_facecolor("#34495e")
        tbl[0, j].set_text_props(color="white", fontweight="bold")
    for i in range(1, len(table_df) + 1):
        if i % 2 == 0:
            for j in range(len(table_df.columns)):
                tbl[i, j].set_facecolor("#f0f3f7")
    ax2.set_title(f"PSD fit summary — {selected_parameter} "
                  f"(M# Dg = mode geo-mean, D50±% = relative uncertainty)",
                  fontweight="bold", pad=12)
    fig2.tight_layout()
    fig2.savefig(OUT_DIR / "psd_summary_table.png", dpi=150, bbox_inches="tight")
    plt.close(fig2)

    # ---- HTML report ---------------------------------------------------------
    skipped_html = ""
    if skipped:
        items = "".join(f"<li><b>{lab}</b> ({fn}): {err}</li>" for lab, fn, err in skipped)
        skipped_html = f"<h2>Skipped files</h2><ul>{items}</ul>"
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>PSD Analysis Report</title>
<style>
body{{font-family:Segoe UI,Arial,sans-serif;margin:2rem;color:#222;background:#fafafa}}
h1{{color:#2c3e50}} h2{{color:#34495e;border-bottom:2px solid #ddd;padding-bottom:.3rem}}
img{{max-width:100%;height:auto;border:1px solid #ddd;border-radius:6px;background:#fff}}
.meta{{color:#666;font-size:.9rem}}
</style></head><body>
<h1>Particle Size Distribution — Batch Fit Report</h1>
<p class="meta">Weighting: <b>{selected_parameter}</b> &middot;
Samples fitted: <b>{len(results)}</b> &middot; Skipped: <b>{len(skipped)}</b> &middot;
Model: sum of Gaussians in ln(d) on raw bins &middot; Source: <code>{config.DATA_DIR}</code></p>
<h2>Summary table</h2>
<img src="psd_summary_table.png" alt="summary table">
<h2>Fitted distributions</h2>
<img src="psd_fits_grid.png" alt="fit grid">
{skipped_html}
</body></html>"""
    (OUT_DIR / "PSD_Report.html").write_text(html, encoding="utf-8")

    print(f"\nFitted {len(results)} samples, skipped {len(skipped)}.")
    print(f"Report written to: {OUT_DIR.resolve()}  (open PSD_Report.html)")


if __name__ == "__main__":
    main()
