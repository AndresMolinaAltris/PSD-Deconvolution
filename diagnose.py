"""
Diagnostic harness (read-only w.r.t. the core modules).

Compares the CURRENT pipeline (lognormal-in-x fit on a cubic-interpolated
log grid) against a PROTOTYPE (sum-of-Gaussians in u = ln(d), fit on the raw
instrument bins, independent non-negative amplitudes, pcov captured).

Goal: produce evidence on three questions raised in review:
  1. Does cubic interpolation invent spurious peaks?  (-> spurious modes)
  2. Is the log-space Gaussian fit better-conditioned / does it give usable
     parameter uncertainties and analytic D-values?
  3. How sensitive are D-values to small parameter errors on a log axis?

Run:  python diagnose.py
Outputs: report_diag/diagnosis.png  + console table.
"""
import re
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import norm
from scipy.special import erf

import config
from data_processing import load_data, find_xlsx_files
from peak_detection import (prepare_distribution_data_interpolation,
                            find_distribution_peaks)
from fitting import fit_multimodal_lognorm, monomodal_lognorm

warnings.filterwarnings("ignore")

PARAMS = {0: "Size distribution Surface weighted [%]",
          1: "Size distribution Volume weighted [%]",
          2: "Size distribution Number weighted [%]"}

# Representative cases: (label, param_fit). 2502-surface is the spurious-mode
# case; 2503/2601 volume are clean/3-mode; 2501 volume is small-mode+main.
CASES = [("2502", 0), ("2503", 1), ("2601", 1), ("2501", 1)]

OUT = Path("report_diag")
OUT.mkdir(exist_ok=True)


def label_of(fn):
    m = re.search(r"Altris-(\d{4})", fn)
    return m.group(1) if m else Path(fn).stem


def path_for(label):
    return next(f for f in find_xlsx_files(config.DATA_DIR) if label_of(f) == label)


# ---- prototype: sum of Gaussians in u = ln(d) -------------------------------
def gauss_sum(u, *p):
    """p = [A1, m1, s1, A2, m2, s2, ...] ; A_i >= 0 independent amplitudes."""
    y = np.zeros_like(u, dtype=float)
    for i in range(0, len(p), 3):
        A, m, s = p[i], p[i + 1], p[i + 2]
        y = y + A * np.exp(-0.5 * ((u - m) / s) ** 2)
    return y


def fit_logspace(u, y, seeds_u):
    """Fit independent Gaussians in log-space; return params, pcov, modes."""
    p0, lo, hi = [], [], []
    span = u.max() - u.min()
    for su in seeds_u:
        p0 += [y.max(), su, 0.35]
        lo += [0.0, u.min() - 0.5 * span, 1e-2]
        hi += [np.inf, u.max() + 0.5 * span, span]
    popt, pcov = curve_fit(gauss_sum, u, y, p0=p0, bounds=(lo, hi), maxfev=20000)
    perr = np.sqrt(np.clip(np.diag(pcov), 0, None))
    modes = []
    for i in range(0, len(popt), 3):
        A, m, s = popt[i:i + 3]
        dA, dm, ds = perr[i:i + 3]
        # analytic mass (area under Gaussian in u) -> relative weight later
        mass = A * s * np.sqrt(2 * np.pi)
        modes.append(dict(A=A, mu=m, sigma=s, dmu=dm, dsigma=ds, mass=mass))
    total = sum(m["mass"] for m in modes) or 1.0
    for m in modes:
        m["weight"] = m["mass"] / total
    return popt, modes


def mixture_cdf(d, modes):
    u = np.log(d)
    c = 0.0
    for m in modes:
        c += m["weight"] * 0.5 * (1 + erf((u - m["mu"]) / (m["sigma"] * np.sqrt(2))))
    return c


def dvals_analytic(modes, ps=(10, 50, 90)):
    """Overall D-values from the analytic mixture CDF (root-find on log grid)."""
    dgrid = np.logspace(-1, 3, 20000)
    cdf = mixture_cdf(dgrid, modes)
    return [float(np.interp(p / 100, cdf, dgrid)) for p in ps]


def r2(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def run_case(ax_row, label, pf):
    swd = PARAMS[pf]
    sp = re.search(r"Size distribution (.*?) \[%\]", swd).group(1)
    path = path_for(label)
    df = load_data(path, preprocess=True, param_fit=swd, min_diameter=config.MIN_DIAMETER)

    # --- raw instrument bins (no interpolation) ---
    d_raw = df["Particle diameter  [µm]"].astype(float).values
    y_raw = df[swd].astype(float).values
    order = np.argsort(d_raw)
    d_raw, y_raw = d_raw[order], y_raw[order]
    u_raw = np.log(d_raw)
    yb = y_raw / np.trapezoid(y_raw, u_raw)  # density per unit ln(d)

    # --- current pipeline: cubic-interp grid + peak detect ---
    xg, yg = prepare_distribution_data_interpolation(df, swd, config.INTERP_POINTS)
    ap = config.PEAK_REL_PROMINENCE * float(np.max(yg))
    n_interp, pk_interp = find_distribution_peaks(xg, yg, prominence=ap,
                                                  distance=config.PEAK_DISTANCE)
    # peaks on the RAW bins (log-space) for comparison
    yb_n = yb / yb.max()
    pk_idx, _ = find_peaks(yb_n, prominence=0.02)
    pk_raw = d_raw[pk_idx]

    # current fit (clamped to 3 as the report does)
    n_cur = min(max(n_interp, 1), 3)
    fp, fdf, r2_cur = fit_multimodal_lognorm(xg, yg, n_cur, pk_interp, swd,
                                             weighting_scheme=None)
    cur_centers = sorted(fp[f"Mean{i+1}"] for i in range(n_cur))
    spurious = [c for c in cur_centers if c > 100]

    # prototype fit in log-space on raw bins
    seeds_u = np.log(pk_raw) if len(pk_raw) else np.array([u_raw[np.argmax(yb)]])
    popt, modes = fit_logspace(u_raw, yb, seeds_u)
    yhat = gauss_sum(u_raw, *popt)
    r2_proto = r2(yb, yhat)
    d10, d50, d90 = dvals_analytic(modes)
    # propagate D50 uncertainty from the dominant mode: dD/D ~ dmu
    dom = max(modes, key=lambda m: m["weight"])
    d50_rel_err = dom["dmu"]  # since dD/D ≈ dmu for the median of dominant mode

    # ---- console summary ----
    print(f"\n=== {label}  ({sp}) ===")
    print(f"  interp peaks detected : {n_interp}  centers(fit, clamped {n_cur}): "
          f"{[round(c,1) for c in cur_centers]}  spurious>100µm: "
          f"{[round(c,1) for c in spurious]}")
    print(f"  raw-bin peaks (logspc): {len(pk_raw)}  at {[round(c,2) for c in pk_raw]}")
    print(f"  R^2   current(lin-x): {r2_cur:.4f}    prototype(log-u): {r2_proto:.4f}")
    print(f"  prototype modes (geo-mean ± dμ, σg, wt):")
    for m in sorted(modes, key=lambda z: z["mu"]):
        print(f"      Dg={np.exp(m['mu']):7.2f} µm  μ={m['mu']:.3f}±{m['dmu']:.3f}"
              f"  σ={m['sigma']:.3f}  wt={m['weight']*100:5.1f}%")
    print(f"  prototype overall D10/D50/D90: {d10:.2f}/{d50:.2f}/{d90:.2f} µm"
          f"   (D50 rel. unc ≈ {d50_rel_err*100:.1f}%)")

    # ---- plots: (a) interp vs raw bins  (b) prototype fit ----
    axa, axb = ax_row
    axa.plot(d_raw, yb / yb.max(), "ko-", ms=3, lw=1, label="raw bins")
    axa.plot(xg, yg / yg.max(), "C1-", lw=1, label="cubic interp")
    for c in spurious:
        axa.axvline(c, color="r", ls=":", lw=1)
    axa.set_xscale("log"); axa.set_xlim(0.3, 1000)
    axa.set_title(f"{label} {sp}: interp vs raw"
                  + (f"  (spurious @{[round(c) for c in spurious]}µm)" if spurious else ""))
    axa.legend(fontsize=7); axa.grid(True, which="both", ls=":", alpha=.4)

    axb.plot(u_raw, yb, "ko", ms=3, label="raw (log-u)")
    uu = np.linspace(u_raw.min(), u_raw.max(), 600)
    axb.plot(uu, gauss_sum(uu, *popt), "C3-", lw=1.5, label=f"proto R²={r2_proto:.3f}")
    for m in modes:
        axb.plot(uu, m["A"] * np.exp(-0.5 * ((uu - m["mu"]) / m["sigma"]) ** 2),
                 "C0--", lw=.8)
    axb.set_title(f"{label}: Gaussian-in-ln(d) fit")
    axb.set_xlabel("u = ln(d)"); axb.legend(fontsize=7)
    axb.grid(True, ls=":", alpha=.4)

    return dict(label=label, sp=sp, r2_cur=r2_cur, r2_proto=r2_proto,
                spurious=spurious, d50_rel_err=d50_rel_err)


def main():
    fig, axes = plt.subplots(len(CASES), 2, figsize=(13, 3.0 * len(CASES)))
    axes = np.atleast_2d(axes)
    rows = [run_case(axes[i], lbl, pf) for i, (lbl, pf) in enumerate(CASES)]
    fig.suptitle("PSD fitting diagnosis: cubic-interp artefacts & log-space fit",
                 fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(OUT / "diagnosis.png", dpi=130)
    print(f"\nFigure -> {(OUT / 'diagnosis.png').resolve()}")


if __name__ == "__main__":
    main()
