"""
PSD fitting core — sum of Gaussians in log-diameter.

A particle-size mode is lognormal in diameter d, i.e. a *Gaussian* in the
natural coordinate u = ln(d). Fitting in u (rather than fitting lognorm.pdf in
linear d on a log grid) keeps the model well-conditioned, makes the residual
weighting uniform across the decades you actually plot, and makes every mode
statistic analytic:

    median   D50_i = exp(mu_i)
    quantile D_p_i = exp(mu_i + sigma_i * z_p),   z_p = Phi^-1(p/100)
    mean         = exp(mu_i + sigma_i**2 / 2)        (arithmetic, linear d)

Because d is exponential in the parameters, dD/D = dmu + z_p*dsigma: a small
*absolute* error in log-space is a constant *relative* error in diameter. We
therefore capture the fit covariance and report every diameter with a
propagated uncertainty instead of hiding it behind extra digits.

The model fits N independent, non-negative Gaussian amplitudes (no 1-w1
simplex, no max(0, .) kinks) so a single code path covers any number of modes.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import ndtri, erf  # inverse normal CDF, error function

SQRT_2PI = np.sqrt(2.0 * np.pi)


def gaussian_sum(u, *params):
    """Sum of Gaussians in u = ln(d).  params = [A1, mu1, s1, A2, mu2, s2, ...]."""
    u = np.asarray(u, dtype=float)
    y = np.zeros_like(u)
    for i in range(0, len(params), 3):
        A, mu, s = params[i], params[i + 1], params[i + 2]
        y = y + A * np.exp(-0.5 * ((u - mu) / s) ** 2)
    return y


def _point_weights(y, weighting):
    """Optional curve_fit sigma array. None = unweighted (recommended)."""
    if weighting is None:
        return None
    y = np.asarray(y, dtype=float)
    ymax = np.max(y) if np.any(y > 0) else 1.0
    floor = 1e-3 * ymax
    if weighting == "sqrt_y":          # emphasise low-amplitude modes
        w = np.sqrt(np.maximum(y, floor))
    elif weighting == "proportional_to_y":
        w = np.maximum(y, floor)
    else:
        raise ValueError(f"unknown weighting: {weighting!r}")
    return 1.0 / np.sqrt(w)


def _mode_statistics(A, mu, sigma, sub_cov, percentiles=(10, 50, 90)):
    """Analytic lognormal statistics for one Gaussian-in-u mode.

    sub_cov is the 2x2 covariance of (mu, sigma) for this mode, used to
    propagate the diameter uncertainties. Returns a plain dict.
    """
    stats = {
        "A": A, "mu": mu, "sigma": sigma,
        "Dg": float(np.exp(mu)),            # geometric mean = median
        "sigma_g": float(np.exp(sigma)),    # geometric standard deviation
        "Mean": float(np.exp(mu + sigma ** 2 / 2.0)),
        "Std": float(np.sqrt(max(np.exp(sigma ** 2) - 1.0, 0.0))
                     * np.exp(mu + sigma ** 2 / 2.0)),
        "mu_err": float(np.sqrt(max(sub_cov[0, 0], 0.0))),
        "sigma_err": float(np.sqrt(max(sub_cov[1, 1], 0.0))),
    }
    for p in percentiles:
        z = ndtri(p / 100.0)
        Dp = np.exp(mu + sigma * z)
        # var(ln Dp) = var(mu) + z^2 var(sigma) + 2 z cov(mu, sigma)
        var_lnDp = (sub_cov[0, 0] + z ** 2 * sub_cov[1, 1]
                    + 2.0 * z * sub_cov[0, 1])
        stats[f"D{p}"] = float(Dp)
        stats[f"D{p}_rel_err"] = float(np.sqrt(max(var_lnDp, 0.0)))  # = dDp/Dp
    return stats


def _mixture_cdf(d, modes):
    u = np.log(np.asarray(d, dtype=float))
    c = np.zeros_like(u)
    for m in modes:
        z = (u - m["mu"]) / (m["sigma"] * np.sqrt(2.0))
        c = c + m["weight"] * 0.5 * (1.0 + erf(z))
    return c


def _mixture_percentiles(modes, percentiles=(10, 50, 90)):
    """Overall D-values of the fitted mixture, by inverting its analytic CDF."""
    if not modes:
        return {p: float("nan") for p in percentiles}
    grid = np.logspace(-2, 3.5, 30000)
    cdf = _mixture_cdf(grid, modes)
    return {p: float(np.interp(p / 100.0, cdf, grid)) for p in percentiles}


def fit_psd(d, y, num_modes, seeds, weighting=None, default_sigma=0.35,
            mc_samples=300, rng_seed=0):
    """
    Fit `num_modes` Gaussians in u = ln(d) to the distribution (d, y).

    Parameters
    ----------
    d, y : array        raw instrument bins (diameter, differential weight).
    num_modes : int     number of modes to fit (operator-controlled).
    seeds : array       seed diameters (e.g. detected peaks). Missing seeds are
                        added in the small-particle region, where the elusive
                        fine mode usually hides.
    weighting : str|None  None (unweighted) | 'sqrt_y' | 'proportional_to_y'.
    mc_samples : int    Monte-Carlo draws for overall D-value CIs (0 disables).

    Returns a dict:
        modes      : list of per-mode stat dicts (sorted by diameter), each with
                     weight, Dg, sigma_g, D10/D50/D90 (+ *_rel_err), Mean, Std.
        overall    : fitted-mixture D10/D50/D90 (+ _ci when mc_samples>0).
        r_squared  : goodness of fit in log-space.
        n_modes    : number of valid (in-range) modes returned.
        dropped    : seed diameters of modes rejected as out-of-range (runaways).
        u, y_norm, y_fit : arrays for plotting (log-diameter space).
    """
    d = np.asarray(d, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(d)
    d, y = d[order], y[order]
    u = np.log(d)
    u_lo, u_hi = float(u.min()), float(u.max())
    span = u_hi - u_lo

    # Normalise to unit area in u so amplitudes are comparable between files.
    area = np.trapezoid(y, u)
    y_norm = y / area if area > 0 else y

    # Build seed list in log-space; pad missing seeds into the fine-particle end.
    seeds_u = sorted(float(np.log(s)) for s in np.atleast_1d(seeds) if s > 0)
    seeds_u = seeds_u[:num_modes]
    if len(seeds_u) < num_modes:
        upper = min(seeds_u) if seeds_u else u[int(np.argmax(y_norm))]
        extra = np.linspace(u_lo, upper, (num_modes - len(seeds_u)) + 2)[1:-1]
        seeds_u = sorted(seeds_u + list(extra))

    ymax = float(np.max(y_norm))
    p0, lo, hi = [], [], []
    for su in seeds_u:
        p0 += [ymax, su, default_sigma]
        lo += [0.0, u_lo - 0.25 * span, 1e-2]
        hi += [np.inf, u_hi + 0.25 * span, max(span, 1e-2)]

    sigma_pts = _point_weights(y_norm, weighting)
    try:
        popt, pcov = curve_fit(gaussian_sum, u, y_norm, p0=p0, bounds=(lo, hi),
                               sigma=sigma_pts, absolute_sigma=False, maxfev=20000)
    except RuntimeError:
        popt = np.asarray(p0, dtype=float)
        pcov = np.full((len(p0), len(p0)), np.nan)

    y_fit = gaussian_sum(u, *popt)
    ss_res = np.sum((y_norm - y_fit) ** 2)
    ss_tot = np.sum((y_norm - np.mean(y_norm)) ** 2)
    r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0

    # Assemble modes; reject any whose centre ran outside the measured range.
    raw_masses, kept, dropped = [], [], []
    for i in range(0, len(popt), 3):
        A, mu, s = popt[i], popt[i + 1], popt[i + 2]
        if not (u_lo <= mu <= u_hi):
            dropped.append(float(np.exp(mu)))
            continue
        sub_cov = pcov[i + 1:i + 3, i + 1:i + 3] if np.all(np.isfinite(pcov)) \
            else np.zeros((2, 2))
        kept.append((A, mu, s, sub_cov))
        raw_masses.append(A * s * SQRT_2PI)        # area under Gaussian in u

    total_mass = sum(raw_masses) or 1.0
    modes = []
    for (A, mu, s, sub_cov), mass in zip(kept, raw_masses):
        st = _mode_statistics(A, mu, s, sub_cov)
        st["weight"] = float(mass / total_mass)
        modes.append(st)
    modes.sort(key=lambda m: m["Dg"])

    overall = _mixture_percentiles(modes)
    if mc_samples and np.all(np.isfinite(pcov)) and modes:
        overall_ci = _mc_overall_ci(popt, pcov, u_lo, u_hi, mc_samples, rng_seed)
        for p in (10, 50, 90):
            overall[f"D{p}_ci"] = overall_ci[p]

    return {
        "modes": modes, "overall": overall, "r_squared": r_squared,
        "n_modes": len(modes), "dropped": dropped,
        "u": u, "y_norm": y_norm, "y_fit": y_fit,
    }


def _mc_overall_ci(popt, pcov, u_lo, u_hi, mc_samples, rng_seed):
    """68% CI on overall D-values by sampling the fit covariance."""
    rng = np.random.default_rng(rng_seed)
    draws = rng.multivariate_normal(popt, pcov, size=mc_samples)
    acc = {10: [], 50: [], 90: []}
    for params in draws:
        modes = []
        masses = []
        for i in range(0, len(params), 3):
            A, mu, s = params[i], params[i + 1], params[i + 2]
            if A <= 0 or s <= 0 or not (u_lo <= mu <= u_hi):
                continue
            modes.append({"mu": mu, "sigma": s})
            masses.append(A * s * SQRT_2PI)
        if not modes:
            continue
        tot = sum(masses) or 1.0
        for m, mass in zip(modes, masses):
            m["weight"] = mass / tot
        q = _mixture_percentiles(modes)
        for p in (10, 50, 90):
            acc[p].append(q[p])
    out = {}
    for p in (10, 50, 90):
        vals = np.asarray(acc[p], dtype=float)
        out[p] = (float(np.percentile(vals, 16)), float(np.percentile(vals, 84))) \
            if vals.size else (float("nan"), float("nan"))
    return out
