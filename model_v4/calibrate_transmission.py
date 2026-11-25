"""Recalibrate the 7 transmission parameters on Starsim 4.

Row 868 (Adam Akullian's exp 024 point fit) was fitted on Starsim 3.7.2. On
Starsim 4 the same parameters give a larger epidemic (stisim's v4 HIV no
longer truncates the acute stage), and women's prevalence -- already too high
on 3.7.2 -- is further above PHIA. This refits the same 7 parameters, with
the exp 027d cascade, PrEP and VMMC inputs held fixed (baseline.cascade_kw()).

Parameters and bounds: exp 024's (hivsim_eswatini experiments/024_hm_wave1/
config.yaml), except beta_m2f, whose floor is lowered from 0.0096 to 0.004
because v4 transmits more at the same per-act probability.

Objective (a sum of squared z-scores; lower is better):
  1. PHIA prevalence, 54 strata (calibration_data/prevalence_by_age_sex.csv:
     year 2007/2011/2016 x sex x 5-year band), sigma = sqrt(se^2 + 0.02^2):
     the survey's standard error plus exp 024's model-discrepancy allowance.
  2. Prevalence 15-49 by sex and year (6), survey-weighted as
     standard_figures.targets_15_49, sigma = sqrt(se^2 + 0.01^2) (exp 025's
     halved allowance).
  3. Women 20-39 by year (3), sigma = sqrt(se^2 + 0.01^2): the researcher's
     priority (women's prevalence too high, worst at 20-40), stated
     explicitly so it can't be traded away against other strata.
  4. SHIMS incidence aggregates (calibration_data/incidence_by_age_sex.csv,
     the 4 rows exp 024 fitted: 2011 18-49 and 2016 15-49 by sex), sigma as
     in that file. The 18-49 rows pro-rate the 15-19 band (x0.4) as
     standard_figures._inc_series does.
Model prevalence in a year is the annual mean (PopByAgeSex stocks are
averaged over each year), as in the figures. 2021 (SHIMS3) is held out.

Search: optuna TPE (multivariate), BATCH trials evaluated in parallel on
SEEDS_SEARCH seeds each, row 868 enqueued as the first trial; the best
TOP_FINAL trials are then re-run on SEEDS_FINAL fresh seeds and the best of
those chosen. Runs stop at 2017 (the last target is 2016).

Run with the hivsim-v4 env:  python model_v4/calibrate_transmission.py
Outputs: ../raw_results/calib_v4/ (all trials, final re-runs) and
params/v4_transmission.json (the chosen parameters).
"""

import os
import sys
import time
from pathlib import Path

import numpy as np
import optuna
import pandas as pd
import sciris as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from baseline import make_baseline_sim, ROW868     # noqa: E402  (chdirs to HERE)
from analyzers import PopByAgeSex                   # noqa: E402
import standard_figures as sf                       # noqa: E402

RAW = HERE.parent / "raw_results" / "calib_v4"
RAW.mkdir(parents=True, exist_ok=True)

BOUNDS = dict(beta_m2f=(0.004, 0.025), rel_beta_f2m=(0.15, 0.60), s_f_young=(0.8, 3.0),
              age_gap_shift=(-2.0, 3.0), age_gap_sd_mult=(0.6, 1.8),
              prop_f0=(0.45, 0.85), prop_m0=(0.40, 0.80))
LOG = {"beta_m2f", "rel_beta_f2m", "s_f_young", "age_gap_sd_mult"}  # sampled on a log scale, as 024
NAMES = list(BOUNDS)
DISC_STRATUM, DISC_AGG = 0.02, 0.01
N_TRIALS, BATCH, SEEDS_SEARCH = 405, 9, 2
TOP_FINAL, SEEDS_FINAL = 10, 8
NCPUS = 9
STOP = 2017

# --- Targets -------------------------------------------------------------------
TG = sf.load_targets()                                    # year, sex, age_low, Count, phia, lb, ub
TG["se"] = (TG.ub - TG.lb) / (2 * 1.96)
TG["sigma"] = np.sqrt(TG.se**2 + DISC_STRATUM**2)


def _aggregate(tg, lo, hi):
    """Survey-weighted prevalence over 5-year bands [lo, hi), by year and sex, with a binomial se."""
    a = tg[tg.age_low.between(lo, hi - 5)]
    g = a.assign(num=a.phia * a.Count).groupby(["year", "sex"])
    out = pd.DataFrame({"phia": g.num.sum() / g.Count.sum(), "n": g.Count.sum()}).reset_index()
    out["sigma"] = np.sqrt(out.phia * (1 - out.phia) / out.n + DISC_AGG**2)
    return out


AGG_15_49 = _aggregate(TG, 15, 50)
WOMEN_20_39 = _aggregate(TG[TG.sex == "f"], 20, 40)
INC = sf.load_inc_targets()
INC = INC[INC.get("fit", True) == True] if "fit" in INC else INC  # noqa: E712
INC = INC[(INC.age_high - INC.age_low) > 20]               # the 4 sex-level aggregates


def model_tables(df):
    """Per-seed model prevalence by stratum, aggregates, and incidence, from one run's annual df."""
    df = df.rename(columns={"time": "timevec"})
    rows = {}
    for _, t in TG.iterrows():
        g = df[df.timevec == t.year]
        inf = g[f"popagesex.n_infected_{t.sex}_{t.age_low}_{t.age_low + 5}"].iloc[0]
        alv = g[f"popagesex.n_alive_{t.sex}_{t.age_low}_{t.age_low + 5}"].iloc[0]
        rows[("prev", t.year, t.sex, t.age_low)] = inf / alv if alv else np.nan

    def band_prev(year, sex, lo, hi):
        g = df[df.timevec == year]
        inf = sum(g[f"popagesex.n_infected_{sex}_{a}_{a + 5}"].iloc[0] for a in range(lo, hi, 5))
        alv = sum(g[f"popagesex.n_alive_{sex}_{a}_{a + 5}"].iloc[0] for a in range(lo, hi, 5))
        return inf / alv

    for _, t in AGG_15_49.iterrows():
        rows[("agg", t.year, t.sex, 15)] = band_prev(t.year, t.sex, 15, 50)
    for _, t in WOMEN_20_39.iterrows():
        rows[("w2039", t.year, "f", 20)] = band_prev(t.year, "f", 20, 40)
    for _, t in INC.iterrows():
        lo = int(t.age_low)
        first = (lo // 5) * 5
        w = (5 - (lo - first)) / 5                         # pro-rate the first band (18-49 -> 0.4 x 15-19)
        g = df[df.timevec == t.year]
        new = alv = inf = 0.0
        for i, a in enumerate(range(first, int(t.age_high), 5)):
            k = w if i == 0 else 1.0
            new += k * g[f"popagesex.new_infections_{t.sex}_{a}_{a + 5}"].iloc[0]
            alv += k * g[f"popagesex.n_alive_{t.sex}_{a}_{a + 5}"].iloc[0]
            inf += k * g[f"popagesex.n_infected_{t.sex}_{a}_{a + 5}"].iloc[0]
        rows[("inc", t.year, t.sex, lo)] = 100.0 * new / (alv - inf)
    return rows


def run_point(pid, x, seed):
    sim = make_baseline_sim(seed=seed, stop=STOP, pars=x, analyzers=[PopByAgeSex()])
    sim.run()
    out = model_tables(sim.to_df(resample="year", use_years=True, sep="."))
    return pid, seed, out


def score(results):
    """results: list of (pid, seed, tables). Mean over seeds, then z-scores."""
    by_pid = {}
    for pid, seed, tab in results:
        by_pid.setdefault(pid, []).append(tab)
    rows = []
    for pid, tabs in by_pid.items():
        m = {k: float(np.nanmean([t[k] for t in tabs])) for k in tabs[0]}
        r = dict(pid=pid)
        z_prev = [((m[("prev", t.year, t.sex, t.age_low)] - t.phia) / t.sigma) for _, t in TG.iterrows()]
        z_agg = [((m[("agg", t.year, t.sex, 15)] - t.phia) / t.sigma) for _, t in AGG_15_49.iterrows()]
        z_w = [((m[("w2039", t.year, "f", 20)] - t.phia) / t.sigma) for _, t in WOMEN_20_39.iterrows()]
        z_inc = [((m[("inc", t.year, t.sex, int(t.age_low))] - t.incidence_pct) / t.sigma) for _, t in INC.iterrows()]
        r["loss_prev"], r["loss_agg"] = float(np.sum(np.square(z_prev))), float(np.sum(np.square(z_agg)))
        r["loss_w2039"], r["loss_inc"] = float(np.sum(np.square(z_w))), float(np.sum(np.square(z_inc)))
        r["loss"] = r["loss_prev"] + r["loss_agg"] + r["loss_w2039"] + r["loss_inc"]
        r["mae_prev"] = float(np.mean([abs(m[("prev", t.year, t.sex, t.age_low)] - t.phia) for _, t in TG.iterrows()]))
        r["bias_prev"] = float(np.mean([m[("prev", t.year, t.sex, t.age_low)] - t.phia for _, t in TG.iterrows()]))
        r["n_within_ci"] = int(sum(t.lb <= m[("prev", t.year, t.sex, t.age_low)] <= t.ub for _, t in TG.iterrows()))
        for _, t in AGG_15_49.iterrows():
            r[f"prev15_49_{t.sex}_{t.year}"] = m[("agg", t.year, t.sex, 15)]
        for _, t in WOMEN_20_39.iterrows():
            r[f"prevw20_39_{t.year}"] = m[("w2039", t.year, "f", 20)]
        for _, t in INC.iterrows():
            r[f"inc_{t.sex}_{t.year}"] = m[("inc", t.year, t.sex, int(t.age_low))]
        rows.append(r)
    return pd.DataFrame(rows)


def evaluate(points, seeds):
    args = [(pid, x, s) for pid, x in points for s in seeds]
    return score(sc.parallelize(run_point, iterarg=args, ncpus=NCPUS))


if __name__ == "__main__":
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    t0 = time.time()
    study = optuna.create_study(direction="minimize", sampler=optuna.samplers.TPESampler(
        seed=20261119, multivariate=True, constant_liar=True, n_startup_trials=45))
    study.enqueue_trial(dict(ROW868))
    records = []
    while len(study.trials) < N_TRIALS:
        trials = [study.ask() for _ in range(BATCH)]
        xs = [{k: t.suggest_float(k, *BOUNDS[k], log=k in LOG) for k in NAMES} for t in trials]
        s = evaluate([(t.number, x) for t, x in zip(trials, xs)], range(1, SEEDS_SEARCH + 1))
        for t, x in zip(trials, xs):
            r = s[s.pid == t.number].iloc[0]
            study.tell(t, float(r.loss))
            records.append({**x, **r.to_dict()})
        best = min(records, key=lambda d: d["loss"])
        print(f"trials {len(study.trials):4d}  {time.time() - t0:6.0f}s  best loss {best['loss']:.1f} "
              f"(prev {best['loss_prev']:.1f}, agg {best['loss_agg']:.1f}, w20-39 {best['loss_w2039']:.1f}, "
              f"inc {best['loss_inc']:.1f}; MAE {best['mae_prev']:.3f})", flush=True)
    search = pd.DataFrame(records)
    search.to_csv(RAW / "search.csv", index=False)

    top = search.nsmallest(TOP_FINAL, "loss")
    pts = [(20_000 + i, {k: float(r[k]) for k in NAMES}) for i, (_, r) in enumerate(top.iterrows())]
    pts.append((29_999, dict(ROW868)))                       # row 868 on the same fresh seeds, for reference
    fin = evaluate(pts, range(101, 101 + SEEDS_FINAL))
    fin = pd.DataFrame([{**x, "pid": pid} for pid, x in pts]).merge(fin, on="pid")
    fin.to_csv(RAW / "final.csv", index=False)
    chosen = fin[fin.pid != 29_999].nsmallest(1, "loss").iloc[0]
    sc.savejson(HERE / "params" / "v4_transmission.json",
                {**{k: float(chosen[k]) for k in NAMES},
                 "loss": float(chosen.loss), "mae_prev": float(chosen.mae_prev),
                 "bias_prev": float(chosen.bias_prev), "n_within_ci": int(chosen.n_within_ci),
                 "fitted_on": "starsim 4.0.0rc5, stisim v4 @ 9fe22a2; 027d cascade inputs fixed",
                 "seeds_final": SEEDS_FINAL})

    cols = NAMES + ["loss", "loss_prev", "loss_agg", "loss_w2039", "loss_inc", "mae_prev", "bias_prev", "n_within_ci"]
    with pd.option_context("display.width", 260, "display.float_format", "{:.4g}".format):
        print("\nFinal (8 fresh seeds); pid 29999 = row 868:")
        print(fin.sort_values("loss")[["pid"] + cols].to_string(index=False))
        ref = fin[fin.pid == 29_999].iloc[0]
        print("\n15-49 prevalence and women 20-39, chosen vs row 868 vs PHIA:")
        for _, t in AGG_15_49.iterrows():
            k = f"prev15_49_{t.sex}_{t.year}"
            print(f"  {t.sex} 15-49 {t.year}: {chosen[k]:.3f}  (row 868 {ref[k]:.3f}; PHIA {t.phia:.3f})")
        for _, t in WOMEN_20_39.iterrows():
            k = f"prevw20_39_{t.year}"
            print(f"  f 20-39 {t.year}: {chosen[k]:.3f}  (row 868 {ref[k]:.3f}; PHIA {t.phia:.3f})")
        for _, t in INC.iterrows():
            k = f"inc_{t.sex}_{t.year}"
            print(f"  incidence {t.sex} {t.year}: {chosen[k]:.2f}  (row 868 {ref[k]:.2f}; SHIMS {t.incidence_pct:.2f})")
    print(f"done in {time.time() - t0:.0f}s")
