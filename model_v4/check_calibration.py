"""Check the Starsim v4 Eswatini model against the calibration targets.

Port of hivsim_eswatini's DTC_check/build_and_check_calibration.py (branch
exp/027-cascade-prep-recal, commit 068dd83). Everything from the "Save
results" cell on is that script's plotting code, unchanged apart from paths;
only this header differs: it builds the model with model_v4's
`baseline.make_baseline_sim()` (Starsim 4) instead of hivsim_eswatini's
run_sims and config switch.

Run with the hivsim-v4 conda env, from anywhere except the IDM root:
    python model_v4/check_calibration.py
Cell-marked (`# %%`) for VS Code's Interactive window.

Starsim 4 differences that matter when reading the figures
------------------------------------------------------------
- `sim.to_df()` names the time column `time`; it is renamed `timevec` below so
  the plotting code (and standard_figures.py) is unchanged.
- The age/sex and cascade counts (analyzers PopByAgeSex, CascadeByAgeSex) are
  recorded at the start of each step and averaged over each year, so a year's
  value is its mean, as in the Starsim 3 figures (v4's default for a stock is
  its value on 1 January).
- The model's own results follow v4 names (e.g. `hiv.on_art`, not `n_on_art`).
- The transmission parameters are the v4 refit (baseline.BEST, from
  calibrate_transmission.py); the testing rates are exp 027d's, fitted on
  Starsim 3.7.2 with row 868 and not refitted.

Data: calibration_data/ (prevalence, incidence and ART-coverage targets, 2021
hold-out) and data/ (SHIMS 90-90-90, SHIMS incidence, Kasturi, eswatini_ppdv,
SHIMS2 incidence by age, PrEP initiations), copied from hivsim_eswatini.
Outputs: figures and the comparison table in ../results/model_v4_check/; the
run itself in ../raw_results/model_v4_check/ (gitignored). With --fsw, the
sex-work fit instead, in ../results/model_v4_fsw_check/.
"""

# %% Imports
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent  # model_v4/ (laid out like hivsim_eswatini)
sys.path.insert(0, str(REPO))

from baseline import make_baseline_sim, cascade_kw, ROW868, SW_BASE  # noqa: E402  (chdirs to REPO)
from analyzers import PopByAgeSex, CascadeByAgeSex        # noqa: E402
import standard_figures as sf                             # noqa: E402

# Which fit to check: "v4" is baseline.BEST; "v4_fsw" is the 10-parameter fit
# with the sex-work parameters free (calibrate_transmission.py --fsw). Run with
# --fsw, or set it here in the Interactive window.
FIT = "v4_fsw" if "--fsw" in sys.argv else "v4"
if FIT == "v4_fsw":
    PARS = {k: v for k, v in sc.loadjson(REPO / "params" / "v4_transmission_fsw.json").items()
            if k in ROW868 or k in SW_BASE}
    CHECK = "model_v4_fsw_check"
    LABEL = "model_v4: v4 refit with sex-work parameters free + 027d inputs on Starsim 4"
else:
    PARS = None  # baseline.BEST
    CHECK = "model_v4_check"
    LABEL = "model_v4: v4-refitted transmission + 027d inputs on Starsim 4"
OUT = REPO.parent / "results" / CHECK
RAW = REPO.parent / "raw_results" / CHECK
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)
MODEL_KW = cascade_kw()  # the 027d settings make_baseline_sim uses (read by the plots, e.g. shims2_year)

N_SEEDS = 10
START, STOP = 1985, 2026


# %% Build + run the model on Starsim 4, one seed per process
def run_one_seed(seed):
    sim = make_baseline_sim(seed=seed, start=START, stop=STOP, pars=PARS,
                            analyzers=[PopByAgeSex(), CascadeByAgeSex()])
    sim.run()
    d = sim.to_df(resample="year", use_years=True, sep=".")
    d = d.rename(columns={"time": "timevec"})  # Starsim 3's name, used by the plots
    d["seed"] = seed
    return d


print(f"Running {N_SEEDS} seeds of the v4 model ({START}-{STOP})...")
dfs = sc.parallelize(run_one_seed, iterarg=list(range(1, N_SEEDS + 1)),
                     ncpus=min(N_SEEDS, 8))
df = pd.concat(dfs, ignore_index=True)
print(f"Done: {len(df)} rows across {N_SEEDS} seeds")

# %% Save results
# sc.saveobj (pickle via sciris), matching run_sims.py's own __main__ convention
# (results/eswatini_sim.df) -- avoids a pyarrow dependency this env doesn't have.
sim_path = RAW / "model_v4_sim.df"
sc.saveobj(sim_path, df)
print(f"Saved sim output to {sim_path}")

# %% Extract and plot age-prevalence curves against PHIA
# standard_figures.py is the repo's single shared implementation of the
# age-stratified prevalence-fit figure -- every experiment's run.py calls this
# same function, so this figure is directly comparable to 024's and 026's.
# Its docstring has the sex-coding note for calibration_data/
# prevalence_by_age_sex.csv: Gender==0 is MALE there (the OPPOSITE of
# stisim's internal 0=female convention) -- plot_prevalence_fit already
# handles this; just don't reach into that CSV's Gender column directly
# without re-reading that note.
fig_path = OUT / "prevalence_fit_v4_check.png"
scorecard = sf.plot_prevalence_fit(df, label=f"v4 check ({LABEL})",
                                   outpath=fig_path)
print(f"Saved age-prevalence fit figure to {fig_path}")
print("Scorecard (lower mae/bias = better fit):")
for k, v in scorecard.items():
    print(f"  {k}: {v}")

# %% Prevalence including 2021 (validation hold-out), alongside the 3 fitted years
# 2021 was NOT used to fit 024's row 868 -- it's CLAUDE.md's validation
# hold-out. Included here as an out-of-sample check only, by extending the
# same shared figure via its `tg` override, so the computation (denominators,
# CI shading, scorecard) is identical to the calibration-only figure above --
# the only thing that changes is which rows are in `tg`.
v21 = pd.read_csv(REPO / "calibration_data" / "prevalence_2021_VALIDATION_ONLY.csv")
v21["sex"] = v21.Gender.map(sf.PHIA_SEX)         # same 0=M,1=F convention as the fitted file
v21["age_low"] = v21["start age"].astype(int)
v21 = v21[["Year", "sex", "age_low", "Count", "NationalPrevalence", "lb", "ub"]].rename(
    columns={"Year": "year", "NationalPrevalence": "phia"})
tg_incl_2021 = pd.concat([sf.load_targets(), v21], ignore_index=True)

fig_path_2021 = OUT / "prevalence_fit_incl_2021_v4_check.png"
scorecard_2021 = sf.plot_prevalence_fit(
    df, label=f"v4 check ({LABEL}) incl. 2021 VALIDATION hold-out",
    outpath=fig_path_2021, tg=tg_incl_2021)
print(f"Saved prevalence-incl-2021 figure to {fig_path_2021}")
print("Scorecard incl. 2021 (those strata were NOT fitted -- read as out-of-sample):")
for k, v in scorecard_2021.items():
    print(f"  {k}: {v}")

# %% Incidence: time series + 2016 age bands, vs the calibration targets
# standard_figures.py's existing shared implementation -- same reuse logic as
# the prevalence figure above. 15-49 time series by sex against the fitted
# targets (SHIMS1 2011 18-49 cohort estimate, SHIMS2 2016 15-49), plus the
# 2016 age-band comparison, which carries CLAUDE.md's false-recency-bias
# caveat on the figure itself. The comparison against the SHIMS2/SHIMS3
# incidence file in SHIMS_Datasets/ is further down (incidence checks).
fig_inc_fit = OUT / "incidence_fit_v4_check.png"
sf.plot_incidence_fit(df, label=f"v4 check ({LABEL})", outpath=fig_inc_fit)
print(f"Saved incidence time-series + 2016 age-band figure to {fig_inc_fit}")

# %% 90-90-90 cascade over time (15-49), by sex
# Three stages, each CONDITIONAL on the previous stage (the standard UNAIDS
# framing):
#   1st: % of PLHIV who know their status    = n_diagnosed     / n_infected
#   2nd: % of diagnosed who are on ART        = n_on_art        / n_diagnosed
#   3rd: % of those on ART who are suppressed = n_effective_art / n_on_art
#
# Data availability -- checked against every file in calibration_data/ and
# external_data/ before writing this, not assumed:
#   1st (diagnosed): no target exists in calibration_data/ or external_data/.
#     Model-only here; compared against DTC_check/SHIMS_Datasets/
#     (SHIMS 2021 + Kasturi's estimates) in checks A-D further down.
#   2nd (on ART): calibration_data/art_coverage_by_age_sex.csv exists, but its
#     denominator is ALL PLHIV (it's what sti.ART's `coverage` input targets --
#     see interventions.py), not diagnosed. That's a DIFFERENT quantity from
#     "of diagnosed", so it's compared on the age-stratified figure below
#     (where that file's age bins live) against the model's matching-
#     denominator line, not against this "of diagnosed" line.
#   3rd (VLS): data/eswatini_ppdv.csv's VLS column is suppression among ALL
#     PLHIV, not among those on ART (see analyzers.Cascade's own docstring on
#     this exact distinction: conflating the two "double-counts the coverage
#     ramp"). So this panel plots BOTH the "of those on ART" line (what was
#     asked for) and a PLHIV-denominator line for a fair data overlay.
#   eswatini_ppdv.csv's ART_prev column is explicitly flagged by
#   ppdv_construction.py as "INDETERMINATE - do not use without checking" and
#   is deliberately not used anywhere in this script.
AGE_BINS = CascadeByAgeSex.AGE_BINS


def cascade_agg(df, lo=15, hi=50):
    """Per-seed, per-year cascade fractions summed over [lo,hi), long-form."""
    bins = [(a, b) for a, b in AGE_BINS if a >= lo and b <= hi]
    counts = {}
    for sex in ("f", "m"):
        for stem in ("n_infected", "n_diagnosed", "n_on_art", "n_effective_art"):
            cols = [f"cascadeagesex.{stem}_{sex}_{a}_{b}" for a, b in bins]
            counts[f"{stem}_{sex}"] = df[cols].sum(axis=1)
    for stem in ("n_infected", "n_diagnosed", "n_on_art", "n_effective_art"):
        counts[f"{stem}_all"] = counts[f"{stem}_f"] + counts[f"{stem}_m"]

    rows = []
    for sex in ("f", "m", "all"):
        def c(stem, sex=sex):
            return counts[f"{stem}_{sex}"]
        rows.append(pd.DataFrame({
            "year": df["timevec"], "seed": df["seed"], "sex": sex,
            "p_dx": c("n_diagnosed") / c("n_infected").replace(0, np.nan),
            "p_art_dx": c("n_on_art") / c("n_diagnosed").replace(0, np.nan),
            "p_art_plhiv": c("n_on_art") / c("n_infected").replace(0, np.nan),
            "p_vls_art": c("n_effective_art") / c("n_on_art").replace(0, np.nan),
            "p_vls_plhiv": c("n_effective_art") / c("n_infected").replace(0, np.nan),
        }))
    return pd.concat(rows, ignore_index=True)


casc = cascade_agg(df)
CASC_COLORS = {"f": "#c0392b", "m": "#2c6fbb", "all": "#444444"}
CASC_NAMES = {"f": "Women", "m": "Men", "all": "Both sexes"}


def _band(ax, g, col, color, label, ls="-", lw=2):
    s = g.groupby("year")[col]
    mean, sd = s.mean(), s.std(ddof=1).fillna(0.0)
    ax.plot(mean.index, mean.values * 100, color=color, lw=lw, ls=ls, label=label)
    if ls == "-":
        ax.fill_between(mean.index, ((mean - sd).clip(lower=0) * 100).values,
                        ((mean + sd).clip(upper=1) * 100).values,
                        color=color, alpha=0.15, lw=0)


fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]
for sex in ("f", "m", "all"):
    _band(ax, casc[casc.sex == sex], "p_dx", CASC_COLORS[sex], CASC_NAMES[sex])
ax.set_title("% of PLHIV diagnosed\n(data comparisons: see *_vs_shims / *_vs_kasturi figures)",
             fontsize=9.5)
ax.set_ylabel("%")

ax = axes[1]
for sex in ("f", "m", "all"):
    _band(ax, casc[casc.sex == sex], "p_art_dx", CASC_COLORS[sex], CASC_NAMES[sex])
ax.set_title("% of diagnosed on ART\n(data uses a different denominator -- "
            "see age-stratified figure)", fontsize=9.5)

ax = axes[2]
for sex in ("f", "m"):
    _band(ax, casc[casc.sex == sex], "p_vls_art", CASC_COLORS[sex],
         f"{CASC_NAMES[sex]} — of those on ART")
    _band(ax, casc[casc.sex == sex], "p_vls_plhiv", CASC_COLORS[sex],
         f"{CASC_NAMES[sex]} — of all PLHIV (matches data)", ls="--", lw=1.4)
ppdv = pd.read_csv(REPO / "data" / "eswatini_ppdv.csv")
for _, r in ppdv.iterrows():
    yerr = ([[(r.VLS - r.VLS_lb) * 100], [(r.VLS_ub - r.VLS) * 100]]
           if pd.notna(r.VLS_lb) else None)
    ax.errorbar(r.Year, r.VLS * 100, yerr=yerr, fmt="o", mfc=CASC_COLORS[r.sex],
                mec="black", mew=1.2, ms=8, ecolor="black", capsize=4, zorder=6)
    ax.annotate(r.AgeCat, (r.Year, r.VLS * 100), textcoords="offset points",
               xytext=(5, 5), fontsize=6.5, color="grey")
ax.errorbar([], [], fmt="o", mfc="grey", mec="black", ms=8, capsize=4,
           label="SHIMS/PHIA VLS among all PLHIV\n(label = survey's own age range)")
ax.set_title("% virally suppressed", fontsize=9.5)

for ax in axes:
    ax.set_xlabel("year"); ax.set_xlim(2000, STOP); ax.set_ylim(0, 105)
    ax.grid(alpha=0.3); ax.legend(fontsize=7)
fig.suptitle(f"90-90-90 cascade over time, ages 15-49 (v4 check, {LABEL})\n"
            "line = mean of 10 seeds, band = +/- 1 SD (diagnosed/on-ART only)",
            y=1.06, fontsize=11)
fig_casc_ts = OUT / "cascade_timeseries_v4_check.png"
fig.savefig(fig_casc_ts, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved 90-90-90 time-series figure to {fig_casc_ts}")

# %% Load the two external 90-90-90 datasets (DTC_check/SHIMS_Datasets/)
# Both are local to DTC_check, not calibration_data/ -- they are comparison
# data for this check, NOT fitting targets (CLAUDE.md excludes ART coverage).
#
# SHIMS_909090_coverage.csv -- all three stages, each CONDITIONAL on the
#   previous one (Diagnosed = of PLHIV; Diagnosed_OnART = of diagnosed;
#   OnART_VLS = of on ART), so it maps 1:1 onto p_dx / p_art_dx / p_vls_art.
#   `Number` is the stage's own denominator (it shrinks stage to stage).
#   Two survey rounds: SHIMS3 (2021) and SHIMS2 (labelled 2017 in the file;
#   fieldwork 2016-17 -- note the prevalence targets file this round as 2016).
#   CIs below are Wilson 95% on the unweighted `Number`, which ignores the
#   survey design effect, so they are narrower than the published CIs.
# aware_and_art_no_interactions_modeled... (Kasturi) -- modelled smooth
#   estimates, 1995-2021, by sex and 5-year band (15-19 .. 80-84), first two
#   stages only. art_extrapolated is ON ART AMONG AWARE, not among PLHIV:
#   it exceeds aware_extrapolated in ~half the rows, which is only possible
#   for the conditional quantity. No CIs are provided.
DATA_DIR = REPO / "data"  # SHIMS 90-90-90, SHIMS incidence, Kasturi (from DTC_check/SHIMS_Datasets)
STAGE_OF_METRIC = {"Diagnosed": "p_dx", "Diagnosed_OnART": "p_art_dx",
                   "OnART_VLS": "p_vls_art"}
STAGE_TITLES = {"p_dx": "% diagnosed (of PLHIV)",
                "p_art_dx": "% on ART (of diagnosed)",
                "p_vls_art": "% suppressed (of on ART)"}
SEX_OF_GENDER = {"female": "f", "male": "m"}  # matched lower-cased


def parse_age_group(s):
    """'15-24' -> (15, 25); '50+' -> (50, 100). Upper bound is exclusive."""
    if s.endswith("+"):
        return int(s[:-1]), 100
    lo, hi = s.split("-")
    return int(lo), int(hi) + 1


def wilson(p, n, z=1.96):
    """Wilson score interval for a proportion p observed on n trials."""
    centre = (p + z**2 / (2 * n)) / (1 + z**2 / n)
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / (1 + z**2 / n)
    return centre - half, centre + half


def load_shims(path):
    """Read the SHIMS 90-90-90 file and add model-facing columns."""
    s = pd.read_csv(path, encoding="utf-8-sig")
    assert s.Year.notna().all(), "SHIMS rows with no year"
    assert (s.Percentage.between(0, 100)).all(), "SHIMS percentage out of range"
    s["stage"] = s.Metric.map(STAGE_OF_METRIC)
    s["sex"] = s.Gender.str.lower().map(SEX_OF_GENDER)
    assert s[["stage", "sex"]].notna().all().all(), "unrecognised Metric/Gender"
    assert not s.duplicated(["Year", "stage", "sex", "Age group"]).any()
    s["p"] = s.Percentage / 100
    s[["age_lo", "age_hi"]] = s["Age group"].apply(
        lambda a: pd.Series(parse_age_group(a)))
    s["lb"], s["ub"] = wilson(s.p, s.Number)
    return s


shims = load_shims(DATA_DIR / "SHIMS_909090_coverage.csv")
SHIMS_YEARS = sorted(shims.Year.unique())
print("SHIMS rows per survey year:", shims.Year.value_counts().sort_index().to_dict())

kas = pd.read_csv(DATA_DIR / "aware_and_art_no_interactions_modeled_output_1995_2021-From_Kasturi.csv",
                  index_col=0)
kas["sex"] = kas.gender.str.lower().map(SEX_OF_GENDER)
kas["age_low"] = kas.agecat.str.split("-").str[0].astype(int)
kas = kas.rename(columns={"aware_extrapolated": "p_dx",
                          "art_extrapolated": "p_art_dx"})

# SHIMS's aggregate rows: (lo, hi) -> label. Used for both the time-series
# (15-49, 15+) and the by-age figure (the four non-overlapping bands).
AGG_GROUPS = {(15, 50): "15-49", (15, 100): "15+"}
SHIMS_AGE_BANDS = [(15, 25), (25, 35), (35, 50), (50, 100)]


def model_summary(df, lo, hi):
    """Mean and SD over seeds of each cascade stage, ages [lo,hi), by year+sex.

    Per-seed ratio-of-sums (cascade_agg), then summarised across seeds --
    the same order of operations as the cascade time-series figure above.
    """
    c = cascade_agg(df, lo, hi)
    stages = ["p_dx", "p_art_dx", "p_vls_art"]
    g = c.groupby(["year", "sex"])[stages]
    return g.mean().join(g.std(ddof=1), rsuffix="_sd").reset_index()


# %% 90-90-90 cascade by age, for the three PHIA survey years
# Columns = survey round, colour = sex, lines = model by 5-year band (mean
# over seeds), from age 0 so the perinatally infected children who drive the
# high values at 10-19 are visible. Four rows, one denominator per panel:
#   1-3: the three conditional 90s, with SHIMS (SHIMS_909090_coverage.csv,
#        circles, Wilson 95% CI) at its four non-overlapping bands -- the
#        horizontal bar is the band's age range. SHIMS has no 2011 round.
#   4:   ART coverage among ALL PLHIV, with
#        calibration_data/art_coverage_by_age_sex.csv (squares; 15-24, 25-34,
#        35-44 only) -- kept in its own row because its denominator differs
#        from row 2's.
# SHIMS2 (fieldwork 2016-17) sits in the column at shims2_year (2017 for the
# 027 configs, the year the model compares it at; 2016 otherwise).
# Model points are 5-year bands and SHIMS bands are 10-15 years wide, so read
# this as a shape check; Check B below compares like with like.
SHIMS2_YEAR = MODEL_KW.get("shims2_year", 2016)
YEARS = (2011, SHIMS2_YEAR, 2021)


def cascade_by_age(df, year):
    """Mean-over-seeds cascade fraction per stage, sex, and 5-year age band."""
    g = df[np.floor(df.timevec) == year]
    rows = []
    for sex in ("f", "m"):
        for lo, hi in AGE_BINS:
            inf = g[f"cascadeagesex.n_infected_{sex}_{lo}_{hi}"]
            dx = g[f"cascadeagesex.n_diagnosed_{sex}_{lo}_{hi}"]
            art = g[f"cascadeagesex.n_on_art_{sex}_{lo}_{hi}"]
            eff = g[f"cascadeagesex.n_effective_art_{sex}_{lo}_{hi}"]
            rows.append(dict(
                sex=sex, age_low=lo,
                p_dx=(dx / inf.replace(0, np.nan)).mean(),
                p_art_dx=(art / dx.replace(0, np.nan)).mean(),
                p_art_plhiv=(art / inf.replace(0, np.nan)).mean(),
                p_vls_art=(eff / art.replace(0, np.nan)).mean()))
    return pd.DataFrame(rows)


art_data = pd.read_csv(REPO / "calibration_data" / "art_coverage_by_age_sex.csv")
# 0=Male, 1=Female -- confirmed against plot_art_by_age_sex.py's own mapping
# for this exact file, not assumed (this project has two sex-coding
# conventions across different files; see experiments/021.../config.yaml).
art_data["sex"] = art_data.Gender.map({0: "m", 1: "f"})
art_data.loc[art_data.Year == 2016, "Year"] = SHIMS2_YEAR
shims_by_age = shims[shims.apply(lambda r: (r.age_lo, r.age_hi) in SHIMS_AGE_BANDS, axis=1)].copy()
shims_by_age.loc[shims_by_age.Year == 2017, "Year"] = SHIMS2_YEAR

STAGE_COLS = [("p_dx", "% diagnosed\n(of PLHIV)"),
              ("p_art_dx", "% on ART\n(of diagnosed)"),
              ("p_vls_art", "% suppressed\n(of on ART)"),
              ("p_art_plhiv", "% on ART\n(of ALL PLHIV)")]

fig, axes = plt.subplots(len(STAGE_COLS), len(YEARS), figsize=(5.2 * len(YEARS), 15.5),
                         sharex=True)
for row, (col, title) in enumerate(STAGE_COLS):
    for j, year in enumerate(YEARS):
        ax = axes[row, j]
        cba = cascade_by_age(df, year)
        for sex in ("f", "m"):
            s = cba[cba.sex == sex].sort_values("age_low")
            ax.plot(s.age_low + 2.5, s[col] * 100, "-o", color=CASC_COLORS[sex], ms=3)
        if col == "p_art_plhiv":
            for _, r in art_data[art_data.Year == year].iterrows():
                lo, hi = (int(x) for x in r.AgeBin.strip("[)").split(","))
                ax.errorbar((lo + hi) / 2, r.NationalARTPrevalence * 100,
                            xerr=(hi - lo) / 2,
                            yerr=[[(r.NationalARTPrevalence - r.lb) * 100],
                                  [(r.ub - r.NationalARTPrevalence) * 100]],
                            fmt="s", mfc=CASC_COLORS[r.sex], mec="black", mew=1.2,
                            ms=7, ecolor="black", elinewidth=1, capsize=3, zorder=6)
        else:
            d = shims_by_age[(shims_by_age.stage == col) & (shims_by_age.Year == year)]
            for _, r in d.iterrows():
                hi = min(r.age_hi, 65)  # 50+ drawn to the axis edge
                ax.errorbar((r.age_lo + hi) / 2, r.p * 100, xerr=(hi - r.age_lo) / 2,
                            yerr=[[(r.p - r.lb) * 100], [(r.ub - r.p) * 100]],
                            fmt="o", mfc=CASC_COLORS[r.sex], mec="black", mew=1.2,
                            ms=7, ecolor="black", elinewidth=1, capsize=3, zorder=6)
        if row == 0:
            ax.set_title(f"{year}" + (" (SHIMS2, fieldwork 2016-17)" if year == SHIMS2_YEAR
                                      else " (SHIMS3)" if year == 2021 else ""), fontsize=11)
        if j == 0:
            ax.set_ylabel(title, fontsize=9.5)
        ax.set_ylim(0, 105); ax.set_xlim(0, 65); ax.grid(alpha=0.3)
        if row == len(STAGE_COLS) - 1:
            ax.set_xlabel("age")

from matplotlib.lines import Line2D  # noqa: E402
handles = [Line2D([], [], color=CASC_COLORS["f"], marker="o", ms=3, label="Women, model"),
           Line2D([], [], color=CASC_COLORS["m"], marker="o", ms=3, label="Men, model"),
           Line2D([], [], ls="", marker="o", mfc="grey", mec="black", ms=7,
                  label="SHIMS 90-90-90 (rows 1-3), colour = sex; bars = age band, 95% CI"),
           Line2D([], [], ls="", marker="s", mfc="grey", mec="black", ms=7,
                  label="PHIA ART coverage among all PLHIV (row 4), colour = sex")]
fig.legend(handles=handles, loc="upper center", ncol=2, fontsize=9,
           bbox_to_anchor=(0.5, 0.995), frameon=False)
fig.suptitle(f"90-90-90 cascade by age, v4 check ({LABEL})", y=1.02, fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.965))
fig_casc_age = OUT / "cascade_by_age_v4_check.png"
fig.savefig(fig_casc_age, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved 90-90-90 by-age figure to {fig_casc_age}")

# %% Kasturi's 5-year bands -> 15-49 and 15+ aggregates
# Kasturi's file has fractions only, no population sizes, so aggregating
# needs weights from somewhere. We use the MODEL's own PLHIV per band
# (mean over seeds) to turn fractions into implied counts:
#   dx_b  = aware_b * PLHIV_b          art_b = art_b|aware * dx_b
#   p_dx  = sum(dx_b) / sum(PLHIV_b)   p_art_dx = sum(art_b) / sum(dx_b)
# So the aggregate data line inherits the model's age/sex distribution of
# PLHIV. If the model's PLHIV age structure is off, this line moves with it --
# the 5-year-band figure further down is the comparison free of that.
# Kasturi stops at 80-84, so its "15+" is 15-84; the model's 15+ is 15-99.
def kasturi_aggregate(kas, df, lo, hi):
    plhiv = []
    for sex in ("f", "m"):
        for a, b in AGE_BINS:
            if a >= lo and b <= hi:
                col = f"cascadeagesex.n_infected_{sex}_{a}_{b}"
                m = df.groupby("timevec")[col].mean()
                plhiv.append(pd.DataFrame({"year": m.index.astype(int), "sex": sex,
                                           "age_low": a, "plhiv": m.values}))
    k = kas.merge(pd.concat(plhiv), on=["year", "sex", "age_low"], how="inner")
    k["dx"] = k.p_dx * k.plhiv
    k["art"] = k.p_art_dx * k.dx
    both = k.assign(sex="all")
    s = pd.concat([k, both]).groupby(["year", "sex"])[["plhiv", "dx", "art"]].sum()
    return pd.DataFrame({"p_dx": s.dx / s.plhiv,
                         "p_art_dx": s.art / s.dx}).reset_index()


kas_agg = {grp: kasturi_aggregate(kas, df, lo, hi) for (lo, hi), grp in AGG_GROUPS.items()}
mod_agg = {grp: model_summary(df, lo, hi) for (lo, hi), grp in AGG_GROUPS.items()}

# %% Check A: 90-90-90 time series vs SHIMS (15-49 and 15+)
# Rows = age group, columns = stage. One point per sex per survey round.
# SHIMS publishes no both-sexes row; the grey model line has no data to
# compare against.
fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
for i, grp in enumerate(AGG_GROUPS.values()):
    for j, stage in enumerate(STAGE_TITLES):
        ax = axes[i, j]
        lo, hi = {v: k for k, v in AGG_GROUPS.items()}[grp]
        c = cascade_agg(df, lo, hi)
        for sex in ("f", "m", "all"):
            _band(ax, c[c.sex == sex], stage, CASC_COLORS[sex], CASC_NAMES[sex])
        d = shims[(shims.stage == stage) & (shims["Age group"] == grp)]
        for _, r in d.iterrows():
            ax.errorbar(r.Year, r.p * 100,
                        yerr=[[(r.p - r.lb) * 100], [(r.ub - r.p) * 100]],
                        fmt="o", mfc=CASC_COLORS[r.sex],
                        mec="black", mew=1.2,
                        ms=8, ecolor="black", capsize=4, zorder=6)
        ax.set_title(f"Ages {grp}: {STAGE_TITLES[stage]}", fontsize=10)
        if j == 0:
            ax.set_ylabel("%")
        ax.set_xlim(2000, STOP); ax.set_ylim(0, 105); ax.grid(alpha=0.3)
        if i == 1:
            ax.set_xlabel("year")
axes[0, 0].errorbar([], [], fmt="o", mfc="grey", mec="black", ms=8, capsize=4,
                    label="SHIMS (Wilson 95% CI, unweighted n)")
axes[0, 0].legend(fontsize=7)
fig.suptitle("90-90-90 time series vs SHIMS, ages 15-49 (top row) and 15+ "
             f"(bottom row), v4 check ({LABEL})\n"
             "line = mean of 10 seeds, band = +/- 1 SD", y=1.01, fontsize=11)
fig_path = OUT / "cascade_timeseries_vs_shims_v4_check.png"
fig.savefig(fig_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved 90-90-90 time series vs SHIMS to {fig_path}")

# %% Check B: 90-90-90 by SHIMS age band, one column per (sex, survey round)
# Rows = cascade stage, columns = sex x survey year (women first, then men),
# so each panel holds a single sex and the two sexes never share a facet.
# The x-axis is CATEGORICAL: SHIMS's own bands (15-24, 25-34, 35-49, 50+),
# with the model's 5-year bands summed into exactly those bands before taking
# ratios (model_summary), so model and data describe the same people.
# Model (open circle, +/- 1 SD over seeds) and SHIMS (filled square, Wilson
# 95% CI) are drawn at the SAME x -- the category's own tick, no offset.
# SHIMS 2017 has no 50+ row, so that category shows the model only.
BAND_LABELS = [f"{lo}-{hi - 1}" if hi < 100 else f"{lo}+" for lo, hi in SHIMS_AGE_BANDS]
mod_bands = {band: model_summary(df, *band) for band in SHIMS_AGE_BANDS}
COLS = [(sex, year) for sex in ("f", "m") for year in SHIMS_YEARS]

fig, axes = plt.subplots(len(STAGE_TITLES), len(COLS),
                         figsize=(4.2 * len(COLS), 11), sharex=True, sharey=True)
for j, (sex, year) in enumerate(COLS):
    for i, stage in enumerate(STAGE_TITLES):
        ax = axes[i, j]
        for x, ((lo, hi), m) in enumerate(mod_bands.items()):
            r = m[(m.year == year) & (m.sex == sex)].iloc[0]
            ax.errorbar(x, r[stage] * 100,
                        yerr=r[f"{stage}_sd"] * 100, fmt="o", mfc="white",
                        mec=CASC_COLORS[sex], ecolor=CASC_COLORS[sex], mew=1.5,
                        ms=8, capsize=3, zorder=5)
            d = shims[(shims.Year == year) & (shims.stage == stage)
                      & (shims.sex == sex) & (shims.age_lo == lo)
                      & (shims.age_hi == hi)]
            for _, d in d.iterrows():
                ax.errorbar(x, d.p * 100,
                            yerr=[[(d.p - d.lb) * 100], [(d.ub - d.p) * 100]],
                            fmt="s", mfc=CASC_COLORS[sex], mec="black", mew=1.2,
                            ms=7, ecolor="black", capsize=3, zorder=6)
        if i == 0:
            ax.set_title(f"{CASC_NAMES[sex]}, {year}", fontsize=11)
        if j == 0:
            ax.set_ylabel(STAGE_TITLES[stage], fontsize=9.5)
        if i == len(STAGE_TITLES) - 1:
            ax.set_xlabel("age group")
        ax.set_xticks(range(len(BAND_LABELS)), BAND_LABELS)
        # Zoomed y-axis: SHIMS strata run ~55-100% (lowest CI ~42%),
        # unreadable on 0-105.
        ax.set_xlim(-0.5, len(BAND_LABELS) - 0.5); ax.set_ylim(40, 102)
        ax.grid(alpha=0.3, axis="y")
ax = axes[0, 0]
ax.plot([], [], "o", mfc="white", mec="grey", mew=1.5, label="model (+/- 1 SD)")
ax.plot([], [], "s", mfc="grey", mec="black", label="SHIMS (Wilson 95% CI)")
ax.legend(fontsize=7, loc="lower right")
fig.suptitle(f"90-90-90 by age group vs SHIMS, v4 check ({LABEL})",
             y=1.01, fontsize=11)
fig_path = OUT / "cascade_by_age_vs_shims_v4_check.png"
fig.savefig(fig_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved 90-90-90 by age vs SHIMS to {fig_path}")

# %% Check C: 90-90-90 time series vs Kasturi, by age group and sex
# Rows = cascade stage, columns = age group x sex (one sex per panel).
# Kasturi covers only the first two stages; the VLS row is model-only, shown
# so all three 90s can be read off one figure. Kasturi's aggregate line
# depends on model PLHIV weights -- see kasturi_aggregate's comment.
# Dashed = Kasturi, solid + band = model.
STAGES_C = list(STAGE_TITLES)
COLS = [(grp, sex) for grp in AGG_GROUPS.values() for sex in ("f", "m")]
GRP_BOUNDS = {grp: lohi for lohi, grp in AGG_GROUPS.items()}
fig, axes = plt.subplots(len(STAGES_C), len(COLS), figsize=(4.2 * len(COLS), 11),
                         sharex=True, sharey=True)
for j, (grp, sex) in enumerate(COLS):
    c = cascade_agg(df, *GRP_BOUNDS[grp])
    k = kas_agg[grp][kas_agg[grp].sex == sex]
    for i, stage in enumerate(STAGES_C):
        ax = axes[i, j]
        _band(ax, c[c.sex == sex], stage, CASC_COLORS[sex], "model")
        if stage in k:
            ax.plot(k.year, k[stage] * 100, "--o", color=CASC_COLORS[sex],
                    ms=3, lw=1.3, label="Kasturi")
        else:
            ax.text(0.97, 0.05, "no Kasturi estimate", transform=ax.transAxes,
                    ha="right", fontsize=8, color="grey")
        if i == 0:
            ax.set_title(f"{CASC_NAMES[sex]}, ages {grp}", fontsize=11)
        if j == 0:
            ax.set_ylabel(STAGE_TITLES[stage], fontsize=9.5)
        if i == len(STAGES_C) - 1:
            ax.set_xlabel("year")
        ax.set_xlim(1995, STOP); ax.set_ylim(0, 105); ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=7, loc="lower right")
fig.suptitle(f"90-90-90 time series vs Kasturi's modelled estimates, v4 check ({LABEL})\n"
             "Kasturi 5-year bands aggregated with MODEL PLHIV weights; "
             "Kasturi '15+' = 15-84", y=1.02, fontsize=11)
fig_path = OUT / "cascade_timeseries_vs_kasturi_v4_check.png"
fig.savefig(fig_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved 90-90-90 time series vs Kasturi to {fig_path}")

# %% Check D: first two 90s by 5-year age band vs Kasturi
# Same 5-year bands on both sides, so no weighting assumption here.
# Model = mean over seeds of per-seed ratios (cascade_by_age, as above).
fig, axes = plt.subplots(2, len(YEARS), figsize=(5.2 * len(YEARS), 8.5),
                         sharex=True, sharey=True)
for j, year in enumerate(YEARS):
    cba = cascade_by_age(df, year)
    for i, stage in enumerate(("p_dx", "p_art_dx")):
        ax = axes[i, j]
        for sex in ("f", "m"):
            s = cba[(cba.sex == sex) & (cba.age_low >= 15) & (cba.age_low < 85)]
            ax.plot(s.age_low + 2.5, s[stage] * 100, "-o", color=CASC_COLORS[sex],
                    ms=3, label=f"{CASC_NAMES[sex]} model")
            k = kas[(kas.year == year) & (kas.sex == sex)].sort_values("age_low")
            ax.plot(k.age_low + 2.5, k[stage] * 100, "--s", color=CASC_COLORS[sex],
                    ms=4, mfc="white", lw=1.3, label=f"{CASC_NAMES[sex]} Kasturi")
        if i == 0:
            ax.set_title(f"{year}", fontsize=11)
        if j == 0:
            ax.set_ylabel(f"{STAGE_TITLES[stage]}")
        if i == 1:
            ax.set_xlabel("age")
        ax.set_xlim(13, 87); ax.set_ylim(0, 105); ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=7)
fig.suptitle("First two 90s by age vs Kasturi's modelled estimates, "
             f"v4 check ({LABEL})", y=1.01, fontsize=11)
fig_path = OUT / "cascade_by_age_vs_kasturi_v4_check.png"
fig.savefig(fig_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved first-two-90s by age vs Kasturi to {fig_path}")

# %% Load SHIMS incidence (DTC_check/SHIMS_Datasets/SHIMS_incidence_data.csv)
# Published annual incidence, % per year, with 95% CIs, by sex (Male, Female,
# and "total" = both sexes) and age group, for SHIMS2 (labelled 2017 in the
# file, as in the 90-90-90 file) and SHIMS3 (2021). Comparison only -- not a
# fitting target; the 2021 rows are CLAUDE.md's validation hold-out year.
# Same caveat as the calibration targets: recency-assay false-recency bias
# inflates older age bands (see incidence_construction.py).
SEX_OF_GENDER_INC = {**SEX_OF_GENDER, "total": "all"}
shims_inc = pd.read_csv(DATA_DIR / "SHIMS_incidence_data.csv", encoding="utf-8-sig")
shims_inc["sex"] = shims_inc.Gender.str.lower().map(SEX_OF_GENDER_INC)
assert shims_inc.sex.notna().all(), "unrecognised Gender in SHIMS incidence"
shims_inc[["age_lo", "age_hi"]] = shims_inc["Age group"].apply(
    lambda a: pd.Series(parse_age_group(a)))
shims_inc = shims_inc.rename(columns={"Percentage annual incidence": "inc",
                                      "CI_lower": "lb", "CI_upper": "ub"})
SHIMS_INC_YEARS = sorted(shims_inc.Year.unique())


def model_incidence(df, sex, lo, hi):
    """Per-seed annual incidence % over ages [lo,hi), long-form (year, seed, inc).

    Same definition as standard_figures._inc_series -- new infections over
    SUSCEPTIBLES (alive minus infected), from PopByAgeSex's 5-year bands --
    extended to sex="all" by summing both sexes' counts before the ratio.
    The STOP year is dropped: the sim ends on its first timestep, so its
    new-infection count covers a fraction of a year and reads as ~0.
    """
    sexes = ("f", "m") if sex == "all" else (sex,)
    new = alv = inf = 0
    for s in sexes:
        for a, b in AGE_BINS:
            if a >= lo and b <= hi:
                new = new + df[f"popagesex.new_infections_{s}_{a}_{b}"]
                alv = alv + df[f"popagesex.n_alive_{s}_{a}_{b}"]
                inf = inf + df[f"popagesex.n_infected_{s}_{a}_{b}"]
    susc = alv - inf
    out = pd.DataFrame({"year": df.timevec.astype(int), "seed": df.seed,
                        "inc": 100.0 * new / susc.where(susc > 0)})
    return out[out.year < STOP]


def model_incidence_at(df, sex, lo, hi, year):
    """(mean, SD) over seeds of incidence % in one year."""
    v = model_incidence(df, sex, lo, hi)
    v = v[v.year == year].inc.dropna()
    return v.mean(), v.std(ddof=1)


# %% Incidence check 1: incidence by SHIMS age group vs SHIMS
# Same layout as the 90-90-90 by-age check: one row, one column per
# (sex, survey year), women first then men; categorical x-axis of SHIMS's own
# age groups, with the model's 5-year bands summed into exactly those groups.
# Model (open circle, +/- 1 SD over seeds) and SHIMS (filled square, 95% CI)
# share the category's tick -- no offset. SHIMS2 has no 50+ row.
COLS = [(sex, year) for sex in ("f", "m") for year in SHIMS_INC_YEARS]
fig, axes = plt.subplots(1, len(COLS), figsize=(4.2 * len(COLS), 4.6),
                         sharey=True)
for ax, (sex, year) in zip(axes, COLS):
    for x, (lo, hi) in enumerate(SHIMS_AGE_BANDS):
        m, sd = model_incidence_at(df, sex, lo, hi, year)
        ax.errorbar(x, m, yerr=sd, fmt="o", mfc="white", mec=CASC_COLORS[sex],
                    ecolor=CASC_COLORS[sex], mew=1.5, ms=8, capsize=3, zorder=5)
        d = shims_inc[(shims_inc.Year == year) & (shims_inc.sex == sex)
                      & (shims_inc.age_lo == lo) & (shims_inc.age_hi == hi)]
        for _, d in d.iterrows():
            ax.errorbar(x, d.inc, yerr=[[d.inc - d.lb], [d.ub - d.inc]],
                        fmt="s", mfc=CASC_COLORS[sex], mec="black", mew=1.2,
                        ms=7, ecolor="black", capsize=3, zorder=6)
    ax.set_title(f"{CASC_NAMES[sex]}, {year}", fontsize=11)
    ax.set_xticks(range(len(BAND_LABELS)), BAND_LABELS)
    ax.set_xlim(-0.5, len(BAND_LABELS) - 0.5); ax.set_ylim(bottom=0)
    ax.set_xlabel("age group"); ax.grid(alpha=0.3, axis="y")
axes[0].set_ylabel("HIV incidence (% per year)")
axes[0].plot([], [], "o", mfc="white", mec="grey", mew=1.5, label="model (+/- 1 SD)")
axes[0].plot([], [], "s", mfc="grey", mec="black", label="SHIMS (95% CI)")
axes[0].legend(fontsize=7, loc="upper right")
fig.suptitle(f"HIV incidence by age group vs SHIMS, v4 check ({LABEL})",
             y=1.03, fontsize=11)
fig_path = OUT / "incidence_by_age_vs_shims_v4_check.png"
fig.savefig(fig_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved incidence by age vs SHIMS to {fig_path}")

# %% Incidence check 2: incidence time series vs SHIMS (15-49 and 15+)
# Rows = age group, columns = women / men / both sexes (SHIMS's "total").
# Unlike incidence_fit_v4_check.png, which compares against the calibration
# targets (SHIMS1 2011 18-49, SHIMS2 2016 15-49), this compares against the
# SHIMS2/SHIMS3 file above. Line = mean of seeds, band = +/- 1 SD.
inc_cols = ("f", "m", "all")
fig, axes = plt.subplots(len(AGG_GROUPS), len(inc_cols),
                         figsize=(5.2 * len(inc_cols), 8.5), sharex=True, sharey=True)
for i, (lohi, grp) in enumerate(AGG_GROUPS.items()):
    for j, sex in enumerate(inc_cols):
        ax = axes[i, j]
        s = model_incidence(df, sex, *lohi).groupby("year").inc
        mean, sd = s.mean(), s.std(ddof=1).fillna(0.0)
        ax.plot(mean.index, mean.values, color=CASC_COLORS[sex], lw=2)
        ax.fill_between(mean.index, (mean - sd).clip(lower=0).values,
                        (mean + sd).values, color=CASC_COLORS[sex], alpha=0.15, lw=0)
        d = shims_inc[(shims_inc["Age group"] == grp) & (shims_inc.sex == sex)]
        for _, r in d.iterrows():
            ax.errorbar(r.Year, r.inc, yerr=[[r.inc - r.lb], [r.ub - r.inc]],
                        fmt="s", mfc=CASC_COLORS[sex], mec="black", mew=1.2,
                        ms=8, ecolor="black", capsize=4, zorder=6)
        ax.set_title(f"{CASC_NAMES[sex]}, ages {grp}", fontsize=11)
        ax.set_xlim(2000, STOP); ax.set_ylim(bottom=0); ax.grid(alpha=0.3)
        if j == 0:
            ax.set_ylabel("HIV incidence (% per year)")
        if i == len(AGG_GROUPS) - 1:
            ax.set_xlabel("year")
axes[0, 0].plot([], [], "-", color="grey", lw=2, label="model (+/- 1 SD)")
axes[0, 0].errorbar([], [], fmt="s", mfc="grey", mec="black", ms=8, capsize=4,
                    label="SHIMS (95% CI)")
axes[0, 0].legend(fontsize=7, loc="upper right")
fig.suptitle(f"HIV incidence over time vs SHIMS, v4 check ({LABEL})",
             y=1.01, fontsize=11)
fig_path = OUT / "incidence_timeseries_vs_shims_v4_check.png"
fig.savefig(fig_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved incidence time series vs SHIMS to {fig_path}")

# %% Intervention counts per year: ART initiations, VMMC, oral and LA PrEP
# Four panels on a shared x-axis, one point per year (mean over seeds, bar =
# +/- 1 SD), no connecting lines. All counts are national (scaled by
# pop_scale) and are yearly sums of the per-step "new_" results:
#   ART  -- hiv.new_agents_on_art: agents whose ART start fell in that step.
#           Includes RE-starts after the ART coverage correction stops
#           someone, so it can exceed first-ever initiations.
#   VMMC -- vmmc.new_circumcisions (programme circumcisions topping each age
#           band up to its data/vmmc_coverage.csv target).
#   PrEP -- prep_oral / prep_la .new_prep (027 only): initiations, which
#           should reproduce data/prep_initiations.csv (overlaid, black x).
# The STOP year is dropped: it holds a single timestep.
COUNT_PANELS = [("hiv.new_agents_on_art", "ART initiations"),
                ("vmmc.new_circumcisions", "VMMC (programme circumcisions)"),
                ("prep_oral.new_prep", "Oral PrEP initiations"),
                ("prep_la.new_prep", "LA-PrEP initiations")]
PREP_DATA = {"prep_oral.new_prep": "OralTDF", "prep_la.new_prep": "LA-PREP_LEN"}
prep_in = pd.read_csv(REPO / "data" / "prep_initiations.csv", encoding="utf-8-sig")

fig, axes = plt.subplots(len(COUNT_PANELS), 1, figsize=(10, 11), sharex=True)
for ax, (col, title) in zip(axes, COUNT_PANELS):
    if col in df:
        s = df[df.timevec < STOP].groupby("timevec")[col]
        ax.errorbar(s.mean().index, s.mean().values, yerr=s.std(ddof=1).values,
                    fmt="o", ms=4, color="#2c6fbb", capsize=2, elinewidth=0.8,
                    label="model (mean of seeds, +/- 1 SD)")
        if col in PREP_DATA:
            d = prep_in[(prep_in.PrEP_type == PREP_DATA[col]) & (prep_in.Year < STOP)]
            ax.plot(d.Year, d.Number, "x", color="black", ms=7, mew=1.5,
                    label="programme data")
        ax.legend(fontsize=8, loc="upper left")
    else:
        ax.text(0.5, 0.5, "not in this model configuration", transform=ax.transAxes,
                ha="center", va="center", color="grey")
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_ylabel("per year"); ax.set_ylim(bottom=0); ax.grid(alpha=0.3)
axes[-1].set_xlabel("year"); axes[-1].set_xlim(1995, STOP)
fig.suptitle(f"Intervention counts per year, national, v4 check ({LABEL})",
             y=1.0, fontsize=11)
fig.tight_layout()
fig_path = OUT / "intervention_counts_v4_check.png"
fig.savefig(fig_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"Saved intervention counts to {fig_path}")

# %% Numbers behind checks A-D and the SHIMS incidence checks, as one table
# One row per (dataset, stage, year, sex, age group): data vs model mean/SD.
# Incidence rows are in % per year; cascade rows are proportions.
# This is the table to read when the figures are too compressed near 90-100%.
rows = []
for _, d in shims.iterrows():
    m = model_summary(df, d.age_lo, d.age_hi)
    r = m[(m.year == d.Year) & (m.sex == d.sex)].iloc[0]
    rows.append(dict(dataset="SHIMS", stage=d.stage, year=d.Year, sex=d.sex,
                     age_group=d["Age group"], data=d.p, data_lb=d.lb,
                     data_ub=d.ub, data_n=d.Number,
                     model=r[d.stage],
                     model_sd=r[f"{d.stage}_sd"]))
for grp in AGG_GROUPS.values():
    k = kas_agg[grp].merge(mod_agg[grp], on=["year", "sex"], suffixes=("", "_mod"))
    for stage in ("p_dx", "p_art_dx"):
        for _, r in k.iterrows():
            rows.append(dict(dataset="Kasturi (model-PLHIV-weighted)", stage=stage,
                             year=r.year, sex=r.sex, age_group=grp, data=r[stage],
                             model=r[f"{stage}_mod"], model_sd=r[f"{stage}_sd"]))
for year in YEARS:
    cba = cascade_by_age(df, year)
    k = kas[kas.year == year].merge(cba, on=["sex", "age_low"], suffixes=("", "_mod"))
    for stage in ("p_dx", "p_art_dx"):
        for _, r in k.iterrows():
            rows.append(dict(dataset="Kasturi", stage=stage, year=year, sex=r.sex,
                             age_group=r.agecat, data=r[stage],
                             model=r[f"{stage}_mod"]))
for _, d in shims_inc.iterrows():
    m, sd = model_incidence_at(df, d.sex, d.age_lo, d.age_hi, d.Year)
    rows.append(dict(dataset="SHIMS incidence (% per year)", stage="incidence",
                     year=d.Year, sex=d.sex, age_group=d["Age group"], data=d.inc,
                     data_lb=d.lb, data_ub=d.ub, model=m, model_sd=sd))
comparison = pd.DataFrame(rows)
comparison["model_minus_data"] = comparison.model - comparison.data
table_path = OUT / "cascade_vs_data_v4_check.csv"
comparison.to_csv(table_path, index=False, float_format="%.4f")
print(f"Saved model-vs-data 90-90-90 table to {table_path}")
with pd.option_context("display.float_format", "{:.3f}".format,
                       "display.width", 160):
    print(comparison[(comparison.dataset == "SHIMS")
                     & comparison.age_group.isin(AGG_GROUPS.values())]
          .drop(columns=["dataset"]).to_string(index=False))

# %%
