"""The calibrated Eswatini HIVsim baseline (starsim 3.7.2 / stisim 1.7.0).

Copied from hivsim_eswatini, branch exp/027-cascade-prep-recal, commit
068dd83: the "027d" configuration of DTC_check/build_and_check_calibration.py.
That is Adam Akullian's 024 point fit (design row 868, 7 transmission
parameters) plus Daniel Citron's exp 027 inputs: age-structured testing and
forced ART, logistic ANC testing, PMTCT, logistic viral suppression, SHIMS2
VMMC totals, and PrEP from programme initiation counts. See README.md.

    from baseline import make_baseline_sim
    sim = make_baseline_sim(seed=1)
    sim.run()

Run directly to simulate one seed and print a quick check against the
reference output in reference/.
"""

import os
from pathlib import Path

import sciris as sc

HERE = Path(__file__).resolve().parent
os.chdir(HERE)  # run_sims.py and interventions.py read data/ relative to cwd

from run_sims import make_sim  # noqa: E402

# 024 design row 868 (hivsim_eswatini experiments/026_scenarios/run.py::BEST)
BEST = dict(beta_m2f=0.014161835611057863, rel_beta_f2m=0.3097261030208422,
            s_f_young=2.224000314771824, age_gap_shift=0.0699801603950227,
            age_gap_sd_mult=1.5001144712607821, prop_f0=0.608751696061895,
            prop_m0=0.4742150061435999)

# Fixed (non-calibrated) network structure, as 024/026
AGE_DIFF_BASE = {"teens": [(7, 3), (6, 3), (5, 1)],
                 "young": [(8, 3), (7, 3), (5, 2)],
                 "adult": [(8, 3), (7, 3), (5, 2)]}
CONC_BASE = {"f1_conc": 0.15, "f2_conc": 0.25, "m1_conc": 0.15, "m2_conc": 0.5}
REL_INIT_PREV = 0.2


def build_pars(p=BEST):
    """Map the 7 calibrated parameters onto make_sim's hiv_pars / network_pars."""
    hiv_pars = dict(
        beta_m2f=p["beta_m2f"], rel_beta_f2m=p["rel_beta_f2m"],
        rel_init_prev=REL_INIT_PREV,
        rel_sus_age=[(15, 25, "f", p["s_f_young"]),
                     (25, 50, "f", 1.0), (15, 50, "m", 1.0)])
    shift, sd_mult = p["age_gap_shift"], p["age_gap_sd_mult"]
    network_pars = dict(
        prop_f0=p["prop_f0"], prop_m0=p["prop_m0"],
        age_diff_pars={g: [(max(m + shift, 1.0), max(s * sd_mult, 0.2))
                           for m, s in v] for g, v in AGE_DIFF_BASE.items()},
        **CONC_BASE)
    return hiv_pars, network_pars


def cascade_kw(path=HERE / "params" / "calib2d_chosen.json"):
    """make_sim keyword arguments for the exp 027d cascade, PrEP and VMMC inputs."""
    c = sc.loadjson(path)
    return dict(
        testing_pars=dict(
            gp_start_m=c["gp_start_m"], gp_start_f=c["gp_start_f"],
            rel_test_age={s: [c[f"rel_{s}_{k}"] for k in ("0", "25", "35", "50")] for s in "mf"},
            anc=True, anc_logistic=True),
        prep=True, shims2_year=2017, pmtct=True, art_vls_coverage="phia_logistic",
        art_coverage="cohorts", vmmc_shims2=True)


def make_baseline_sim(seed=1, start=1985, stop=2026, verbose=0, **kwargs):
    """The calibrated baseline. kwargs pass through to run_sims.make_sim
    (e.g. analyzers=, extra_interventions=) and override the defaults here."""
    hiv_pars, network_pars = build_pars()
    kw = dict(seed=seed, start=start, stop=stop, verbose=verbose,
              hiv_pars=hiv_pars, network_pars=network_pars, **cascade_kw())
    kw.update(kwargs)
    return make_sim(**kw)


if __name__ == "__main__":
    sim = make_baseline_sim(seed=1)
    sim.run()
    df = sim.to_df(resample="year", use_years=True, sep=".")
    # The reference is the Starsim 3.7.2 run (model_3_7_2): its time column is
    # 'timevec' (v4: 'time'), and v3 result names are mapped to v4's here. It is
    # a different random stream, and v4 records stocks at the start of a step,
    # so compare over several seeds, not run for run.
    ref = sc.loadobj(HERE / "reference" / "baseline_027d_10seeds.df")
    ref = ref.groupby("timevec").mean(numeric_only=True)
    pairs = {"hiv.prevalence": "hiv.prevalence", "hiv.new_infections": "hiv.new_infections",
             "hiv.new_deaths": "hiv.new_deaths", "hiv.n_on_art": "hiv.on_art"}
    for v3_col, v4_col in pairs.items():
        if v4_col in df and v3_col in ref:
            print(v4_col, {y: (round(float(df.loc[df.time == y, v4_col].iloc[0]), 4),
                               round(float(ref.loc[y, v3_col]), 4))
                           for y in (2007, 2017, 2021)}, "(this run, 3.7.2 reference mean)")
