# Eswatini HIVsim baseline — starsim 3.7.2 snapshot

A frozen copy of the calibrated Eswatini model as it runs on **starsim 3.7.2 /
stisim 1.7.0**. It's the reference point for porting the model to Starsim v4.

## Provenance

- **Source:** [`hivsim_eswatini`](https://github.com/starsimhub/hivsim_eswatini),
  branch `exp/027-cascade-prep-recal`, commit `068dd83` (2026-11-19).
- **Model:** the `"027d"` configuration of `DTC_check/build_and_check_calibration.py`:
  - Adam Akullian's 024 point fit (design row 868; 7 transmission parameters);
  - Daniel Citron's exp 027 inputs (below).
- **Dependencies:** starsim `main` @ `2a99fc85` (3.7.2) and stisim `main` @
  `3730ddc` (1.7.0), from the editable checkouts in `IDM/`. Python 3.13,
  conda env `hivsim-dev`.

## What's in the model (exp 027 additions)

| Component | Setting | Code |
|---|---|---|
| General-population HIV testing | separate start year by sex; level by sex × age cohort (0–24, 25–34, 35–49, 50+; children share 0–24), fitted to SHIMS2 % diagnosed | `interventions.py`, `params/calib2d_chosen.json` |
| ANC testing | logistic: starts 1992.5, midpoint 2003, ~max (90%) by 2010 | `interventions.anc_test_prob_logistic` |
| PMTCT | prophylaxis for diagnosed pregnant women not on ART, 2003 → 2010 | `pmtct.py` |
| ART | forced coverage among PLHIV by sex × the same cohorts, from SHIMS 90-90-90 (0–14 take 15–24) | `data/art_coverage_cohorts.csv` |
| Viral suppression on ART | logistic through SHIMS2 (2017) and SHIMS3 (2021) | `vls_construction.to_vls_logistic`, `vls_stock_target.py` |
| VMMC | forced circumcision coverage by age; SHIMS2 medical + non-medical totals at 2017, SHIMS3 at 2021 | `data/vmmc_coverage.csv`, `interventions.SHIMS2_VMMC` |
| PrEP | oral TDF (PrEPWatch counts, 75% efficacy) and LA-PrEP (MOH counts, 92%); 1-year courses; annual initiations; FSW, clients, >1 partner, no stable partner | `prep_initiations.py`, `data/prep_initiations.csv` |

Full rationale, fits and caveats are in hivsim_eswatini's
`experiments/027_cascade_prep_recal/README.md`.

## Use

```python
from baseline import make_baseline_sim
sim = make_baseline_sim(seed=1)          # 1985-2026; kwargs pass to run_sims.make_sim
sim.run()
```

`python baseline.py` runs seed 1 and compares prevalence, new infections and
deaths with `reference/baseline_027d_10seeds.df` (seeds 1–10, annual
`sim.to_df(resample="year", use_years=True, sep=".")`). Seed 1 reproduces it
exactly on the versions above. That reference is what a v4 port should be
checked against.

## Files

- `baseline.py` — the calibrated parameters and `make_baseline_sim()`.
- `run_sims.py`, `interventions.py`, `analyzers.py`, `pmtct.py`,
  `prep_initiations.py`, `vls_construction.py`, `vls_stock_target.py` —
  copied unchanged from the source commit.
- `data/` — exactly the files a run reads (traced from a full run).
- `params/calib2d_chosen.json` — the fitted testing inputs.
- `reference/` — 10-seed output of this configuration.

The code reads `data/` relative to the working directory; `baseline.py`
changes into this folder on import.
