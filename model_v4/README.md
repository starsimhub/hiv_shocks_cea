# Eswatini HIVsim baseline — Starsim v4 port (in progress)

The `model_3_7_2` snapshot being ported to **Starsim 4.0.0rc5 / stisim `v4`**
(conda env `hivsim-v4`; stisim from the `IDM/stisim_v4` worktree, branch
`v4` @ `9fe22a2`). `model_3_7_2/` is the frozen reference it is checked against.

## Port status (2026-11-19)

- **Runs end to end on v4.** Changes so far: Starsim's migration script
  (`ss.MaternalNet` → `ss.PrenatalNet`, `ss.lognorm_ex` → `ss.lognormal`,
  pregnancy state names, `sim.pop_scale`); ANC first-trimester eligibility
  rewritten for v4's pregnancy module (time in `P` under 13 weeks); `baseline.py`'s
  comparison maps v3 result names to v4's.
- **Transmission recalibrated on v4** (`calibrate_transmission.py` →
  `params/v4_transmission.json`; now the default, `baseline.BEST`; Adam's row
  868 is `baseline.ROW868`). Row 868's parameters, fitted on 3.7.2, gave too large
  an epidemic on v4 (women 15–49 in 2016: 47% vs PHIA 36%), as stisim's v4
  changelog predicts. The same 7 parameters were refitted with optuna TPE (405
  trials × 2 seeds, best 10 re-run on 8 fresh seeds) against PHIA prevalence by
  age and sex, 15–49 by sex, women 20–39 and the SHIMS incidence aggregates,
  with the 027d cascade inputs fixed and 2021 held out. Result (10 seeds):
  prevalence-by-age MAE 0.035, bias −0.009, 32/54 within PHIA's CI (3.7.2
  fit: 0.038, +0.008, ~28/54); with 2021, 44/74. Largest parameter changes:
  young women's extra susceptibility `s_f_young` 2.22 → 1.03; narrower partner
  age gaps (`age_gap_sd_mult` 1.50 → 0.70, `age_gap_shift` +0.07 → −0.65);
  more men in the lowest risk group (`prop_m0` 0.47 → 0.61); `rel_beta_f2m`
  0.31 → 0.42. The testing rates (027d) were not refitted.
- **Analyzers ported** (`analyzers.py`): `PopByAgeSex` and `CascadeByAgeSex`
  (moved here from hivsim_eswatini's check script) record counts at the start of
  each step and new infections at its end, with counts averaged over each year
  when resampled, as in 3.7.2; `hiv_epi` keeps its end-of-step timing. All
  results declare `kind=`.
- **Calibration check:** `check_calibration.py` runs 10 seeds and draws the
  hivsim_eswatini check figures (same plotting code) into
  `../results/model_v4_check/`, against `calibration_data/` (prevalence,
  incidence and ART-coverage targets, 2021 hold-out, copied from
  hivsim_eswatini) and the SHIMS / Kasturi data in `data/`. With row 868 on
  v4: MAE 0.068, bias +0.052; with the v4 refit: 0.035, −0.009.

The timeline and component tables below describe the model as inherited from
3.7.2 (hivsim_eswatini commit `068dd83`; full provenance in `model_3_7_2/README.md`).

## Timeline of model updates

| When | Who | What changed |
|---|---|---|
| Mar 2026 | Robyn Stuart | Initial repository: the Eswatini model, data pipeline and plotting scripts. |
| Apr–May 2026 (exp 001–007) | Adam Akullian | Model development: ART coverage forced by age, sex and year; DHS sexual-debut ages; women:men incidence balance; VMMC as a circumcision-coverage target (SDHS 2007, SHIMS2, SHIMS3); Bellan acute-infection parameters; initial prevalence. |
| May–Aug 2026 (exp 008–016) | Adam Akullian | Calibration targets set: PHIA/SHIMS prevalence by age and sex (2007, 2011, 2016), UNAIDS AIDS deaths, UNAIDS age structure; 2021 held out for validation. VMMC fix kept in the repo; double-counted mortality removed. |
| Aug–Sep 2026 (exp 017–023) | Adam Akullian | stisim upgraded to 1.5.11; stisim's built-in PrEP default removed (it had silently added PrEP for FSW from 2004); age-dependent survival on ART; viral suppression on ART taken from SHIMS and held as a target; 7 transmission and network parameters chosen for calibration. |
| Sep 2026 (exp 024) | Adam Akullian | History matching wave 1: 1,000 points, 7 parameters, 48 targets. **Design row 868** fits all 48 within 3σ and becomes the calibrated point fit used here. |
| Sep 2026 (exp 025–026) | Adam Akullian | Wave 2 designed (not run). Decision scenarios (long-acting PrEP vs ART scale-up) on row 868. |
| Oct 2026 | Daniel Citron | Reviewed Adam's calibration with a check dashboard: prevalence by age, the 90-90-90 cascade, incidence, all against SHIMS and Kasturi's modelled series. Found that the model's diagnosis rates came from hand-set testing ramps, not data, and were far too high before 2020. starsim/stisim updated to 3.7.2/1.7.0. |
| 28 Oct 2026 (exp 027, steps 1–2b) | Daniel Citron | PrEP added back, from programme data (oral and LA). Testing split by sex and fitted to SHIMS2. SHIMS2-based inputs moved to early 2017 (mid-fieldwork). |
| 18 Nov 2026 (step 2c) | Daniel Citron | ANC testing schedule (1992/2003/2010), PMTCT (2003→2010) and logistic viral suppression added; separate testing start year for men and women; Kasturi's series used as a guide on the men:women balance of diagnosis. |
| 19 Nov 2026 (step 2d) | Daniel Citron | Testing and forced ART split by age cohort (0–24, 25–34, 35–49, 50+; children as 15–24); ANC made logistic, starting 1992.5; VMMC corrected to SHIMS2 totals; testing recalibrated by age and sex. **This snapshot** (commit `068dd83`). Adam's 7 transmission parameters are unchanged throughout. |

## What's in the model (exp 027 additions)

**Daniel:** marks a decision made by Daniel Citron. *Implementation:* marks a
choice made while coding it, open to revision.

| Component | Setting | Data provenance | Notes | Code |
|---|---|---|---|---|
| General-population HIV testing | separate start year by sex; level by sex × age cohort (0–24, 25–34, 35–49, 50+; children share 0–24), fitted to SHIMS2 % diagnosed | **Calibrated, not observed.** Fitted to SHIMS2 2016–17 (Final Report Table 10.3.B, via `SHIMS_909090_coverage.csv`): % diagnosed and % on ART among diagnosed, by sex, for 15–24, 25–34, 35–49 and 15+. Ramp shape from Adam's original testing ramp (2020 and 2040 levels). | **Daniel:** testing rates are fixed inputs fitted to SHIMS, not part of Adam's calibration; split by these age cohorts; under-15s take 15–24 rates (no data for children); SHIMS2 compared at early 2017; SHIMS3 (2021) held out; Kasturi's modelled series is a guide on the men:women balance, not a target. *Implementation:* 10 parameters fitted with optuna; FSW and CD4<200 testing kept as Adam's ramps. | `interventions.py`, `params/calib2d_chosen.json` |
| ANC testing | logistic: starts 1992.5, midpoint 2003, ~max (90%) by 2010 | Researcher-specified schedule. The 90% maximum is carried over from the earlier EMOD calibration, not a survey. | **Daniel:** start 1992.5, midpoint 2003, maximum by 2010; start low enough to bring early diagnosis down. *Implementation:* logistic curve reaching 95% of the maximum by 2010. Previously ~100% of pregnant women were tested from 1985. | `interventions.anc_test_prob_logistic` |
| PMTCT | prophylaxis for diagnosed pregnant women not on ART, 2003 → 2010 | Researcher-specified programme timing. Efficacy 0.96 is stisim's default for mothers on ART. | **Daniel:** begins 2003, maximum by 2010. *Implementation:* maximum set to 90% (not specified; matches ANC); applies only to women not already on ART; pre-2010 short-course regimens were likely less effective. Affects child infections only, so no fitted adult target. | `pmtct.py` |
| ART | forced coverage among PLHIV by sex × the same cohorts, from SHIMS 90-90-90 (0–14 take 15–24) | SHIMS1 2011 (from Adam's `art_coverage.csv`, remapped to the new cohorts); SHIMS2 2016–17 Table 10.3.B; SHIMS3 2021. Coverage = % diagnosed × % on ART among diagnosed. | **Daniel:** same age cohorts as testing; children take 15–24 values; SHIMS2 placed at 2017. *Implementation:* SHIMS2 has no 50+ band, so 50+ for 2017 was backed out of the 15+ and 15–49 rows using unweighted sample sizes (approximate). ART only reaches diagnosed people. Yearly "initiations" include restarts as the model keeps each band on target. | `data/art_coverage_cohorts.csv` |
| Viral suppression on ART | logistic through SHIMS2 (2017) and SHIMS3 (2021) | SHIMS2 Table 10.3.B (15+) and SHIMS3 Table 9.1.B (15+), by sex. | **Daniel:** suppression should rise gently along a logistic curve matching the 2017→2021 increase. *Implementation:* the 80% floor (early ART era) and 98% ceiling are assumptions; one curve per sex for all ages, although SHIMS2 shows much lower suppression in young men (58%). | `vls_construction.to_vls_logistic`, `vls_stock_target.py` |
| VMMC | forced circumcision coverage by age; SHIMS2 medical + non-medical totals at 2017, SHIMS3 at 2021 | 1990: assumed low baseline; 2007: SDHS 2006–07 (50+ extrapolated); 2017: SHIMS2 Final Report Table 15.5.A; 2021: SHIMS3 Table 12.5. Ages 10–14 have no survey data (extrapolated). | **Daniel:** replace the old 2016 values with SHIMS2 medical + non-medical totals (40–49 and 60–64 had been 5–8 points too low) and place them at 2017. Efficacy 0.6 is stisim's default from the trials. Coverage is joined by straight lines between survey years, so yearly circumcisions jump at those years. | `data/vmmc_coverage.csv`, `interventions.SHIMS2_VMMC` |
| PrEP | oral TDF (PrEPWatch counts, 75% efficacy) and LA-PrEP (MOH counts, 92%); 1-year courses; annual initiations; FSW, clients, >1 partner, no stable partner | Oral: PrEPWatch annual distribution, 2019–2025. LA-PrEP: Ministry of Health (personal communication), 4,000/year from 2024. Oral efficacy: Baeten et al., NEJM 2012 (doi:10.1056/NEJMoa1108524). Values after 2025 are carried forward. | **Daniel:** counts are initiations per year, spread over the year; 1-year courses; 75% and 92% efficacy; target FSW, clients of FSW, people with multiple partners and unmarried people; HIV-negative only; no second initiation while on PrEP. *Implementation:* "unmarried" = no current stable partner (the model has no marriage); sexually active 15+ only; users picked at random within the eligible group; re-initiation allowed after a course ends. 75% is the trial's TDF-FTC arm (TDF alone: 67%). | `prep_initiations.py`, `data/prep_initiations.csv` |

Full rationale, fits and caveats are in hivsim_eswatini's
`experiments/027_cascade_prep_recal/README.md`.

## Use

```python
from baseline import make_baseline_sim
sim = make_baseline_sim(seed=1)          # 1985-2026; kwargs pass to run_sims.make_sim
sim.run()
```

`python baseline.py` runs seed 1 and prints prevalence, new infections, HIV
deaths and people on ART beside the mean of the 3.7.2 reference
(`reference/baseline_027d_10seeds.df`, v3 result names mapped to v4's). v4 uses
different random streams and records stocks at the start of each step, so the
two agree statistically at best, not run for run.

## Files

- `baseline.py` — the calibrated parameters and `make_baseline_sim()`.
- `run_sims.py`, `interventions.py`, `analyzers.py`, `pmtct.py`,
  `prep_initiations.py`, `vls_construction.py`, `vls_stock_target.py` —
  ported to Starsim v4 (see the status note at the top).
- `art_coverage_cohort_construction.py` — rebuilds `data/art_coverage_cohorts.csv`
  from `data/SHIMS_909090_coverage.csv` and `data/art_coverage.csv`
  (paths relative to this folder; verified to reproduce the shipped file exactly).
- `params/calib2d_chosen.json` — the fitted testing inputs.
- `reference/` — 10-seed output of the 3.7.2 configuration.

### `data/`

Model inputs (read on every run):

| File | What it is |
|---|---|
| `art_coverage_cohorts.csv` | Forced ART coverage by sex × age cohort (built by `art_coverage_cohort_construction.py`) |
| `prep_initiations.csv` | Oral and LA-PrEP initiations per year; byte-identical to hivsim_eswatini's `DTC_check/SHIMS_Datasets/SWZ_PrEP_Distribution_Data.csv` |
| `vmmc_coverage.csv` | Circumcision coverage by age (2016 row replaced by SHIMS2 totals in code, `interventions.SHIMS2_VMMC`) |
| `condom_use.csv`, `init_prev_hiv.csv`, `eswatini_hiv_calib.csv` | Network condom use, initial prevalence, UNAIDS series attached to the sim |
| `eswatini_age_1985.csv`, `eswatini_asfr.csv`, `eswatini_deaths.csv`, `eswatini_migration.csv` | Demographics |

Source and reference data (not read by a baseline run), copied unchanged from
hivsim_eswatini's `DTC_check/SHIMS_Datasets/` and `data/`:

| File | What it is | Used for |
|---|---|---|
| `SHIMS_909090_coverage.csv` | SHIMS2 (2017) and SHIMS3 (2021) % diagnosed, on ART of diagnosed, VLS of on ART, by sex and age | Source of the ART input; calibration targets for the testing rates (2017; 2021 held out) |
| `aware_and_art_no_interactions_modeled_output_1995_2021-From_Kasturi.csv` | Kasturi's modelled % aware and % on ART among aware, 1995–2021, by sex and 5-year age | Guide on the men:women balance of diagnosis in the testing calibration; comparison plots |
| `SHIMS_incidence_data.csv` | SHIMS2 and SHIMS3 incidence by sex and age | Comparison plots only |
| `art_coverage.csv` | Adam's original ART input ([15,25) … 45+), used up to exp 027c | 2011 (SHIMS1) rows of the cohort ART input; `vls_construction.derive_given_art_by_age`; `make_sim(art_coverage=None)` |
| `SWZ_art_coverage.csv` | An older (Nov 2023) ART-coverage table | Not used by any code; kept for reference |

The code reads `data/` relative to the working directory; `baseline.py`
changes into this folder on import. Paths in comments and docstrings such as
`DTC_check/...`, `calibration_data/...` and `experiments/...` refer to the
source repository, hivsim_eswatini, and are kept as provenance.
