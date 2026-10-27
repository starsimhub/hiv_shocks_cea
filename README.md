# hiv_shocks_cea

**Cost-effectiveness of HIV interventions in the face of external shocks.**

## Question

How does the cost-effectiveness ranking of HIV interventions change under different future funding and supply scenarios? The hypothesis is that one-off interventions such as voluntary medical male circumcision (VMMC) look relatively better in worlds where funding for ongoing treatment is cut. Interventions that need continued supply, such as condoms, PrEP and ART, have benefits that depend on that funding continuing.

- **Setting:** Eswatini, using the calibrated HIVsim model in [`hivsim_eswatini`](https://github.com/starsimhub/hivsim_eswatini) (Adam Akullian).
- **Purpose:** a pilot project that fulfils one of our grant aims and shows what HIVsim can do as a research tool.
- **Lead:** Daniel Citron.

## Where code lives

| Layer | Repo | Holds |
|---|---|---|
| Disease model | `stisim` | Generic costing and DALY analyzer, usable by any HIVsim model |
| Localization | `hivsim_eswatini` | Calibrated Eswatini model and its parameters |
| Analysis | this repo | Shock and intervention scenarios, cost inputs, runs, figures, manuscript |

## Layout

- `docs/plan.md`: sprint plan
- `results/`: committed; everything the plotting scripts need
- `raw_results/`, `raw_data/`: gitignored

## Setup

Follow the setup in `hivsim_eswatini`. Versions will be pinned here once the baseline is reproduced (Sprint 3).
