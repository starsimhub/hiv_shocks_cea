# Sprint plan

*Weekly sprints, Tuesday to Monday, on the same cycle as `hivsim_paper`. Phase lengths follow Daniel's proposal. Runs on Starsim v4.*

## Phases

| Phase | Length | Goal |
|---|---|---|
| 1. Compete preparing HIVSim baseline model | October-November 2026 | Cost-effectiveness framework built into HIVsim; Eswatini cost inputs; a costed status-quo baseline |
| 2. Develop scenarios and prototype analysis | November 2026 | Baseline; interventions; shocks; post-processing and analysis |
| 3. Execute analysis and generate figures | December 2026 | Run planned simulations and prepare results |
| 4. Cost-effectiveness framework built into Starsim | December 2026-January 2027 | Add functions for incorporating cost and calculating DALYS to calculate ICERs, CEACs, and other standard health economics analysis |
| 5. Prepare Manuscript | January-February 2027 | Draft, co-author round, submission |

## Sprints to 15 Dec

| Sprint | Dates | Work | Deliverable |
|---|---|---|---|
| **2** | 27 Oct–2 Nov | Scoping doc: interventions, shocks (timing, depth and duration), outcomes, and the cost-effectiveness framework (Daniel's choice) | `docs/scoping.md`, agreed with Robyn |
| **3** | 3–16 Nov | Update Adam's calibrated Eswatini model to include 90-90-90 art coverage metrics as calibration targets, PrEP, ANC | Convert Eswatini model to run on Starsim v4, alongside the HIVsim v4 port; pin versions. Design the costing analyzer (as a stisim issue). Eswatini baseline on v4; analyzer design issue |
| **4** | 17–30 Nov | Finalize baseline SWZ model | Push baseline model to this repo |
| **5** | 1–7 Dec | Implement shock mechanics, prototype Disruptions scenarios and analysis, Assemble Eswatini cost inputs with sources | `results/cost_inputs.csv`; first ICER table and ranking; preliminary results and review with NYU team |
| **6** | 8–15 Dec | Prepare and execute revised scenario | Revised results with ICER tables, plots, rankings **Phase 1 done** |

## After 15 Dec (to be refined at Sprint 8)

- **Phase 2:** full grid on a VM; uncertainty analysis; how stable the rankings are across time horizons; when one-off interventions overtake ongoing ones.
- **Phase 2, continued:** Build the costing analyzer on a stisim branch
- **Phase 3:** manuscript, co-author round, submission.

## Undecided

- How will we "formally" calibrate the Eswatini model? If we decide to use history matching, add to sprint 5 in phase 1 above.

## Decided

- **Starsim v4.**
- **No `starsim_econ`.**
- **No `stisim/logistics`** (likely to be removed).
- **The cost-effectiveness framework is Daniel's call.**

## To settle in Sprint 2

1. Shock definitions: stylised, or tied to real scenarios (e.g. the PEPFAR transition).
2. Target journal (Medical Decision Making)
