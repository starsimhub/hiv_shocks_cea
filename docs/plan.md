# Sprint plan

*Weekly sprints, Tuesday to Monday, on the same cycle as `hivsim_paper`. Phase lengths follow Daniel's proposal. Runs on Starsim v4.*

## Phases

| Phase | Length | Goal |
|---|---|---|
| 1. Costing integration | ~1 month | Cost-effectiveness framework built into HIVsim; Eswatini cost inputs; a costed status-quo baseline |
| 2. Scenarios | ~2 months | Design, run and analyse the intervention × shock grid, with uncertainty |
| 3. Manuscript | ~1 month | Draft, co-author round, submission |

## Sprints to 15 Dec

| Sprint | Dates | Work | Deliverable |
|---|---|---|---|
| **2** | 27 Oct–2 Nov | Scoping doc: interventions, shocks (timing, depth and duration), outcomes, and the cost-effectiveness framework (Daniel's choice) | `docs/scoping.md`, agreed with Robyn |
| **3** | 3–9 Nov | Get Adam's calibrated Eswatini model running on Starsim v4, alongside the HIVsim v4 port; pin versions. Design the costing analyzer (as a stisim issue) | Eswatini baseline on v4; analyzer design issue |
| **4** | 10–16 Nov | Build the costing analyzer on a stisim branch | Analyzer PR open, with tests |
| **5** | 17–23 Nov | Assemble Eswatini cost inputs with sources; merge the analyzer; run and sanity-check the costed baseline | `results/cost_inputs.csv`; baseline costs and outcomes; **Phase 1 done** |
| **6** | 24–30 Nov | Implement the shock mechanics (e.g. ART funding loss, commodity supply disruption) | Shock scenarios run end to end |
| **7** | 1–7 Dec | Implement the intervention scenarios; run a pilot grid locally (all interventions × one shock, few draws) | First ICER table and ranking |
| **8** | 8–15 Dec | Review the pilot; finalise the scenario grid; plan the full runs and uncertainty analysis | **Interim deliverable (15 Dec):** costed baseline, pilot ranking under one shock, plan for the rest of Phase 2 |

## After 15 Dec (to be refined at Sprint 8)

- **Phase 2, continued:** full grid on a VM; uncertainty analysis; how stable the rankings are across time horizons; when one-off interventions overtake ongoing ones.
- **Phase 3:** manuscript, co-author round, submission.

## Decided

- **Starsim v4.**
- **No `starsim_econ`.**
- **No `stisim/logistics`** (likely to be removed).
- **The cost-effectiveness framework is Daniel's call.**

## To settle in Sprint 2

1. Shock definitions: stylised, or tied to real scenarios (e.g. the PEPFAR transition).
2. Target journal.
