# Sprint plan

*Weekly sprints running Tuesday to Monday, on the same cycle as `hivsim_paper`. The phase lengths follow Daniel's proposal.*

## Phases

| Phase | Length | Goal |
|---|---|---|
| 1. Costing integration | ~1 month | Costing and DALY analyzer in stisim; Eswatini cost inputs; a costed status-quo baseline |
| 2. Scenarios | ~2 months | Design, run and analyse the intervention × shock grid, with uncertainty |
| 3. Manuscript | ~1 month | Draft, co-author round, submission |

## Sprints to 15 Dec

| Sprint | Dates | Work | Deliverable |
|---|---|---|---|
| **2** | 27 Oct–2 Nov | Write a scoping doc covering: interventions (VMMC, condoms, oral PrEP, long-acting PrEP, testing, ART); shocks (ART funding loss, commodity supply disruption, testing disruption, each with timing, depth and duration); outcomes (costs, DALYs, ICERs, ranking by time horizon); perspective; discount rate; horizon | `docs/scoping.md`, agreed with Robyn |
| **3** | 3–9 Nov | Reproduce Adam's calibrated Eswatini runs from this repo and pin versions. Design the costing and DALY analyzer (as a stisim issue) | Baseline reproduced; analyzer design issue |
| **4** | 10–16 Nov | Build the analyzer on a stisim branch off the current release candidate: cost per intervention delivered (using `stisim/logistics` unit costs where they fit), DALYs from GBD disability weights, and discounting | Analyzer PR open, with tests |
| **5** | 17–23 Nov | Assemble Eswatini cost inputs with their sources; merge the analyzer PR; run the costed baseline and sanity-check it against published Eswatini HIV spending | `results/cost_inputs.csv`; baseline costs and DALYs; **Phase 1 done** |
| **6** | 24–30 Nov | Implement the shock mechanics: ART funding loss (as a coverage or capacity cut) and supply disruption (via `SuppliedIntervention`) | Shock scenarios run end to end |
| **7** | 1–7 Dec | Implement the intervention scale-up scenarios; run a pilot grid locally (all interventions × one shock, small number of draws) | First ICER table and ranking |
| **8** | 8–15 Dec | Review the pilot; finalise the scenario grid; plan the full runs and uncertainty analysis | **Interim deliverable (15 Dec):** costed baseline, pilot ranking under one shock, plan for the rest of Phase 2 |

## After 15 Dec (to be refined at Sprint 8)

- **Phase 2, continued:**
  - full grid on a VM;
  - probabilistic sensitivity analysis;
  - ranking stability across time horizons and discount rates;
  - the headline question: when do one-off interventions overtake ongoing ones?
- **Phase 3:** manuscript, co-author round, submission.

## Decisions (to settle in Sprint 2)

1. **Starsim version:** v3, matching the existing calibration, or v4 alongside `hivsim_paper`.
2. **Shock definitions:** stylised (e.g. a 50% ART funding cut from 2026) or tied to real scenarios (e.g. the PEPFAR transition).
3. **Target journal.**
4. **Use of `starsim_econ`:** probably not, since `stisim/logistics` already covers the costs needed.
