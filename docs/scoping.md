# Evaluating the cost-effectiveness of HIV interventions in the face of exogeneous shocks and disruptions

## Background

* Eswatini is a country in southern Africa with the highest HIV prevalence in the world
    * Population: 1.27M
    * PLHIV: 220,000
* Eswatini’s HIV control program has achieved tremendous success
* First country to reach 95-95-95 coverage with ART
* NIH grant (R01HL169021-03) led by Anna Bershteyn (NYU) and Paul Revill (U of York) to assist with planning a health benefits package (HBP) for determining which healthcare-related investments should be prioritized in the future

## Project Goals

* Develop Health Benefits Package that is robust to exogeneous shocks
* Cost-effectiveness evaluations often assume that the status quo will remain stable for many years into the future
  * 2025 saw a massive disruption to HIV services globally
  * Robust HBP planning will account for exogenous shocks to the healthcare system
* Re-evaluate cost-effectiveness estimates of HIV interventions if we assume that disruptions to HIV services are likely to occur in the future
* Pilot project for developing HIVSim and establishing it as a credible HIV modeling tool

## Hypothesis

The value of VMMC will be the most robust in the presence of exogeneous shocks and investing in VMMC today will serve as an “insurance policy” against exogenous shocks in the future


## Methods

* Extend previous HIV service disruptions work
* How do our cost effectiveness evaluations change if we believe that in the future HIV services and availability of interventions will become disrupted?


* Estimate incremental cost-effectiveness ratio (ICER, dollars per DALY averted) for one person-year of ART, VMMC, Oral PrEP, and long-acting PrEP
* Implement future disruptions to HIV services: Vary duration, magnitude, and frequency
* Obtain distribution of ICER values across distribution of possible HIV service disruptions


* Modeled scenarios: range of possible future disruptions varying in
  * Duration 
  * Magnitude
  * Frequency (or, how soon they happen)
* Modeled outcomes:
  * HIV acquisitions, deaths, DALYs
  * Cost of interventions
  * ICER + uncertainty ranges from range of possible future disruptions, Measured as impact of intervention vs. cost of intervention

### Other software development work

Write out starsim module dedicated to costing analysis, calculating ICERs and other health economic metrics