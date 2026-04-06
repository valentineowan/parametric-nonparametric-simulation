# Parametric vs Non-Parametric Simulation Study

## Overview
This repository contains R code and outputs for a Monte Carlo simulation study examining whether violations of normality affect statistical conclusions when using parametric and non-parametric tests.

## Methods
The study compares:
- t-test vs Mann–Whitney U test
- ANOVA vs Kruskal–Wallis test

across:
- Sample sizes: 5 to 10,000
- Effect sizes: 0, 0.2, 0.5
- Distributions: normal, skewed, heavy-tailed, uniform

## Outputs
- `SimulationStudy.R` — main script
- `FINAL_simulation_results.csv` — simulation results
- Figures:
  - Agreement plots
  - Power differences
  - Heatmap

## Key Findings
- High agreement across most conditions
- Minimal impact of non-normality alone
- Heavy-tailed distributions lead to meaningful differences
- Non-parametric tests perform better under heavy tails

## Reproducibility
Run `SimulationStudy.R` to reproduce all results.

## Requirements
- R version 4.5.2
- Packages: dplyr, ggplot2, officer, flextable

## Author
Valentine Joseph Owan

## License
For academic use