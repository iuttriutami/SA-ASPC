# SA-ASPC: calibrating the spatial penalty in fuzzy clustering of multivariate time trajectories

Code, data and results for

> Utami IT, Daniaty D. *When does a spatial penalty pay? Calibrating the contiguity coefficient in
> fuzzy clustering of multivariate time trajectories.* Submitted to *Advances in Data Analysis and Classification*.

SA-ASPC selects the spatial penalty β of DTW fuzzy c-medoids clustering with a weighted ideal-point rule over
four objectives (fuzzy Moran index, fuzzy silhouette, Xie–Beni, perturbation stability) and certifies the
selection with an exact Hoeffding bound and a Monte Carlo confidence set.

## Repository layout

| Folder | Content |
|---|---|
| `R/` | `SA_ASPC_pipeline_v2.R`: core functions (DTW, safeguarded fuzzy c-medoids with the corrected penalty weight, objectives, bank, selection, certificate) |
| `scripts/application1/` | Central Java, 35 regencies/cities, 2022–2025 (Sections 5.1–5.6, Table 5) |
| `scripts/application2/` | Java, 114 regencies/cities, 2020–2024 (Section 5.7), contiguity matrix and deep choice of C |
| `scripts/simulation/` | Sections 6.1–6.5 and Appendix C |
| `scripts/validation/` | Certificate coverage (Table S20), near-tie validation (Table S21), exact optimum at β = 0 (Table S22), change-point design (Table S23) |
| `app/` | `app_SA_ASPC.R`: Shiny explorer of the selection, the certificate and the weight regions |
| `data/` | Assembled panels and frozen contiguity matrices |
| `results/` | CSV files behind the tables of the paper and of Online Resource 1 |

## Reproducing the results

All scripts are written for R ≥ 4.3 and were run on Kaggle notebooks (4 CPU cores). Each script reads its
inputs from `/kaggle/input` and writes checkpoints and outputs to `/kaggle/working`; to run locally, change
the `WORK` and `inp()` definitions at the top of the script to point to this repository.
Long runs save a checkpoint after every block, so an interrupted session can be resumed by attaching the
previous output as input.

Required packages: those loaded by `R/SA_ASPC_pipeline_v2.R`, plus `readxl`, `Rcpp` (exact optimum), `sf` (maps and contiguity), `shiny` and `ggplot2` (app).

| Paper | Script | Output |
|---|---|---|
| Sections 5.1–5.5, Tables 1–4 | `scripts/application1/R_aplikasi_pertama_C.R` (`C_FIXED <- 8`) | `app1_out_C8/` |
| Section 5.6, Table 5 | same script with `C_FIXED <- 4` and `8` | `results/application1/table5_C4_vs_C8.csv` |
| Section 5.7, Tables S4–S17 | `scripts/application2/R_buat_W_jawa.R`, then `R_aplikasi_kedua.R` | `results/application2/` |
| Choice of C under the deep protocol (Table S19) | `scripts/application2/R_deep_Cselect.R` | `tableS19*.csv` |
| Section 5.4, Table S20 | `scripts/validation/R_coverage_certificate.R` | `results/validation/` |
| Section 5.4, Table S21 | `scripts/validation/R_near_tie_certificate.R` | `results/validation/` |
| Section 4, Table S22 | `scripts/validation/R_exact_optimum.R` | `results/validation/` |
| Section 6.5, Table S23 | `scripts/validation/R_sim_changepoint.R` | `results/validation/` |

## Data sources

Development indicators: Statistics Indonesia (BPS). Administrative boundaries: Geospatial Information Agency
(BIG) for the first application and GADM version 3.6 for the second. See `data/README.md`.

## Citation

See `CITATION.cff`. A versioned archive of this repository is available on Zenodo: [DOI to be added].

## Licence

Code: MIT licence (`LICENSE`). Data assembled from public sources; please cite the original providers.
