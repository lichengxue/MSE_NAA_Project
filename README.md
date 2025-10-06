# Random effects on numbers-at-age transitions implicitly account for movement dynamics and improve stock assessment and management

This repository contains the code, configuration files, and supporting materials for the manuscript:

> Li, C., Deroba, J. J., Berger, A. M., Goethel, D. R., Langseth, B. J., Schueller, A. M., & Miller, T. J.  
> *Random effects on numbers-at-age transitions implicitly account for movement dynamics and improve stock assessment and management.*  
> *Canadian Journal of Fisheries and Aquatic Sciences* (in press). 

---

## 📦 Purpose

This repository accompanies the accepted CJFAS manuscript and provides the complete set of model code, configurations, and outputs used to evaluate how random effects on numbers-at-age transitions can implicitly represent movement dynamics in spatial stock assessment and management simulations.

The analyses were conducted using the **Woods Hole Assessment Model (WHAM)** and an **MSE (Management Strategy Evaluation)** framework implemented in R and TMB.

---

## 📂 Repository structure

| Folder | Description |
|---|---|
| `Code-CJFAS/` | Main R and TMB scripts used to generate simulations, estimation models, and manuscript figures |
| `Final/` | Final accepted version of the manuscript (PDF / Rmd) and associated figures and tables |
| `Revision/` | Revised manuscript files and response-to-reviewers materials |
| `wham-CJFAS/wham/` | Archived version of the WHAM package used to run simulations for this study (backup only) |
| `whamMSE-CJFAS/whamMSE/` | Archived version of the whamMSE package used to run simulations for this study (backup only) |
| `README.md` | This file describing repository content and structure |

---

## ⚙️ Reproducibility

Analyses were conducted in **R (≥4.3)** using **Template Model Builder (TMB)**.  
To reproduce results:

1. Clone this repository:
   ```bash
   git clone https://github.com/lichengxue/MSE_NAA_Project.git
   cd MSE_NAA_Project
