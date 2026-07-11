# FlucoShiny

**Bayesian dose optimization of intravenous fluconazole for critically ill (digestive ICU) patients.**

FlucoShiny is an R [Shiny](https://shiny.posit.co/) application that performs *a posteriori* (Bayesian) individual pharmacokinetic estimation of fluconazole from measured plasma concentrations, and proposes personalized dosing regimens based on a PK/PD (AUC-based) target. Population PK is simulated with [mrgsolve](https://mrgsolve.org/) and individual (conditional posterior) estimation is carried out with [saemix](https://github.com/saemixdevelopment/saemixextension).

> ⚠️ **Research prototype — not a medical device.** This tool is intended to *support*, never to *replace*, the judgment of a qualified physician or clinical pharmacologist. Every dose it produces must be reviewed and validated by a competent healthcare professional in the full clinical context of the patient. See the [Disclaimer](#disclaimer) below.

---

## Features

- **A priori simulation** — population PK predictions from patient covariates (weight, CKD-EPI renal function, cirrhosis).
- **A posteriori Bayesian estimation** — individual CL and V estimated from measured concentrations using the conditional posterior distribution (SAEM / MCMC sampling), following the Levoshiny approach (Lemaitre et al.).
- **AUC-based PK/PD targeting** — computes AUC₀–₂₄, steady-state AUC, and Probability of Target Attainment (PTA) against a user-defined AUC (or AUC/MIC) target.
- **Dose recommendation engine** — scans candidate daily doses and recommends the lowest regimen reaching the target PTA.
- **Concentration–time profiles** — MAP curve, conditional mean, and 95% credible interval plotted against observed samples.
- **Time-varying renal function** — per-observation CKD-EPI values interpolated across the dosing history.
- **PDF report generation** via `rmarkdown` / `knitr`.

## Pharmacokinetic model

A 1-compartment IV-infusion model estimated in a digestive ICU population:

| Parameter | Value | Description |
|-----------|-------|-------------|
| Typical CL | 0.70 L/h | Clearance (reference) |
| Typical V | 63.74 L | Volume of distribution (reference) |
| CKD effect on CL | power, β = 0.88 | CKD-EPI covariate (ref 73.7 mL/min) |
| Cirrhosis effect on CL | β = −0.44 | Binary covariate |
| Weight effect on V | power, β = 1.66 | Total body weight (ref 70.5 kg) |
| IIV | ω(CL)=0.23, ω(V)=0.30 | Inter-individual variability |

Model reference: *Bendjilali-Sabiani et al. — 1-compartment IV fluconazole PK model, digestive ICU population, 2026.*

## Requirements

- **R** ≥ 4.1
- A C++ toolchain (Rtools on Windows, Xcode CLT on macOS) — required by `mrgsolve` to compile the ODE model.
- R packages:

```r
install.packages(c("shiny", "bslib", "rmarkdown", "knitr", "mrgsolve", "saemix"))
```

## Running the app

Clone the repository and launch the app from R:

```r
# from the repo root
shiny::runApp("app.R")
```

or from the command line:

```bash
Rscript -e 'shiny::runApp("app.R", launch.browser = TRUE)'
```

On first launch, `mrgsolve` compiles the model (this can take a moment). A default access password is defined in the source (`AP_PASSWORD`) — change it before any shared deployment.

## Project structure

```
FlucoShiny/
├── app.R          # the complete Shiny application (UI + server + PK engine)
├── README.md
├── LICENSE
└── .gitignore
```

## Disclaimer

This software is provided **"AS IS"**, without warranty of any kind, express or implied. In no event shall the authors, contributors, or their affiliated institutions be liable for any claim, damages, or other liability arising from the software, its use, misuse, modification, or any clinical or other decision made on the basis of its outputs.

This tool is a research prototype intended to **support — never to replace —** the clinical judgment of a qualified physician or pharmacologist. All dose recommendations must be reviewed, validated, and adjusted by a competent healthcare professional taking into account the full clinical context of the patient (comorbidities, concomitant treatments, organ function, microbiology, and local guidelines). Use of this software implies acceptance of these terms.

*(Version française du disclaimer disponible en tête du fichier `app.R`.)*

## License

Released under the [GNU Affero General Public License v3.0 (AGPL-3.0)](LICENSE).

This means anyone may use, study, and modify the software, but **any modified version — including one made available to users over a network (e.g. a hosted Shiny app) — must also be published under the AGPL-3.0 with its complete source code.** The AGPL governs reuse of the code; it does **not** waive the medical disclaimer above.

## Citation

If you use FlucoShiny in academic work, please cite:

> Bendjilali-Sabiani et al. *FlucoShiny: Bayesian dose optimization of intravenous fluconazole in the digestive ICU.* 2026.
