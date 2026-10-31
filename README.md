# Fluconazole population PK model — code

Model code accompanying:

> Bendjilali-Sabiani et al. *Impact of liver cirrhosis on fluconazole
> pharmacokinetics in critically ill patients: model-informed precision dosing
> through a validated R-Shiny application.* Computer Methods and Programs in
> Biomedicine, 2026.

This repository contains the **population pharmacokinetic model** used in the
study, in the three environments used in the paper (Monolix, mrgsolve, saemix).
It is sufficient to reproduce the model, the a priori Monte-Carlo simulations,
and the a posteriori Bayesian individual estimation. The R-Shiny application
itself (user interface, server logic, deployment) is not part of this code
release; the application is freely usable online at
<https://flucoshinyicu.shinyapps.io/FlucoShiny/>.

## Model

One-compartment model with first-order elimination and IV (30-min) infusion.

```
CL = 0.70 · (CKD-EPI / 73.7)^0.88 · exp(-0.44 · cirrhosis)      [L/h]
V  = 63.74 · (TBW / 70.5)^1.66                                   [L]
```

Between-subject variability (log-normal): omega_CL = 0.23, omega_V = 0.30
(no CL–V correlation). Residual error: proportional, b = 0.12.
CKD-EPI is a time-varying regressor (last value carried forward). Full
parameter estimates with uncertainty are in Table 2 of the manuscript.

## Files

| File | Purpose |
|------|---------|
| `fluconazole_model_monolix.txt` | Mlxtran structural model (parameter estimation) |
| `fluconazole_model_mrgsolve.R`  | mrgsolve model (a priori Monte-Carlo simulations) |
| `fluconazole_model_saemix.R`    | saemix structural model + a posteriori MAP setup |

## Requirements

R (>= 4.2) with `mrgsolve` and `saemix`; Monolix (>= 2024R1) for the Mlxtran
model.

## Target

PK/PD target: AUC0-24h / MIC >= 100. Worst-case clinical breakpoint MIC = 4 mg/L.
