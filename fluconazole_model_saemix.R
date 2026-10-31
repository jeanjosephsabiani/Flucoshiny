# =============================================================================
# Fluconazole population PK model - saemix implementation (a posteriori / MAP)
# One-compartment, first-order elimination, IV infusion
# CKD-EPI and cirrhosis on CL; total body weight (TBW) on V
#
# Companion model code for:
#   Bendjilali-Sabiani et al. "Impact of liver cirrhosis on fluconazole
#   pharmacokinetics in critically ill patients: model-informed precision
#   dosing through a validated R-Shiny application." CMPB, 2026.
#
# This file contains ONLY the structural model and the a posteriori (Bayesian
# MAP / conditional-mode) estimation setup. Population parameters are held
# fixed to the final values reported in Table 2 of the manuscript; saemix is
# used to estimate individual parameters from measured concentrations.
# =============================================================================

library(saemix)

# ---- Final population parameters (Table 2) ---------------------------------
CL_pop   <- 0.70    # L/h
V_pop    <- 63.74   # L
beta_CKD <- 0.88    # power, CKD-EPI on CL
beta_cir <- -0.44   # cirrhosis on CL
beta_TBW <- 1.66    # power, TBW on V
CKDref   <- 73.7    # mL/min
TBWref   <- 70.5    # kg
omega_CL <- 0.23
omega_V  <- 0.30
b_prop   <- 0.12    # proportional residual error

# ---- Covariate-adjusted typical parameters ---------------------------------
typ_CL <- function(ckd, cirr) CL_pop * (ckd/CKDref)^beta_CKD * exp(beta_cir*cirr)
typ_V  <- function(tbw)       V_pop  * (tbw/TBWref)^beta_TBW

# ---- Structural model: 1-cmt IV, 30-min infusion, q24h superposition -------
# psi    : matrix with columns (CL, V), one row per observation (saemix-expanded)
# id     : subject index per observation
# xidep  : data.frame with columns
#            time  : time since first dose (h)
#            dose  : amount per administration (mg)
#            tinf  : infusion duration (h)
#            tau   : dosing interval (h)
model_fluco <- function(psi, id, xidep) {
  CL <- psi[, 1]; V <- psi[, 2]
  t  <- xidep$time; dose <- xidep$dose; tinf <- xidep$tinf; tau <- xidep$tau
  k  <- CL / V; R <- dose / tinf
  cp <- numeric(length(t))
  for (i in seq_along(t)) {
    nd <- floor(t[i] / tau[i]) + 1L          # doses started on or before t[i]
    acc <- 0
    for (d in 0:(nd - 1L)) {
      x <- t[i] - d * tau[i]
      if (x <= 0) next
      if (x <= tinf[i]) acc <- acc + (R[i]/CL[i]) * (1 - exp(-k[i]*x))
      else acc <- acc + (R[i]/CL[i]) * (1 - exp(-k[i]*tinf[i])) * exp(-k[i]*(x - tinf[i]))
    }
    cp[i] <- acc
  }
  cp
}

# ---- saemix model object (population parameters fixed) ---------------------
# psi0 holds the TYPICAL CL and V; per-subject covariate adjustment is applied
# to psi0 before fitting (see the per-patient wrapper below), and the fixed
# effects are not re-estimated (a posteriori use).
saemix_model <- saemixModel(
  model         = model_fluco,
  description   = "Fluconazole 1-cmt IV, a posteriori MAP",
  psi0          = matrix(c(CL_pop, V_pop), ncol = 2,
                         dimnames = list(NULL, c("CL", "V"))),
  transform.par = c(1, 1),                       # log-normal CL, V
  covariance.model = diag(2),                    # diagonal Omega
  omega.init    = diag(c(omega_CL^2, omega_V^2)),
  error.model   = "proportional",
  error.init    = c(0, b_prop)
)

# ---- A posteriori individual estimation ------------------------------------
# For each patient, set psi0 to the covariate-adjusted typical CL and V, then
# run saemix with the population kept fixed (nbiter.saemix = c(0, 0)) so that
# only individual parameters (conditional mode = MAP) are estimated from the
# measured concentrations.
#
# saemix_data <- saemixData(obs_df,
#   name.group = "id", name.predictors = c("time","dose","tinf","tau"),
#   name.response = "conc")
# opt <- saemixControl(nbiter.saemix = c(0, 0), nb.chains = 1,
#                      map = TRUE, fim = FALSE, ll.is = FALSE,
#                      displayProgress = FALSE, save = FALSE, save.graphs = FALSE)
# fit <- saemix(saemix_model, saemix_data, opt)
# individual MAP estimates: fit@results@map.psi    (CL, V per subject)
# individual AUC0-24 at steady state: dose / CL_i
#
# Note: CKD-EPI enters as a time-varying regressor; when CKD changes over the
# dosing history, CL is piecewise-constant (last value carried forward) and
# the structural model is integrated across CKD breakpoints (see the Monolix
# model file and the manuscript, section 2.4).
