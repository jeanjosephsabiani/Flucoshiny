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
# Structural model and a posteriori (Bayesian MAP = conditional mode) estimation
# of individual CL and V from measured concentrations. Population parameters are
# held fixed to the final values reported in Table 2.
# Simplified example: regular dosing (same dose every tau hours) and covariates
# constant over the dosing history.
# =============================================================================

library(saemix)

# ---- Final population parameters (Table 2), all fixed -----------------------
CL_pop   <- 0.70    # L/h (CKD-EPI = CKDref, no cirrhosis)
V_pop    <- 63.74   # L   (TBW = TBWref)
beta_CKD <- 0.88    # power, CKD-EPI on CL
beta_cir <- -0.44   # cirrhosis on CL
beta_TBW <- 1.66    # power, TBW on V
CKDref   <- 73.7    # mL/min/1.73 m2
TBWref   <- 73.75   # kg
omega_CL <- 0.23
omega_V  <- 0.30
b_prop   <- 0.127   # proportional residual error

# ---- Structural model --------------------------------------------------------
# psi   : one row per subject, columns (CL, V) at the reference covariates
# id    : subject index of each observation
# xidep : time (h since first dose), dose (mg), tinf (h), tau (h),
#         ckd (mL/min/1.73 m2), tbw (kg), cirr (0/1)
model_fluco <- function(psi, id, xidep) {
  CL <- psi[id, 1] * (xidep$ckd / CKDref)^beta_CKD * exp(beta_cir * xidep$cirr)
  V  <- psi[id, 2] * (xidep$tbw / TBWref)^beta_TBW
  t <- xidep$time; tinf <- xidep$tinf; tau <- xidep$tau
  k <- CL / V; R <- xidep$dose / tinf
  cp <- numeric(length(t))
  for (i in seq_along(t)) {
    for (d in 0:floor(t[i] / tau[i])) {        # superposition of all doses given before t
      x <- t[i] - d * tau[i]
      if (x <= 0) next
      cp[i] <- cp[i] + if (x <= tinf[i]) (R[i] / CL[i]) * (1 - exp(-k[i] * x))
                       else (R[i] / CL[i]) * (1 - exp(-k[i] * tinf[i])) * exp(-k[i] * (x - tinf[i]))
    }
  }
  cp
}

saemix_model <- saemixModel(
  model = model_fluco, description = "Fluconazole 1-cmt IV infusion",
  psi0 = matrix(c(CL_pop, V_pop), ncol = 2, dimnames = list(NULL, c("CL", "V"))),
  transform.par = c(1, 1),                      # log-normal CL and V
  covariance.model = diag(2),                   # no CL-V correlation
  omega.init = diag(c(omega_CL^2, omega_V^2)),
  error.model = "proportional", error.init = c(0, b_prop))

# ---- A posteriori estimation (MAP) -------------------------------------------
# obs_df: one row per measured concentration, columns
#   id, time, dose, tinf, tau, ckd, tbw, cirr, conc (mg/L)
estimate_MAP <- function(obs_df) {
  dat <- saemixData(name.data = obs_df, name.group = "id",
                    name.predictors = c("time", "dose", "tinf", "tau", "ckd", "tbw", "cirr"),
                    name.response = "conc", verbose = FALSE)
  # saemix needs a fitted object: a short run is made, then the population
  # parameters are reset to the fixed final values before the MAP step.
  fit <- saemix(saemix_model, dat, list(nbiter.saemix = c(5, 5), map = FALSE, fim = FALSE,
                ll.is = FALSE, displayProgress = FALSE, print = FALSE, save = FALSE,
                save.graphs = FALSE, seed = 632545))
  fit@results@fixed.psi     <- c(CL_pop, V_pop)
  fit@results@fixed.effects <- c(CL_pop, V_pop)
  fit@results@betas         <- t(c(log(CL_pop), log(V_pop)))
  fit@results@omega         <- diag(c(omega_CL^2, omega_V^2))
  fit@results@respar        <- c(0, b_prop)
  fit <- saemix::estimateMeanParametersNewdata(fit)
  fit <- saemix::estimateIndividualParametersNewdata(fit, type = "mode")
  psi <- fit@results@map.psi                     # CL, V at the reference covariates
  cov <- obs_df[!duplicated(obs_df$id), ]
  data.frame(id   = cov$id,
             CL_i = psi[, 1] * (cov$ckd / CKDref)^beta_CKD * exp(beta_cir * cov$cirr),
             V_i  = psi[, 2] * (cov$tbw / TBWref)^beta_TBW)
}
# Individual steady-state AUC0-24 for a daily dose D: D / CL_i

# ---- Example -----------------------------------------------------------------
# obs <- data.frame(id = c(1, 1), time = c(23, 71), dose = 400, tinf = 0.5, tau = 24,
#                   ckd = 30, tbw = 80, cirr = 1, conc = c(6.1, 14.2))
# estimate_MAP(obs)
