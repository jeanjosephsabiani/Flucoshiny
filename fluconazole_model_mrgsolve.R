# =============================================================================
# Fluconazole population PK model - mrgsolve implementation
# One-compartment, first-order elimination, IV infusion
# Covariates: CKD-EPI and cirrhosis on CL; total body weight (TBW) on V
#
# Companion model code for:
#   Bendjilali-Sabiani et al. "Impact of liver cirrhosis on fluconazole
#   pharmacokinetics in critically ill patients: model-informed precision
#   dosing through a validated R-Shiny application." CMPB, 2026.
#
# This file contains ONLY the structural PK model and its parameters
# (as reported in Table 2 of the manuscript). It is sufficient to reproduce
# the a priori Monte-Carlo simulations.
# =============================================================================

library(mrgsolve)

code_fluconazole <- '
$PARAM @annotated
TVCL       : 0.70  : Typical clearance (L/h)
TVV        : 63.74 : Typical volume of distribution (L)
BETA_CKD   : 0.88  : Power coefficient, CKD-EPI on CL
BETA_CIRRH : -0.44 : Effect of cirrhosis on CL
BETA_WT    : 1.66  : Power coefficient, TBW on V
CKDref     : 73.7  : Reference CKD-EPI (mL/min)
WTref      : 70.5  : Reference total body weight (kg)

$PARAM @annotated @covariates
CKD   : 90 : CKD-EPI (mL/min)
WT    : 70 : Total body weight (kg)
CIRRH :  0 : Cirrhosis (0 = no, 1 = yes)

$CMT @annotated
CENT : Central compartment (mg) [ADM] [OBS]

$OMEGA @annotated @block
ETA_CL : 0.0529 : IIV on CL (variance; omega_CL = 0.23)
ETA_V  : 0.00   0.09 : IIV on V (variance; omega_V = 0.30)

$SIGMA @annotated
PROP : 0.0144 : Proportional residual error (variance; b = 0.12)

$MAIN
double CL = TVCL * pow(CKD/CKDref, BETA_CKD) * exp(BETA_CIRRH*CIRRH) * exp(ETA_CL);
double V  = TVV  * pow(WT/WTref,   BETA_WT)  * exp(ETA_V);

$ODE
dxdt_CENT = -(CL/V) * CENT;

$TABLE
double CP = CENT / V;
double DV = CP * (1 + PROP);

$CAPTURE CP DV CL V
'

mod_fluco <- mcode("fluconazole", code_fluconazole)

# -----------------------------------------------------------------------------
# Minimal example: steady-state AUC0-24 for a 400 mg q24h, 30-min infusion
# in a 70 kg, CKD-EPI 90, non-cirrhotic patient (n = 1000 virtual subjects).
# -----------------------------------------------------------------------------
# ev <- ev(amt = 400, rate = 400/0.5, ii = 24, addl = 13, cmt = "CENT")
# out <- mod_fluco %>%
#   param(CKD = 90, WT = 70, CIRRH = 0) %>%
#   mrgsim(events = ev, nid = 1000, end = 24*14, delta = 0.1)
