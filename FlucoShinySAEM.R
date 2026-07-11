# ==============================================================================
# FlucoShiny v6.0 - mrgsolve + saemix engine (a posteriori)
# ==============================================================================
#
# DISCLAIMER / AVERTISSEMENT
# ---------------------------------------------------------------------------
# This software is provided "AS IS", without warranty of any kind, express or
# implied, including but not limited to the warranties of merchantability,
# fitness for a particular purpose, and noninfringement. In no event shall the
# authors, contributors, or their affiliated institutions be liable for any
# claim, damages, or other liability — whether in an action of contract, tort,
# or otherwise — arising from, out of, or in connection with the software, its
# use, misuse, modification, or any decision (clinical or otherwise) made on
# the basis of its outputs.
#
# This tool is a research prototype intended to support — never to replace —
# the clinical judgment of a qualified physician or pharmacologist. All
# dose recommendations produced by this software must be reviewed,
# validated, and adjusted by a competent healthcare professional taking into
# account the full clinical context of the patient, including comorbidities,
# concomitant treatments, organ function, microbiology, and local guidelines.
# The authors decline all responsibility for any adverse outcome, treatment
# failure, injury, or death that may occur in connection with the use of
# this software.
#
# Use of this software implies acceptance of these terms.
# ---------------------------------------------------------------------------
#
# Ce logiciel est fourni « EN L'ÉTAT », sans garantie d'aucune sorte, expresse
# ou implicite, y compris mais sans s'y limiter les garanties de qualité
# marchande, d'adéquation à un usage particulier et de non-contrefaçon. En
# aucun cas les auteurs, contributeurs ou leurs établissements de
# rattachement ne sauraient être tenus responsables d'une quelconque
# réclamation, dommage ou autre obligation — qu'il s'agisse d'une action
# contractuelle, délictuelle ou autre — résultant de, ou en lien avec ce
# logiciel, son utilisation, sa mauvaise utilisation, sa modification ou
# toute décision (clinique ou autre) prise sur la base de ses résultats.
#
# Cet outil est un prototype de recherche destiné à appuyer — jamais à
# remplacer — le jugement clinique d'un médecin ou d'un pharmacologue
# qualifié. Toute recommandation de dose produite par ce logiciel doit être
# revue, validée et adaptée par un professionnel de santé compétent prenant
# en compte le contexte clinique global du patient (comorbidités, traitements
# concomitants, fonctions d'organe, microbiologie, recommandations locales).
# Les auteurs déclinent toute responsabilité en cas d'événement indésirable,
# d'échec thérapeutique, de blessure ou de décès survenant en lien avec
# l'usage de ce logiciel.
#
# L'utilisation de ce logiciel implique l'acceptation de ces conditions.
# ==============================================================================

library(shiny)
library(bslib)
library(rmarkdown)
library(knitr)
library(mrgsolve)
library(saemix)

# ==============================================================================
# MRGSOLVE MODEL - Fluconazole 1-compartment IV infusion (ODE)
# ==============================================================================
code_fluconazole <- '
$PARAM @annotated
TVCL : 0.70   : Typical CL (L/h)
TVV  : 63.74  : Typical V (L)
BETA_CKD   : 0.88  : Power covariate CKD on CL
BETA_CIRRH : -0.44 : Effect cirrhosis on CL
BETA_WT    : 1.66  : Power covariate WT on V
CKDref : 73.7 : Reference CKD-EPI (mL/min)
WTref  : 70.5 : Reference weight (kg)
ETA1 : 0 : IIV on CL
ETA2 : 0 : IIV on V

$PARAM @annotated @covariates
CKD   : 90 : CKD-EPI (mL/min)
WT    : 70 : Total body weight (kg)
CIRRH :  0 : Cirrhosis (0/1)

$CMT @annotated
CENT : Central compartment (mg) [ADM] [OBS]

$OMEGA @block
0.0484
0 0.0841

$SIGMA
0.04
0

$MAIN
double CKDc = CKD; if(CKDc < 15) CKDc = 15; if(CKDc > 175) CKDc = 175;
double WTc  = WT;  if(WTc  < 42) WTc = 42; if(WTc  > 123) WTc  = 123;
double CL = TVCL * pow(CKDc / CKDref, BETA_CKD) * exp(BETA_CIRRH * CIRRH) * exp(ETA1 + ETA(1));
double V  = TVV  * pow(WTc  / WTref,  BETA_WT)  * exp(ETA2 + ETA(2));
if(CL < 0.001) CL = 0.001;
if(V  < 0.1)   V  = 0.1;

$ODE
dxdt_CENT = -(CL / V) * CENT;

$TABLE
double CP = CENT / V;
double DV = CP * (1 + EPS(1)) + EPS(2);

$CAPTURE CP DV CL V
'

options(shiny.sanitize.errors = FALSE)
mod_fluco <- mcode("fluconazole_v6safe", code_fluconazole, compile = TRUE, soloc = tempdir())

# ==============================================================================
# CONSTANTS
# ==============================================================================
CLpop <- 0.70; Vpop <- 63.74; beta_ckd <- 0.88; beta_cirrh <- -0.44; beta_wt <- 1.66
CKDref <- 73.7; WTref <- 70.5; omega_CL <- 0.23; omega_V <- 0.30; rho_CL_V <- 0.0  # No CL/V correlation in original Monolix model
NSIM <- 10000L; PTA_OPTIMAL <- 0.90; PTA_ACCEPTABLE <- 0.75
WT_LOWER <- 42.1; WT_UPPER <- 123; CKD_LOWER <- 15; CKD_UPPER <- 175; AGE_ADULT <- 18
AP_PASSWORD <- "fluco2026"
NPOST_SAMP <- 4000L; POST_PROFILE_GRID <- 300L
OMEGA_MAT <- matrix(c(omega_CL^2, rho_CL_V*omega_CL*omega_V, rho_CL_V*omega_CL*omega_V, omega_V^2), nrow=2)
COL_MAP <- "#1565C0"; COL_CM <- "#AD1457"; COL_MEDIAN <- "#37474F"; COL_BAND <- "#90CAF9"
COL_TARGET <- "#B71C1C"; COL_OPT_M <- "#1B5E20"

# ── PK MODEL REFERENCE ──────────────────────────────────────────────────────
# TODO: replace with full citation once published
# Format: Authors. Title. Journal. Year;Vol(Issue):Pages.
PK_MODEL_REF <- "Bendjilali-Sabiani et al. — 1-compartment IV fluconazole PK model, digestive ICU population, 2026."
# ────────────────────────────────────────────────────────────────────────────

# ==============================================================================
# PURE FUNCTIONS
# ==============================================================================
predict_CL_typ <- function(ckd, cirrh) { ckd <- max(5, ckd); CLpop * (ckd/CKDref)^beta_ckd * exp(beta_cirrh*cirrh) }
predict_V_typ <- function(wt) { wt <- max(20, wt); Vpop * (wt/WTref)^beta_wt }

simulate_params <- function(wt, ckd, cirrh) {
  # Clamp to avoid numerical issues
  wt <- max(WT_LOWER, min(wt, WT_UPPER)); ckd <- max(CKD_LOWER, min(ckd, CKD_UPPER))
  set.seed(123)
  sim <- tryCatch(
    mod_fluco %>%
      param(CKD = ckd, WT = wt, CIRRH = cirrh) %>%
      mrgsim(nid = NSIM, end = 0) %>%
      as.data.frame(),
    error = function(e) { message("simulate_params error: ",e$message); NULL })
  CLtyp <- predict_CL_typ(ckd, cirrh); Vtyp <- predict_V_typ(wt)
  if(is.null(sim) || nrow(sim)==0) {
    # Fallback: manual sampling
    eta <- rmvnorm_chol(NSIM, c(0,0), OMEGA_MAT)
    return(list(CL=CLtyp*exp(eta[,1]), V=Vtyp*exp(eta[,2]), CLtyp=CLtyp, Vtyp=Vtyp))
  }
  list(CL = sim$CL, V = sim$V, CLtyp = CLtyp, Vtyp = Vtyp)
}

auc_0_24_inf <- function(dose, CL, V, tinf, tend=24) {
  tinf <- pmax(tinf,1e-6); tend <- pmax(tend,1e-6); k <- CL/V; R <- dose/tinf
  if (tinf >= tend) (R/CL)*(tend-(1-exp(-k*tend))/k) else {
    auc1 <- (R/CL)*(tinf-(1-exp(-k*tinf))/k); auc2 <- (R/CL)*(1-exp(-k*tinf))*(1-exp(-k*(tend-tinf)))/k; auc1+auc2 }
}

# AUC over interval [t_start, t_end] with repeated daily dose (superposition)
# day = number of doses already given before this interval (0=first dose)
auc_loading_day <- function(dose, CL, V, tinf, day_idx) {
  # day_idx: 1=AUC 0-24h (first dose), 2=AUC 24-48h (2 doses given), 3=AUC 48-72h (3 doses given)
  tau <- 24
  dose_times <- seq(0, by=tau, length.out=day_idx)  # doses at 0, 24, 48, ...
  dose_amounts <- rep(dose, day_idx)
  t_start <- (day_idx - 1) * tau
  t_end   <- day_idx * tau
  auc_model_interval(t_start, t_end, dose_times, dose_amounts, CL, V, tinf)
}
auc_ss_24 <- function(dose, CL, V, tinf, tau=24) {
  # AUC at steady-state by analytical integral (14 doses to reach SS)
  n_d <- 14; dt <- seq(0, by=tau, length.out=n_d); da <- rep(dose, n_d)
  auc_model_interval((n_d-1)*tau, n_d*tau, dt, da, CL, V, tinf)
}

summ_auc <- function(mat) {
  if (is.vector(mat)) mat <- matrix(mat, ncol=1); if (!is.matrix(mat)) mat <- as.matrix(mat)
  data.frame(mean=colMeans(mat), median=apply(mat,2,median), p2.5=apply(mat,2,quantile,0.025), p97.5=apply(mat,2,quantile,0.975))
}
summ_quant <- function(x) { q <- quantile(x, probs=c(0.025,0.50,0.975), na.rm=TRUE); c(mean=mean(x,na.rm=TRUE), median=as.numeric(q[2]), p2.5=as.numeric(q[1]), p97.5=as.numeric(q[3])) }

loading_regimens <- function(wt) data.frame(DoseLabel=c("600 mg/day","800 mg/day","1000 mg/day","1200 mg/day","1400 mg/day","1600 mg/day","1800 mg/day","12 mg/kg/day"), Dose=c(600,800,1000,1200,1400,1600,1800,12*wt), stringsAsFactors=FALSE)
maintenance_regimens <- function(wt) data.frame(DoseLabel=c("50 mg/day","100 mg/day","200 mg/day","400 mg/day","600 mg/day","800 mg/day","1000 mg/day","1200 mg/day","6 mg/kg/day","8 mg/kg/day"), Dose=c(50,100,200,400,600,800,1000,1200,6*wt,8*wt), stringsAsFactors=FALSE)

calc_age_years <- function(dob, ref_date=Sys.Date()) { if (is.na(dob)) return(NA_real_); as.numeric(difftime(ref_date, as.Date(dob), units="days"))/365.25 }

ss_conc_1cpt_infusion <- function(time, dose, CL, V, tinf, tau) {
  tinf <- pmax(tinf,1e-8); tau <- pmax(tau,tinf+1e-8); k <- CL/V; R <- dose/tinf
  conc <- ifelse(time<=tinf, (R/CL)*(1-exp(-k*time))/(1-exp(-k*tau)), (R/CL)*(1-exp(-k*tinf))*exp(-k*(time-tinf))/(1-exp(-k*tau)))
  pmax(conc, 1e-12)
}

build_obs_df <- function(n, input) {
  out <- do.call(rbind, lapply(seq_len(n), function(i) data.frame(time=as.numeric(input[[paste0("tconc_",i)]]), conc=as.numeric(input[[paste0("cconc_",i)]]), stringsAsFactors=FALSE)))
  out <- out[is.finite(out$time) & is.finite(out$conc),,drop=FALSE]; rownames(out) <- NULL; out
}

# Build obs from dates: compute time_post_dose from date schedule
build_obs_from_dates <- function(n, input, date_first, heure_first=0, tau=24) {
  d1 <- as.Date(date_first)
  rows <- lapply(seq_len(n), function(i) {
    d <- tryCatch(as.Date(input[[paste0("dconc_",i)]]), error=function(e) NA)
    h <- as.numeric(input[[paste0("hconc_",i)]])
    conc <- as.numeric(input[[paste0("cconc_",i)]])
    ckd_i <- as.numeric(input[[paste0("ckdconc_",i)]])
    if (is.na(d) || is.na(h) || is.na(conc)) return(NULL)
    if (is.na(ckd_i) || !is.finite(ckd_i)) ckd_i <- as.numeric(input$ckd)
    abs_h <- as.numeric(difftime(d, d1, units="hours")) + h - heure_first
    last_dose_h <- floor(abs_h / tau) * tau
    if (last_dose_h >= abs_h) last_dose_h <- last_dose_h - tau  # strictly before obs
    if (last_dose_h < 0) last_dose_h <- 0
    time_post_dose <- abs_h - last_dose_h
    data.frame(date=format(d,"%d/%m/%Y"), heure=h, conc=conc, abs_h=abs_h, time=time_post_dose, ckd=ckd_i, stringsAsFactors=FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out)==0) return(data.frame(date=character(),heure=numeric(),conc=numeric(),abs_h=numeric(),time=numeric(),ckd=numeric(),stringsAsFactors=FALSE))
  out[is.finite(out$conc) & out$conc>0,,drop=FALSE]
}
rmvnorm_chol <- function(n, mu, Sigma) { L <- chol(Sigma); Z <- matrix(rnorm(n*length(mu)), ncol=length(mu)); sweep(Z %*% L, 2, mu, "+") }

# Model-based AUC: analytical integral of 1-cpt IV infusion over [t_start, t_end]
# Sums contributions from all doses (superposition)
auc_model_interval <- function(t_start, t_end, dose_times, dose_amounts, CL, V, tinf) {
  k <- CL/V; auc <- 0
  for (i in seq_along(dose_times)) {
    td <- dose_times[i]; D <- dose_amounts[i]; R <- D/tinf
    ds <- t_start - td; de <- t_end - td
    if (de <= 0) next  # dose not yet given
    ds <- max(ds, 0)   # clamp
    # Infusion phase contribution [max(ds,0), min(de,tinf)]
    if (ds < tinf) {
      t1 <- ds; t2 <- min(de, tinf)
      auc <- auc + (R/CL) * ((t2-t1) - (1/k)*(exp(-k*t1) - exp(-k*t2)))
    }
    # Post-infusion phase contribution [max(ds,tinf), de]
    if (de > tinf) {
      t1 <- max(ds, tinf); t2 <- de
      auc <- auc + (R/CL)*(1-exp(-k*tinf))*(1/k)*(exp(-k*(t1-tinf)) - exp(-k*(t2-tinf)))
    }
  }
  auc
}

# ==============================================================================
# DOSE RECOMMENDATION ENGINE
# ==============================================================================
simulate_dose_candidates <- function(post_CL, post_V, tinf, tau, target_auc, doses) {
  dose_times <- seq(0, by = tau, length.out = 14)
  do.call(rbind, lapply(doses, function(d) {
    dose_amounts <- rep(d, length(dose_times))
    auc_vec <- mapply(function(cl, v) {
      auc_model_interval((length(dose_times)-1)*tau, length(dose_times)*tau,
                         dose_times, dose_amounts, cl, v, tinf)
    }, post_CL, post_V)
    q <- quantile(auc_vec, c(0.025,0.50,0.975), na.rm=TRUE)
    pta <- mean(auc_vec >= target_auc, na.rm=TRUE)
    data.frame(dose=d, median=as.numeric(q[2]), p2.5=as.numeric(q[1]), p97.5=as.numeric(q[3]), pta=pta)
  }))
}
build_dose_recommendation <- function(ap, target_auc, mic) {
  post <- ap$posterior; current_daily <- ap$dose*(24/ap$tau); q_auc_curr <- summ_quant(post$AUC24); pta_current <- mean(post$AUC24>=target_auc, na.rm=TRUE)
  sim <- simulate_dose_candidates(post$CL, post$V, ap$tinf, ap$tau, target_auc, seq(25,2000,by=25))
  opt_idx <- which(sim$pta>=PTA_OPTIMAL); sub_idx <- which(sim$pta>=PTA_ACCEPTABLE)
  list(current_daily=current_daily, auc_current=q_auc_curr, pta_current=pta_current, status=if(q_auc_curr["median"]>=target_auc)"above" else "below",
       target_auc=target_auc, mic=mic, sim=sim, rec_optimal=if(length(opt_idx)>0) sim$dose[min(opt_idx)] else NA_real_, rec_suboptimal=if(length(sub_idx)>0) sim$dose[min(sub_idx)] else NA_real_, candidates=seq(25,2000,by=25))
}

# ==============================================================================
# PLOT HELPERS
# ==============================================================================
plot_fun <- function(sum, pta, target, labels, rec_idx, ylab, main) {
  x <- seq_along(sum$median); cols <- ifelse(pta >= PTA_OPTIMAL, "#2E7D32", "#B0BEC5")
  plot(x, sum$median, pch=19, col=cols, xaxt="n", ylim=range(c(sum$p2.5, sum$p97.5, target)), ylab=ylab, main=main)
  axis(1, at=x, labels=labels); arrows(x, sum$p2.5, x, sum$p97.5, angle=90, code=3, length=0.04)
  abline(h=target, col="red", lty=2, lwd=2)
  if (!is.na(rec_idx)) points(rec_idx, sum$median[rec_idx], pch=21, bg="#1565C0", cex=1.6)
}

plot_ct_profile <- function(prof, obs, tinf, title="Posterior concentration-time profile") {
  ylim_max <- max(c(prof$p97.5, obs$conc, prof$map), na.rm=TRUE)*1.18
  par(mar=c(4.5,4.5,3.5,1.5), bg="white")
  plot(prof$time, prof$map, type="n", xlab="Time since last administration (h)", ylab="Fluconazole concentration (mg/L)", main=title, font.main=2, ylim=c(0,ylim_max), las=1)
  grid(col="#ECEFF1", lty=1)
  polygon(c(prof$time,rev(prof$time)), c(prof$p2.5,rev(prof$p97.5)), col=grDevices::adjustcolor(COL_BAND,0.40), border=NA)
  lines(prof$time, prof$map, lwd=3, col=COL_MAP); lines(prof$time, prof$condmean, lwd=2, col=COL_CM, lty=5)
  abline(v=tinf, lty=3, col="#78909C", lwd=1.5)
  points(obs$time, obs$conc, pch=21, cex=1.6, bg="#FF6F00", col="white", lwd=1.5)
  legend("topright", bty="n", cex=0.82, legend=c("MAP curve","Conditional mean","95% credible CI","Observed concentrations"), lty=c(1,5,NA,NA), lwd=c(3,2,NA,NA), pch=c(NA,NA,15,21), pt.cex=c(NA,NA,2,1.4), pt.bg=c(NA,NA,grDevices::adjustcolor(COL_BAND,.5),"#FF6F00"), col=c(COL_MAP,COL_CM,grDevices::adjustcolor(COL_BAND,.6),"#FF6F00"))
}

plot_dose_recommendation <- function(sim, current_daily, target_auc, rec_optimal, rec_suboptimal) {
  sim_d <- sim[sim$dose>=50 & sim$dose<=1600,]; ymax <- max(c(sim_d$p97.5, target_auc), na.rm=TRUE)*1.15
  par(mar=c(5,5,3,10), bg="white", xpd=FALSE)
  plot(sim_d$dose, sim_d$median, type="n", xlab="Daily dose (mg/day)", ylab="AUC (mg\u00b7h/L)", main="", ylim=c(0,ymax), las=1)
  grid(col="#ECEFF1", lty=1)
  polygon(c(sim_d$dose,rev(sim_d$dose)), c(sim_d$p2.5,rev(sim_d$p97.5)), col=grDevices::adjustcolor("#B3E5FC",0.55), border=NA)
  lines(sim_d$dose, sim_d$median, lwd=3, col=COL_MAP); abline(h=target_auc, col=COL_TARGET, lty=2, lwd=2)
  abline(v=current_daily, col="#546E7A", lty=4, lwd=2)
  if (!is.na(rec_optimal)) abline(v=rec_optimal, col=COL_OPT_M, lty=1, lwd=2.5)
  par(xpd=TRUE)
  legend_e <- c("Post. AUC median","Post. 95% CI","PK/PD target","Current dose")
  legend_lt <- c(1,NA,2,4); legend_lw <- c(3,NA,2,2); legend_pc <- c(NA,15,NA,NA)
  legend_co <- c(COL_MAP,grDevices::adjustcolor("#B3E5FC",.7),COL_TARGET,"#546E7A")
  if (!is.na(rec_optimal)) { legend_e <- c(legend_e,"Recommended"); legend_lt <- c(legend_lt,1); legend_lw <- c(legend_lw,2.5); legend_pc <- c(legend_pc,NA); legend_co <- c(legend_co,COL_OPT_M) }
  legend("right", inset=c(-0.18,0), bty="o", box.col="#E0E0E0", bg="white", cex=0.78, legend=legend_e, lty=legend_lt, lwd=legend_lw, pch=legend_pc, pt.cex=2, col=legend_co)
  par(xpd=FALSE)
}

# ==============================================================================
# SAEMIX ENGINE — individual estimation (a posteriori)
# Remplace l'ancien moteur mapbayr (estimation MAP via mrgsolve).
# Approche IDENTIQUE a Levoshiny (Lemaitre et al.):
#   1) saemixPredictNewdata avec nsamp=1000 -> distribution conditionnelle MCMC
#   2) AUC24 = Dose/CL calcule sur chacun des echantillons
#   3) IC95% = quantiles [2.5% ; 97.5%] de la distribution
# Aucune modification des cibles PK/PD, du modele de population, ou de l'UI
# ==============================================================================

# ---- Helper Levoshiny: saemixPredictNewdata --------------------------------
# Issu integralement de Levoshiny (Lemaitre et al.). Sera implemente dans
# une future version officielle du package saemix. Permet la prediction
# bayesienne (icpred + parSample) pour de nouveaux sujets.
saemixPredictNewdata <- function(saemixObject, newdata,
                                 type = c("ipred","ypred","ppred","icpred"),
                                 nsamp = 1) {
  set.seed(12345)
  saemixObject <- saemix::replaceData.saemixObject(saemixObject, newdata)
  if (sum(is.na(saemixObject["data"]["data"][, saemixObject["data"]["name.response"]])) > 0) {
    if (saemixObject["model"]["modeltype"] == "likelihood") {
      if (saemixObject["options"]$warnings) cat("Please provide values of the response to obtain predictions for a model defined by loglikelihood\n")
      return(NULL)
    }
    type <- type[type != "ipred" & type != "icpred"]
  }
  if (length(type) == 0) type <- "ppred"

  saemixObject <- saemix::estimateMeanParametersNewdata(saemixObject)
  newdata <- saemixObject["data"]
  chdat   <- saemixObject["rep.data"]
  NM   <- chdat["NM"]
  IdM  <- chdat["dataM"]$IdM
  yM   <- chdat["dataM"]$yM
  XM   <- chdat["dataM"][, c(newdata["name.predictors"], newdata["name.cens"],
                             newdata["name.mdv"], newdata["name.ytype"]), drop = FALSE]
  mean.phi <- saemixObject["results"]["mean.phi"]
  psiM <- saemix::transphi(mean.phi, saemixObject["model"]["transform.par"])
  fpred <- saemixObject["model"]["model"](psiM, IdM, XM)
  colnames(psiM) <- saemixObject["model"]["name.modpar"]
  predictions <- data.frame(IdM, XM, ppred = fpred)
  colnames(predictions)[1] <- newdata["name.group"]
  parameters <- list(id = unique(newdata["data"][, newdata["name.group"]]), population = psiM)
  saemixObject["results"]["ppred"] <- fpred

  if (length(grep("ypred", type)) > 0) {
    ind.eta <- saemixObject["model"]["indx.omega"]
    nb.etas <- length(ind.eta)
    omega   <- saemixObject["results"]["omega"]
    chol.omega <- try(chol(omega[ind.eta, ind.eta]), silent = TRUE)
    etaM   <- matrix(0, nrow = NM, ncol = nb.etas)
    ypred  <- matrix(0, nrow = dim(XM)[1], ncol = saemixObject["options"]$nb.sim)
    mean.phiM <- do.call(rbind, rep(list(mean.phi), 1))
    phiMc <- mean.phiM
    for (isim in 1:saemixObject["options"]$nb.sim) {
      etaMc <- 0.5 * matrix(rnorm(NM * nb.etas), ncol = nb.etas) %*% chol.omega
      phiMc[, ind.eta] <- mean.phiM[, ind.eta] + etaMc
      psiMc <- saemix::transphi(phiMc, saemixObject["model"]["transform.par"])
      fpred <- saemixObject["model"]["model"](psiMc, IdM, XM)
      ypred[, isim] <- fpred
    }
    ypred <- rowMeans(ypred)
    predictions$ypred <- ypred
    saemixObject["results"]["ypred"] <- ypred
  }
  if (sum(!is.na(match(c("icpred","ipred"), type))) == 0) {
    return(list(param = parameters, predictions = predictions, object = saemixObject))
  }

  ctype <- c()
  if (length(grep("ipred", type)) == 1)  ctype <- c(ctype, "mode")
  if (length(grep("icpred", type)) == 1) ctype <- c(ctype, "mean")
  saemixObject <- saemix::estimateIndividualParametersNewdata(saemixObject, type = ctype, nsamp = nsamp)

  if (length(grep("icpred", type)) == 1) {
    psiM <- parameters$cond.mean.psi <- saemixObject["results"]["cond.mean.psi"]
    parameters$cond.var.phi  <- saemixObject["results"]["cond.var.phi"]
    parameters$cond.mean.phi <- saemixObject["results"]["cond.mean.phi"]
    fpred <- saemixObject["model"]["model"](psiM, IdM, XM)
    predictions <- cbind(predictions, icpred = fpred)
    saemixObject["results"]["icpred"] <- fpred
    if (nsamp > 1) {
      samp.pred <- array(dim = c(length(fpred), nsamp))
      samp.par  <- array(dim = c(dim(psiM), nsamp))
      for (isamp in 1:nsamp) {
        phiM <- saemixObject["results"]["phi.samp"][, , isamp]
        if (is.null(dim(phiM))) phiM <- matrix(phiM, nrow = 1)
        psiM <- samp.par[, , isamp] <- saemix::transphi(phiM, saemixObject["model"]["transform.par"])
        fpred <- saemixObject["model"]["model"](psiM, IdM, XM)
        samp.pred[, isamp] <- fpred
      }
    }
  }
  if (length(grep("ipred", type)) == 1) {
    psiM <- parameters$map.psi <- saemixObject["results"]["map.psi"]
    fpred <- saemixObject["model"]["model"](psiM, IdM, XM)
    predictions <- cbind(predictions, ipred = fpred)
    saemixObject["results"]["ipred"] <- fpred
  }
  saemixObject["results"]["predictions"] <- predictions
  rlist <- list(param = parameters, predictions = predictions, object = saemixObject)
  if (length(grep("icpred", type)) == 1 & nsamp > 1)
    rlist <- list(param = parameters, predictions = predictions, object = saemixObject,
                  parSample = samp.par, predSample = samp.pred)
  rlist
}

# Number of MCMC samples drawn from the conditional posterior
# (Levoshiny standard = 200; ici 1000 pour distribution plus fine)
FLUCO_NSAMP_COND <- 1000L

fluco_current_dose_times   <- numeric(0)
fluco_current_dose_amounts <- numeric(0)
fluco_current_tinf         <- 1
fluco_saemix_fit_global    <- NULL  # build on first use

fluco_modelFluconazole <- function(psi, id, xidep) {
  CL_v <- psi[id, 1]; V_v <- psi[id, 2]
  t_v  <- xidep[, 1]
  n <- length(t_v); ypred <- numeric(n)
  dt_g <- fluco_current_dose_times; da_g <- fluco_current_dose_amounts; tinf <- fluco_current_tinf
  for (i in seq_len(n)) {
    CL <- CL_v[i]; V <- V_v[i]; t_obs <- t_v[i]; k <- CL / V
    c_tot <- 0
    for (j in seq_along(dt_g)) {
      td <- dt_g[j]; D <- da_g[j]; R <- D / tinf
      dt_diff <- t_obs - td
      if (dt_diff <= 0) next
      if (dt_diff <= tinf) c_tot <- c_tot + (R / CL) * (1 - exp(-k * dt_diff))
      else                 c_tot <- c_tot + (R / CL) * (1 - exp(-k * tinf)) * exp(-k * (dt_diff - tinf))
    }
    ypred[i] <- max(c_tot, 1e-12)
  }
  ypred
}

fluco_build_saemix_fit <- function() {
  if (!requireNamespace("saemix", quietly = TRUE)) {
    install.packages("saemix", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  fluco_current_dose_times   <<- 0
  fluco_current_dose_amounts <<- 400
  fluco_current_tinf         <<- 1
  fake_data <- data.frame(
    Id = c(1,1,2,2,3,3,4,4,5,5),
    Time = c(2, 22, 2, 22, 2, 22, 2, 22, 2, 22),
    Concentration = c(15, 5, 12, 4, 18, 6, 10, 3, 16, 5.5),
    stringsAsFactors = FALSE)
  saemix.data <- saemix::saemixData(name.data = fake_data, header = TRUE, sep = " ", na = NA,
                                     name.group = "Id", name.predictors = "Time",
                                     name.response = "Concentration",
                                     units = list(x = "h", y = "mg/L"),
                                     name.X = "Time", verbose = FALSE)
  saemix.model <- saemix::saemixModel(
    model = fluco_modelFluconazole, modeltype = "structural",
    description = "Fluconazole 1cpt IV (Bendjilali-Sabiani 2026)",
    psi0 = matrix(c(CLpop, Vpop), nrow = 1, ncol = 2,
                  dimnames = list(NULL, c("CL", "V"))),
    transform.par = c(1, 1),
    covariance.model = diag(2),
    omega.init = diag(c(omega_CL^2, omega_V^2)),
    error.model = "proportional",
    error.init = c(0, 0.20),
    verbose = FALSE)
  saemix.options <- list(seed = 632545, save = FALSE, save.graphs = FALSE,
                         nbiter.saemix = c(5, 5), fim = FALSE, map = FALSE,
                         displayProgress = FALSE, print = FALSE)
  fit0 <- suppressMessages(suppressWarnings(saemix::saemix(saemix.model, saemix.data, saemix.options)))
  fit0@results@fixed.psi     <- c(CLpop, Vpop)
  fit0@results@fixed.effects <- c(CLpop, Vpop)
  fit0@results@omega         <- diag(c(omega_CL^2, omega_V^2))
  fit0@results@respar        <- c(0, 0.20)
  fit0@results@betas         <- t(c(log(CLpop), log(Vpop)))
  fit0 <- saemix::estimateMeanParametersNewdata(fit0)
  fit0@options$warnings        <- FALSE
  fit0@options$displayProgress <- FALSE
  fit0
}

fluco_saemix_estimate <- function(nmtran) {
  if (is.null(nmtran) || nrow(nmtran) == 0) stop("Empty NM-TRAN")
  doses <- nmtran[nmtran$evid == 1, , drop = FALSE]
  obs   <- nmtran[nmtran$evid == 0, , drop = FALSE]
  if (nrow(obs) == 0) stop("No observation in NM-TRAN")
  if (nrow(doses) == 0) stop("No dose in NM-TRAN")

  bw_v    <- median(nmtran$WT, na.rm = TRUE);    if (!is.finite(bw_v))    bw_v    <- 70.5
  ckd_v   <- median(nmtran$CKD, na.rm = TRUE);   if (!is.finite(ckd_v))   ckd_v   <- 73.7
  cirrh_v <- median(nmtran$CIRRH, na.rm = TRUE); if (!is.finite(cirrh_v)) cirrh_v <- 0
  bw_v    <- max(WT_LOWER,  min(bw_v,  WT_UPPER))
  ckd_v   <- max(CKD_LOWER, min(ckd_v, CKD_UPPER))
  cirrh_v <- ifelse(cirrh_v >= 0.5, 1, 0)

  d_times <- doses$time
  d_amts  <- doses$amt
  d_rate  <- doses$rate
  d_tinf  <- ifelse(d_rate > 0, d_amts / d_rate, 1)
  tinf_typ <- median(d_tinf[d_tinf > 0])
  if (!is.finite(tinf_typ) || tinf_typ <= 0) tinf_typ <- 1

  CL_typ <- predict_CL_typ(ckd_v, cirrh_v)
  V_typ  <- predict_V_typ(bw_v)

  # NOTE IMPORTANTE: les globals fluco_current_dose_times/amounts/tinf sont
  # ECRASES par fluco_build_saemix_fit() (qui les met a 0/400/1 pour son
  # fit bidon). Il FAUT donc les ecrire APRES la construction du fit, juste
  # avant les appels saemix qui utilisent la fonction modele.

  newdata <- data.frame(Id = 1L, Time = obs$time, Concentration = obs$DV,
                        stringsAsFactors = FALSE)

  # Construction d'un fit FRAIS a chaque appel (evite toute contamination
  # entre estimations successives qui partageraient un meme objet saemix)
  fit_pat <- fluco_build_saemix_fit()

  # MAINTENANT (et pas avant!) on ecrit les vraies doses du patient
  fluco_current_dose_times   <<- d_times
  fluco_current_dose_amounts <<- d_amts
  fluco_current_tinf         <<- tinf_typ

  fit_pat@results@fixed.psi     <- c(CL_typ, V_typ)
  fit_pat@results@fixed.effects <- c(CL_typ, V_typ)
  fit_pat@results@betas         <- t(c(log(CL_typ), log(V_typ)))
  fit_pat <- suppressMessages(suppressWarnings(saemix::estimateMeanParametersNewdata(fit_pat)))

  # ---- Prediction bayesienne facon Levoshiny: nsamp echantillons MCMC --
  # de la distribution conditionnelle a posteriori (vraie posterior, pas
  # un sampling normal autour du MAP via la matrice OMEGA prior)
  xpred <- tryCatch(
    suppressMessages(suppressWarnings(
      saemixPredictNewdata(fit_pat, newdata,
                           type = c("ipred","ypred","ppred","icpred"),
                           nsamp = FLUCO_NSAMP_COND))),
    error = function(e) { message("saemixPredictNewdata failed: ", e$message); NULL })

  # ---- Recuperation du MAP (mode de la conditionnelle) --------------------
  CL_map <- NA_real_; V_map <- NA_real_
  if (!is.null(xpred) && !is.null(xpred$param$map.psi)) {
    psi_map <- xpred$param$map.psi
    if (is.matrix(psi_map)) { CL_map <- as.numeric(psi_map[1, 1]); V_map <- as.numeric(psi_map[1, 2]) }
    else                    { CL_map <- as.numeric(psi_map[1]);    V_map <- as.numeric(psi_map[2])    }
  }
  if (!is.finite(CL_map) || CL_map <= 0) CL_map <- CL_typ
  if (!is.finite(V_map)  || V_map  <= 0) V_map  <- V_typ

  # ---- Recuperation de la moyenne conditionnelle --------------------------
  CL_cm <- CL_map; V_cm <- V_map
  if (!is.null(xpred) && !is.null(xpred$param$cond.mean.psi)) {
    psi_cm <- xpred$param$cond.mean.psi
    if (is.matrix(psi_cm)) { CL_cm <- as.numeric(psi_cm[1, 1]); V_cm <- as.numeric(psi_cm[1, 2]) }
    else                   { CL_cm <- as.numeric(psi_cm[1]);    V_cm <- as.numeric(psi_cm[2])    }
    if (!is.finite(CL_cm) || CL_cm <= 0) CL_cm <- CL_map
    if (!is.finite(V_cm)  || V_cm  <= 0) V_cm  <- V_map
  }

  # ---- Extraction des nsamp echantillons MCMC pour CL et V ----------------
  # parSample est un array de dimension [N=1, n_param=2, nsamp]
  CL_samp <- rep(CL_map, FLUCO_NSAMP_COND)
  V_samp  <- rep(V_map,  FLUCO_NSAMP_COND)
  if (!is.null(xpred) && !is.null(xpred$parSample)) {
    pr <- xpred$parSample
    # pr[1, 1, ] = CL pour patient 1 sur les nsamp echantillons
    # pr[1, 2, ] = V  pour patient 1 sur les nsamp echantillons
    CL_samp <- as.numeric(pr[1, 1, ])
    V_samp  <- as.numeric(pr[1, 2, ])
    CL_samp[!is.finite(CL_samp) | CL_samp <= 0] <- CL_map
    V_samp[!is.finite(V_samp)   | V_samp  <= 0] <- V_map
  }

  list(
    param = c(CL = CL_map, V = V_map),
    cond_mean = c(CL = CL_cm, V = V_cm),
    eta = c(ETA1 = log(CL_map / CL_typ), ETA2 = log(V_map / V_typ)),
    CLtyp = CL_typ, Vtyp = V_typ,
    CL_samp = CL_samp, V_samp = V_samp,   # distribution conditionnelle MCMC
    n_cond  = length(CL_samp),
    opt.value = list(hessian = NULL),
    .saemix_compat = TRUE
  )
}

fluco_get_CL <- function(est) {
  if (is.list(est) && isTRUE(est$.saemix_compat)) return(as.numeric(est$param[["CL"]]))
  stop("fluco_get_CL: object not recognized")
}
fluco_get_V <- function(est) {
  if (is.list(est) && isTRUE(est$.saemix_compat)) return(as.numeric(est$param[["V"]]))
  stop("fluco_get_V: object not recognized")
}
fluco_get_etas <- function(est) {
  if (is.list(est) && isTRUE(est$.saemix_compat)) return(est$eta)
  stop("fluco_get_etas: object not recognized")
}
# Recupere les nsamp echantillons MCMC de la distribution conditionnelle
# Retourne une liste(CL = vecteur, V = vecteur) ou NULL si indisponible
fluco_get_cond_samples <- function(est) {
  if (is.list(est) && isTRUE(est$.saemix_compat) &&
      !is.null(est$CL_samp) && !is.null(est$V_samp)) {
    return(list(CL = est$CL_samp, V = est$V_samp))
  }
  NULL
}
# ==============================================================================

# ==============================================================================
# A POSTERIORI - Bayesian estimation using saemix (conditional posterior)
# ==============================================================================

# Build NM-TRAN dataset from Shiny inputs (like prof's data.frame)
# Time 0 = datetime of first dose. All times in hours relative to t0.
build_nmtran_data <- function(date_first, heure_first, date_maint, heure_maint, date_last, dose_first, dose_maint, tau, tinf, obs_full, wt, ckd, cirrh) {
  wt <- max(WT_LOWER, min(wt, WT_UPPER)); ckd <- max(CKD_LOWER, min(ckd, CKD_UPPER)); tinf <- max(0.05, tinf)
  # Convert all to POSIXct
  t0 <- as.POSIXct(paste(as.Date(date_first, "%d/%m/%Y"), sprintf("%02d:%02d:00", floor(heure_first), round((heure_first%%1)*60))), format="%Y-%m-%d %H:%M:%S")
  t_maint <- as.POSIXct(paste(as.Date(date_maint, "%d/%m/%Y"), sprintf("%02d:%02d:00", floor(heure_maint), round((heure_maint%%1)*60))), format="%Y-%m-%d %H:%M:%S")
  t_last <- as.POSIXct(paste(as.Date(date_last, "%d/%m/%Y"), sprintf("%02d:%02d:00", floor(heure_maint), round((heure_maint%%1)*60))), format="%Y-%m-%d %H:%M:%S")
  # Loading dose at t=0
  dose_times_h <- 0; dose_amounts <- dose_first
  # Maintenance doses from t_maint every tau hours until t_last
  h_maint_start <- as.numeric(difftime(t_maint, t0, units="hours"))
  h_last <- as.numeric(difftime(t_last, t0, units="hours"))
  if (h_maint_start > 0 && h_last >= h_maint_start) {
    mt <- seq(h_maint_start, h_last, by=tau)
    dose_times_h <- c(dose_times_h, mt)
    dose_amounts <- c(dose_amounts, rep(dose_maint, length(mt)))
  }
  # Time-varying CKD at dose times: linear interpolation between obs CKD values
  # If obs_full has CKD per obs, interpolate at each dose time
  if (!is.null(obs_full) && nrow(obs_full)>0 && !is.null(obs_full$ckd) && any(is.finite(obs_full$ckd))) {
    obs_ckd_t <- obs_full$abs_h; obs_ckd_v <- obs_full$ckd
    valid_ckd <- is.finite(obs_ckd_v)
    if (sum(valid_ckd) >= 2) {
      ckd_at_dose <- approx(obs_ckd_t[valid_ckd], obs_ckd_v[valid_ckd], xout=dose_times_h, rule=2)$y
    } else {
      ckd_at_dose <- rep(if(sum(valid_ckd)>=1) obs_ckd_v[valid_ckd][1] else ckd, length(dose_times_h))
    }
    ckd_at_dose <- pmax(CKD_LOWER, pmin(ckd_at_dose, CKD_UPPER))
  } else {
    ckd_at_dose <- rep(ckd, length(dose_times_h))
  }
  dose_rows <- data.frame(
    ID = 1, time = dose_times_h, evid = 1, cmt = 1,
    amt = dose_amounts, rate = dose_amounts / tinf,
    mdv = 1, DV = 0,
    CKD = ckd_at_dose, WT = wt, CIRRH = cirrh,
    stringsAsFactors = FALSE)
  if (!is.null(obs_full) && nrow(obs_full) > 0) {
    # Time-varying CKD: use per-observation CKD if available, else global ckd
    ckd_per_obs <- if(!is.null(obs_full$ckd) && all(is.finite(obs_full$ckd))) obs_full$ckd else rep(ckd, nrow(obs_full))
    obs_rows <- data.frame(
      ID = 1, time = obs_full$abs_h, evid = 0, cmt = 1,
      amt = 0, rate = 0, mdv = 0, DV = obs_full$conc,
      CKD = ckd_per_obs, WT = wt, CIRRH = cirrh,
      stringsAsFactors = FALSE)
  } else { return(NULL) }
  nmtran <- rbind(dose_rows, obs_rows)
  nmtran <- nmtran[order(nmtran$time, -nmtran$evid), ]
  rownames(nmtran) <- NULL
  message("=== NM-TRAN (",nrow(dose_rows)," doses, ",nrow(obs_rows)," obs) ==="); print(nmtran)
  nmtran
}

estimate_aposteriori <- function(obs_full, wt, ckd, cirrh, total_daily_dose, n_admin, tinf,
                                  date_first, heure_first, date_maint, heure_maint, date_last, tau, dose_first, dose_maint) {
  # CRITICAL: clamp to model development range to prevent mrgsolve segfault
  wt <- max(WT_LOWER, min(wt, WT_UPPER))   # 42.1 - 123 kg
  ckd <- max(CKD_LOWER, min(ckd, CKD_UPPER))  # 15 - 175 mL/min
  tinf <- max(0.05, tinf)
  CLtyp <- predict_CL_typ(ckd, cirrh); Vtyp <- predict_V_typ(wt)
  dose <- dose_maint
  nmtran <- build_nmtran_data(date_first, heure_first, date_maint, heure_maint, date_last, dose_first, dose_maint, tau, tinf, obs_full, wt, ckd, cirrh)
  if (is.null(nmtran)) stop("Pas de donnees NM-TRAN valides")
  message("=== NM-TRAN data ==="); print(nmtran)

  # MAP estimation - wrapped (saemix engine)
  est <- tryCatch(fluco_saemix_estimate(nmtran), error = function(e) {
    message("saemix failed: ", e$message); stop(paste0("saemix: ", e$message))
  })
  message("=== saemix OK ===")
  CL_map <- tryCatch(fluco_get_CL(est), error = function(e) CLtyp)
  V_map <- tryCatch(fluco_get_V(est), error = function(e) Vtyp)
  etas <- tryCatch(fluco_get_etas(est), error = function(e) c(ETA1=0, ETA2=0))
  eta1 <- as.numeric(etas[1]); eta2 <- as.numeric(etas[2])
  if(!is.finite(eta1)) eta1 <- 0; if(!is.finite(eta2)) eta2 <- 0
  # Conditional mean (Levoshiny: la moyenne conditionnelle, distincte du MAP)
  CL_cm <- if (!is.null(est$cond_mean)) as.numeric(est$cond_mean["CL"]) else CL_map
  V_cm  <- if (!is.null(est$cond_mean)) as.numeric(est$cond_mean["V"])  else V_map
  if (!is.finite(CL_cm) || CL_cm <= 0) CL_cm <- CL_map
  if (!is.finite(V_cm)  || V_cm  <= 0) V_cm  <- V_map

  # ---- Posterior sampling — APPROCHE LEVOSHINY -----------------------------
  # Distribution conditionnelle a posteriori issue de SAEMIX (nsamp=1000).
  # On reechantillonne avec remise pour atteindre NPOST_SAMP afin de garder
  # la compatibilite avec le code aval (rapports, plots, simulations PTA).
  cond_s <- fluco_get_cond_samples(est)
  if (!is.null(cond_s) && length(cond_s$CL) > 0 && length(cond_s$V) > 0) {
    set.seed(12345)
    idx <- sample.int(length(cond_s$CL), size = NPOST_SAMP, replace = TRUE)
    cl_v <- pmax(as.numeric(cond_s$CL[idx]), 0.001)
    v_v  <- pmax(as.numeric(cond_s$V[idx]),  0.1)
    post <- data.frame(
      etaCL = log(cl_v / CLtyp),
      etaV  = log(v_v  / Vtyp),
      CL    = cl_v,
      V     = v_v
    )
  } else {
    # Fallback (saemix sampling indispo): ancien comportement (sampling normal)
    samp_eta <- tryCatch(rmvnorm_chol(NPOST_SAMP, c(eta1, eta2), OMEGA_MAT), error = function(e) {
      matrix(rnorm(NPOST_SAMP*2, sd=0.2), ncol=2)
    })
    post <- data.frame(etaCL = samp_eta[, 1], etaV = samp_eta[, 2])
    post$CL <- pmax(CLtyp * exp(post$etaCL), 0.001)
    post$V  <- pmax(Vtyp  * exp(post$etaV),  0.1)
  }
  dose_times_auc <- seq(0, by = tau, length.out = 14)
  dose_amounts_auc <- rep(dose, length(dose_times_auc))
  post$AUC24 <- mapply(function(cl, v) {
    auc_model_interval((length(dose_times_auc)-1)*tau, length(dose_times_auc)*tau,
                       dose_times_auc, dose_amounts_auc, cl, v, tinf)
  }, post$CL, post$V)

  # SS profile - wrapped
  # Posterior credible band: 500 profiles simulated from MCMC samples
  # (post$CL, post$V come from saemix conditional distribution).
  # Vectorized via idata_set() for efficiency (single mrgsim call).
  prof <- tryCatch({
    tgrid <- seq(0, tau, length.out = POST_PROFILE_GRID)

    # ---- (a) MAP profile (central line) -------------------------------------
    mod_patient <- mod_fluco %>%
      param(ETA1 = eta1, ETA2 = eta2, CKD = ckd, WT = wt, CIRRH = cirrh) %>%
      zero_re()
    sim_prof <- mod_patient %>%
      ev(amt = dose, rate = dose / tinf, ii = tau, addl = 13, cmt = 1) %>%
      mrgsim(end = 14 * tau, delta = tau / POST_PROFILE_GRID) %>%
      as.data.frame()
    ss_start <- 13 * tau
    sim_ss   <- sim_prof[sim_prof$time >= ss_start &
                         sim_prof$time <= (ss_start + tau), ]
    dup      <- !duplicated(sim_ss$time - ss_start)
    prof_map <- approx((sim_ss$time - ss_start)[dup],
                       sim_ss$CP[dup], xout = tgrid, rule = 2)$y
    if (any(is.na(prof_map))) prof_map[is.na(prof_map)] <- 0

    # ---- (b) 95% credible band from posterior MCMC samples ------------------
    N_BAND   <- min(500L, nrow(post))
    set.seed(12345)
    band_idx <- sample.int(nrow(post), N_BAND)
    cl_b     <- post$CL[band_idx]
    v_b      <- post$V[band_idx]
    eta_cl_b <- log(pmax(cl_b, 1e-6) / max(CLtyp, 1e-6))
    eta_v_b  <- log(pmax(v_b,  1e-6) / max(Vtyp,  1e-6))

    idata_band <- data.frame(
      ID    = seq_len(N_BAND),
      ETA1  = eta_cl_b,
      ETA2  = eta_v_b,
      WT    = wt,
      CKD   = ckd,
      CIRRH = cirrh
    )
    band_sim <- mod_fluco %>%
      zero_re() %>%
      idata_set(idata_band) %>%
      ev(amt = dose, rate = dose / tinf, ii = tau, addl = 13, cmt = 1) %>%
      mrgsim(end = 14 * tau, delta = tau / POST_PROFILE_GRID) %>%
      as.data.frame()
    band_ss <- band_sim[band_sim$time >= ss_start &
                        band_sim$time <= (ss_start + tau), ]

    band_mat <- vapply(seq_len(N_BAND), function(j) {
      sub <- band_ss[band_ss$ID == j, ]
      if (nrow(sub) < 2) return(rep(NA_real_, length(tgrid)))
      dup_j <- !duplicated(sub$time - ss_start)
      approx((sub$time - ss_start)[dup_j], sub$CP[dup_j],
             xout = tgrid, rule = 2)$y
    }, numeric(length(tgrid)))

    q025 <- apply(band_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
    q500 <- apply(band_mat, 1, quantile, probs = 0.500, na.rm = TRUE)
    q975 <- apply(band_mat, 1, quantile, probs = 0.975, na.rm = TRUE)

    data.frame(time = tgrid, map = prof_map, condmean = prof_map,
               p2.5 = q025, p50 = q500, p97.5 = q975)
  }, error = function(e) {
    message("SS profile (credible band) failed: ", e$message)
    tgrid <- seq(0, tau, length.out = POST_PROFILE_GRID)
    data.frame(time = tgrid, map = rep(0, length(tgrid)), condmean = rep(0, length(tgrid)),
               p2.5 = rep(0, length(tgrid)), p50 = rep(0, length(tgrid)), p97.5 = rep(0, length(tgrid)))
  })

  mod_patient <- tryCatch(
    mod_fluco %>% param(ETA1=eta1, ETA2=eta2, CKD=ckd, WT=wt, CIRRH=cirrh) %>% zero_re(),
    error = function(e) NULL)

  obs_compat <- data.frame(time = obs_full$time, conc = obs_full$conc)
  list(CLtyp = CLtyp, Vtyp = Vtyp, dose = dose, tau = tau, tinf = tinf, obs = obs_compat,
       map = c(etaCL = eta1, etaV = eta2),
       map_psi = c(CL = CL_map, V = V_map), cond_mean_psi = c(CL = CL_cm, V = V_cm),
       posterior = post, profile = prof, convergence = 0, objective = 0,
       nmtran = nmtran, saemix_est = est, mod_patient = mod_patient)
}

# ==============================================================================
# HTML REPORT - shared CSS
# ==============================================================================
report_css <- function() {
'<style>
@import url("https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700;800;900&display=swap");
body { font-family: Inter, Arial, sans-serif; margin: 0; padding: 0; color: #1a1a2e; background: #f5f7fa; line-height: 1.5; }
.title { display: none; }
.chu-header { display: none; }

/* ── REPORT HEADER BANNER ── */
.report-banner {
  background: linear-gradient(135deg, #003a75 0%, #00539F 60%, #0077cc 100%);
  max-width: 100%; margin: 0 0 0 0; padding: 32px 40px 26px 40px;
  position: relative; overflow: hidden;
}
.report-banner::after {
  content: ""; position: absolute; bottom: 0; left: 0; right: 0; height: 4px;
  background: linear-gradient(to right, #FFD700, #FFA500, #FFD700);
}
.report-banner .rpt-label {
  font-size: 11px; font-weight: 700; letter-spacing: 2.5px; color: rgba(255,255,255,0.65);
  text-transform: uppercase; margin: 0 0 8px 0;
}
.report-banner h1 {
  font-size: 28px; font-weight: 900; color: #ffffff; margin: 0 0 10px 0;
  letter-spacing: -0.5px; line-height: 1.15;
}
.report-banner .rpt-drug {
  font-size: 16px; font-weight: 600; color: rgba(255,255,255,0.92); margin: 0 0 6px 0;
}
.report-banner .rpt-engine {
  display: inline-block; background: rgba(255,255,255,0.15); border: 1px solid rgba(255,255,255,0.3);
  border-radius: 20px; padding: 3px 12px; font-size: 12px; font-weight: 600;
  color: rgba(255,255,255,0.90); margin: 0 0 12px 0; letter-spacing: 0.3px;
}
.report-banner .rpt-date {
  font-size: 13px; font-weight: 400; color: rgba(255,255,255,0.70);
  border-top: 1px solid rgba(255,255,255,0.20); padding-top: 10px; margin-top: 4px;
}
.report-banner .rpt-date strong { color: rgba(255,255,255,0.95); font-weight: 700; }

/* ── PAGE BODY ── */
.page-body { max-width: 900px; margin: 0 auto; padding: 24px 0 24px 0; }

/* ── PATIENT ROW ── */
.patient-row { display: flex; gap: 14px; margin: 0 0 18px 0; }
.patient-box { flex: 1; background: white; border: 1.5px solid #d0d8e8; border-radius: 10px; padding: 14px 18px; box-shadow: 0 1px 4px rgba(0,83,159,0.07); }
.patient-box .pat-name { font-size: 16px; font-weight: 800; color: #00539F; margin-bottom: 6px; }
.patient-box .pat-line { font-size: 13.5px; line-height: 1.85; color: #2a2a3a; }
.patient-box .pat-label { font-weight: 700; color: #444; }
.service-box { flex: 1; background: white; border: 1.5px solid #d0d8e8; border-radius: 10px; padding: 14px 18px; display: flex; flex-direction: column; justify-content: center; box-shadow: 0 1px 4px rgba(0,83,159,0.07); }
.service-box .svc-line { font-size: 13.5px; color: #2a2a3a; text-align: left; line-height: 1.7; }
.service-box .svc-line:first-child { font-size: 15px; font-weight: 700; color: #00539F; margin-bottom: 4px; }
.service-box .svc-addr { font-size: 13px; color: #666; }

/* ── SECTIONS ── */
.section { margin: 0 0 16px 0; }
.section h2 { color: #00539F; font-size: 15px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.8px; border-bottom: 2px solid #00539F; padding-bottom: 5px; margin-top: 24px; margin-bottom: 12px; }
.param-box { background: white; border: 1.5px solid #d0d8e8; border-radius: 10px; padding: 14px 18px; margin-bottom: 14px; font-size: 14px; line-height: 1.8; box-shadow: 0 1px 4px rgba(0,83,159,0.06); }
.param-box .lab { font-weight: 700; color: #003a75; }

/* ── REC BOXES ── */
.rec-ok { background: #f0faf2; border-left: 5px solid #2E7D32; border-radius: 10px; padding: 14px 18px; margin-bottom: 10px; font-size: 14px; box-shadow: 0 1px 3px rgba(46,125,50,0.10); }
.rec-warn { background: #fffbf0; border-left: 5px solid #FB8C00; border-radius: 10px; padding: 14px 18px; margin-bottom: 10px; font-size: 14px; box-shadow: 0 1px 3px rgba(251,140,0,0.10); }
.rec-fail { background: #fff5f5; border-left: 5px solid #C62828; border-radius: 10px; padding: 14px 18px; margin-bottom: 10px; font-size: 14px; box-shadow: 0 1px 3px rgba(198,40,40,0.10); }
.alert-box { background: #fff5f5; border-left: 5px solid #C62828; border-radius: 10px; padding: 12px 16px; margin-bottom: 12px; font-weight: 600; color: #C62828; font-size: 14px; }

/* ── TABLES ── */
table { width: 100%; border-collapse: collapse; margin: 10px 0; font-size: 13px; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 4px rgba(0,0,0,0.07); }
th { background: #00539F; color: white; padding: 10px 7px; text-align: center; font-weight: 600; font-size: 12px; letter-spacing: 0.3px; }
td { padding: 8px 7px; text-align: center; border-bottom: 1px solid #edf0f5; background: white; }
tr:nth-child(even) td { background: #f7f9fd; }

/* ── FOOTER ── */
.footer { margin: 28px 0 0 0; border-top: 1px solid #dde3ee; padding-top: 10px; font-size: 10.5px; color: #999; line-height: 1.6; }
</style>'
}

report_header_chunk <- function() {
''
}

report_patient_block <- function() {
'```{r patient_block, results="asis"}
age_txt <- tryCatch({
  age <- as.numeric(difftime(Sys.Date(), as.Date(clin$dob, "%d/%m/%Y"), units="days")) / 365.25
  paste0(round(age), " years")
}, error = function(e) "--")
cirrh_txt <- if (!is.null(clin$cirrh) && clin$cirrh == 1) "Yes" else "No"
sexe_code <- if (!is.null(p$sexe) && p$sexe == "Female") "F" else "M"
svc_lines <- paste0("<div class=\\"svc-line\\">", p$service, "</div>")
if (!is.null(p$uf) && nzchar(p$uf)) svc_lines <- paste0(svc_lines, "<div class=\\"svc-line\\">UF : ", p$uf, "</div>")
if (!is.null(p$hopital) && nzchar(p$hopital)) svc_lines <- paste0(svc_lines, "<div class=\\"svc-line\\">", p$hopital, "</div>")
if (!is.null(p$adresse) && nzchar(p$adresse)) svc_lines <- paste0(svc_lines, "<div class=\\"svc-line svc-addr\\">", p$adresse, "</div>")
cat(paste0(
  "<div class=\\"patient-row\\">",
  "<div class=\\"patient-box\\">",
  "<div class=\\"pat-name\\">", p$last_name, " ", p$first_name, "</div>",
  "<div class=\\"pat-line\\"><span class=\\"pat-label\\">Date of birth:</span> ", p$dob, " (", age_txt, ")</div>",
  "<div class=\\"pat-line\\"><span class=\\"pat-label\\">IPP :</span> ", p$ipp, "</div>",
  "<div class=\\"pat-line\\"><span class=\\"pat-label\\">Sex:</span> ", sexe_code, "</div>",
  "<div class=\\"pat-line\\"><span class=\\"pat-label\\">Weight:</span> ", clin$wt, " kg</div>",
  "<div class=\\"pat-line\\"><span class=\\"pat-label\\">CKD-EPI :</span> ", clin$ckd, " mL/min/1.73m²</div>",
  "<div class=\\"pat-line\\"><span class=\\"pat-label\\">Cirrhosis:</span> ", cirrh_txt, "</div>",
  "</div>",
  "<div class=\\"service-box\\">",
  svc_lines,
  "</div></div>"
))
```'
}

report_footer <- function() {
'```{r footer_valid, results="asis"}
pharma_txt <- if (!is.null(params$pharmacologue) && nzchar(params$pharmacologue)) params$pharmacologue else "___________________"
cat(paste0(
  "<div style=\\"max-width:900px;margin:20px auto;padding:14px 18px;border:2px solid #00539F;border-radius:8px;font-size:14px;\\">",
  "<b>Validated by</b> ", pharma_txt, ", Clinical Pharmacologist<br>",
  "<b>Validation date:</b> ", format(Sys.Date(), "%d/%m/%Y"),
  "</div>"))
```
<div style="max-width:900px;margin:20px auto;padding:12px 16px;border:1px solid #B71C1C;border-left:5px solid #B71C1C;border-radius:6px;background:#FFF8F8;font-size:12px;color:#444;line-height:1.5;">
<b style="color:#B71C1C;">DISCLAIMER / AVERTISSEMENT</b><br>
This software is provided "AS IS", without warranty of any kind. It is a research prototype intended to support &mdash; never to replace &mdash; the clinical judgment of a qualified physician or pharmacologist. All dose recommendations must be reviewed, validated, and adjusted by a competent healthcare professional taking into account the full clinical context of the patient. The authors and their affiliated institutions decline all responsibility for any adverse outcome, treatment failure, injury, or death that may occur in connection with the use of this software or any decision based on its outputs.
</div>
<div class="footer">
<b>PK Model:</b> 1-compartment IV, adult digestive ICU population.
CL = 0.70&middot;(CKD/73.7)<sup>0.88</sup>&middot;exp(&minus;0.44&middot;cirrh) L/h ;
V = 63.74&middot;(WT/70.5)<sup>1.66</sup> L.
&omega;<sub>CL</sub>=0.23, &omega;<sub>V</sub>=0.30, &rho;=0.30.<br>
<b>Model reference:</b> <em>`r PK_MODEL_REF`</em><br>
<em>PK/PD references:</em> Ashbee <em>et al.</em> JAC 2014 ; G&oacute;mez-L&oacute;pez. CMI 2020 ; EUCAST v6.0, 2025.<br>
<em>Report generated by</em> <b>FlucoShiny</b> &mdash; Clinical Pharmacology Unit.
</div>'
}

# ==============================================================================
# REPORT WRITERS
# ==============================================================================
write_rmd_apriori <- function(path, logo_path = NULL) {
  con <- file(path, open = "wt", encoding = "UTF-8"); on.exit(close(con))

  cat('---
title: " "
output:
  html_document:
    self_contained: true
    toc: false
    df_print: kable
params:
  patient: NULL
  res: NULL
  logo_path: NULL
  logo_html: ""
  clinical: NULL
  pharmacologue: NULL
  pk_model_ref: ""
  nsim: 10000
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
r <- params$res; p <- params$patient; clin <- params$clinical
PK_MODEL_REF <- params$pk_model_ref
NSIM <- params$nsim
logo_html_str <- if(nzchar(params$logo_html)) params$logo_html else ""
today_txt <- format(Sys.Date(), "%d/%m/%Y")
PTA_OPT <- 0.90
plot_fun <- function(sum, pta, target, labels, rec_idx, ylab, main) {
  x <- seq_along(sum$median)
  cols <- ifelse(pta >= PTA_OPT, "#2E7D32", "#B0BEC5")
  plot(x, sum$median, pch=19, col=cols, xaxt="n",
       ylim=range(c(sum$p2.5, sum$p97.5, target)),
       ylab=ylab, main=main)
  axis(1, at=x, labels=labels)
  arrows(x, sum$p2.5, x, sum$p97.5, angle=90, code=3, length=0.04)
  abline(h=target, col="red", lty=2, lwd=2)
  if (!is.na(rec_idx)) points(rec_idx, sum$median[rec_idx], pch=21, bg="#1565C0", cex=1.6)
}
logo_b64 <- if (!is.null(params$logo_path) && file.exists(params$logo_path)) knitr::image_uri(params$logo_path) else ""
rec_cls <- function(case) switch(case, optimal="rec-ok", suboptimal="rec-warn", "rec-fail")
rec_icon <- function(case) switch(case, optimal="&#10003;", suboptimal="&#9888;", "&#10007;")
```

', file = con)
  cat(report_css(), "\n\n", file = con)
  cat(report_header_chunk(), "\n\n", file = con)

  # Logo chunk written separately to avoid quote conflicts
  cat('<div class="report-banner">
`r logo_html_str`
<div class="rpt-label">Clinical Pharmacology &mdash; PK/PD Simulation Report</div>
<h1>Clinical Pharmacology Report</h1>
<div class="rpt-drug">Fluconazole &mdash; Individual <em>a priori</em> simulation</div>
<div class="rpt-engine">&#9881; Engine: mrgsolve + saemix &nbsp;&nbsp;&#124;&nbsp;&nbsp; n = `r format(NSIM, big.mark=" ")` IIV draws per patient</div>
<div class="rpt-date">Report date: <strong>`r today_txt`</strong></div>
</div>
<div class="page-body">

', file = con)
  cat(report_patient_block(), "\n\n", file = con)

  cat('<div class="section">

<h2>Pharmacokinetic Parameters</h2>
<div class="param-box">
<span class="lab">Target AUC:</span> `r r$target_auc` mg&middot;h/L<br>
<span class="lab">Infusion duration:</span> `r r$tinf` h<br>
<span class="lab">Typical CL:</span> `r round(r$CLtyp, 3)` L/h &nbsp;&nbsp;|&nbsp;&nbsp;
<span class="lab">Typical V:</span> `r round(r$Vtyp, 1)` L<br>
<span class="lab">Target MIC:</span> `r r$mic_label`<br>
<span class="lab">PK model:</span> <em>`r PK_MODEL_REF`</em><br>
<span style="font-size:12px;color:#666;font-style:italic;">Monte Carlo simulation: n = `r format(NSIM, big.mark=" ")` IIV draws (&eta;<sub>CL</sub>, &eta;<sub>V</sub>) on individual patient PK parameters.</span>
</div>

```{r extrap, results="asis"}
if (isTRUE(r$wt_extrap)) cat("<div class=\\"alert-box\\">&#9888; Poids hors plage du mod&egrave;le (42.1&ndash;123 kg).</div>")
if (isTRUE(r$ckd_extrap)) cat("<div class=\\"alert-box\\">&#9888; CKD-EPI hors plage du mod&egrave;le (15&ndash;175 mL/min).</div>")
```

<h2>Dosage Recommendation</h2>

```{r recs, results="asis"}
dose_choice_text <- function(label, wt) {
  if (grepl("mg/kg/day", label, fixed=TRUE)) {
    mgkg <- suppressWarnings(as.numeric(sub(" .*","",label)))
    if (is.na(mgkg)) return(label)
    paste0(label, " (", round(mgkg*wt), " mg/day)")
  } else label
}
load_label <- switch(r$load_case, optimal="Optimal (PTA &ge; 90%)", suboptimal="Suboptimal (PTA &ge; 75%)", "Failure (PTA < 75%)")
maint_label <- switch(r$maint_case, optimal="Optimal (PTA &ge; 90%)", suboptimal="Suboptimal (PTA &ge; 75%)", "Failure (PTA < 75%)")
ld <- if (!is.na(r$rec_load_idx)) dose_choice_text(r$load_df$DoseLabel[r$rec_load_idx], r$wt) else "None"
md <- if (!is.na(r$rec_maint_idx)) dose_choice_text(r$maint_df$DoseLabel[r$rec_maint_idx], r$wt) else "None"
day_labels   <- c("Day 1 (0\u201324h)","Day 2 (24\u201348h)","Day 3 (48\u201372h)")
cum_ranges   <- c("0\u201324h","0\u201348h","0\u201372h")
best_lbl     <- if(!is.null(r$best_load_day)) day_labels[min(r$best_load_day,3L)] else "Day 1 (0\u201324h)"
n_days_rmd   <- if(!is.null(r$best_load_day)) r$best_load_day else 1L
cum_rng      <- cum_ranges[min(n_days_rmd,3L)]
need_repeat  <- !is.null(r$best_load_day) && r$best_load_day > 1
if(need_repeat) {
  subopt_warn <- if(identical(r$load_case,"suboptimal")) paste0(
    "<br><span style=\'color:#E65100;font-size:12px;\'>&#9888; Suboptimal: PTA between 75\u201390% even after ",
    n_days_rmd," days. Consider a higher dose or an alternative antifungal.</span>") else ""
  cat(paste0(
    "<div style=\'padding:12px 16px;border-radius:8px;background:#FFF3CD;",
    "border:1px solid #E65100;border-left:5px solid #E65100;margin-bottom:10px;font-size:13px;\'>",
    "<b style=\'color:#B71C1C;font-size:14px;\'>&#9888;&nbsp;Repeated loading dose required &mdash; ",
    n_days_rmd," consecutive day",if(n_days_rmd>1L)"s" else "","</b><br>",
    "<span style=\'color:#555;\'>",
    "No loading dose achieves PTA &ge;&nbsp;90% for the PK/PD target ",
    "(AUC &ge;&nbsp;",r$target_auc,"&nbsp;mg&middot;h/L) after a <u>single</u> administration (Day 1 alone).",
    "</span><br>",
    "A loading dose of <b>",ld,"</b> repeated over ",
    "<b>",n_days_rmd,"&nbsp;consecutive days</b> (<em>",cum_rng,"</em>) is required.<br>",
    "&#10003;&nbsp;PTA &ge;&nbsp;90% is first achieved on <b>",best_lbl,"</b>.",
    subopt_warn,
    "</div>"
  ))
} else {
  cat(paste0("<div class=\\"",rec_cls(r$load_case),"\\">",rec_icon(r$load_case),
    " <b>Loading dose:</b> ",ld," &mdash; ",load_label,"</div>"))
}
cat(paste0("<div class=\\"", rec_cls(r$maint_case), "\\">", rec_icon(r$maint_case), " <b>Maintenance dose:</b> ", md, " &mdash; ", maint_label, "</div>"))
```

<h2>Loading Dose &mdash; AUC 0&ndash;24 h (Day 1)</h2>

```{r plot_load_d1, fig.height=4.5, fig.width=9}
plot_fun(r$load_sum_d1, r$load_pta_d1, r$target_auc, r$load_df$DoseLabel, r$rec_load_idx, "AUC 0-24 h (mg\u00b7h/L)", "Loading dose \u2013 Day 1 (0\u201324h)")
```

```{r table_load_d1}
knitr::kable(data.frame("Dose"=r$load_df$DoseLabel, "AUC median"=round(r$load_sum_d1$median,1), "95% CI"=paste0(round(r$load_sum_d1$p2.5,1),"\u2013",round(r$load_sum_d1$p97.5,1)), "PTA \u226590%"=round(r$load_pta_d1,3), check.names=FALSE), align="c")
```

```{r multiday_load, results="asis"}
if(isTRUE(r$need_multiday) && !is.null(r$load_sum_d2)){
  cat("<h2>Loading Dose &mdash; AUC 24&ndash;48 h (Day 2, if repeated)</h2>")
}
```

```{r plot_load_d2, fig.height=4.5, fig.width=9}
if(isTRUE(r$need_multiday) && !is.null(r$load_sum_d2)){
  plot_fun(r$load_sum_d2, r$load_pta_d2, r$target_auc, r$load_df$DoseLabel, r$rec_load_idx, "AUC 24-48 h (mg\u00b7h/L)", "Loading dose \u2013 Day 2 (24\u201348h)")
}
```

```{r table_load_d2}
if(isTRUE(r$need_multiday) && !is.null(r$load_sum_d2)){
  knitr::kable(data.frame("Dose"=r$load_df$DoseLabel, "AUC median"=round(r$load_sum_d2$median,1), "95% CI"=paste0(round(r$load_sum_d2$p2.5,1),"\u2013",round(r$load_sum_d2$p97.5,1)), "PTA \u226590%"=round(r$load_pta_d2,3), check.names=FALSE), align="c")
}
```

```{r multiday_load3, results="asis"}
if(isTRUE(r$need_multiday) && !is.null(r$load_sum_d3)){
  cat("<h2>Loading Dose &mdash; AUC 48&ndash;72 h (Day 3, if repeated)</h2>")
}
```

```{r plot_load_d3, fig.height=4.5, fig.width=9}
if(isTRUE(r$need_multiday) && !is.null(r$load_sum_d3)){
  plot_fun(r$load_sum_d3, r$load_pta_d3, r$target_auc, r$load_df$DoseLabel, r$rec_load_idx, "AUC 48-72 h (mg\u00b7h/L)", "Loading dose \u2013 Day 3 (48\u201372h)")
}
```

```{r table_load_d3}
if(isTRUE(r$need_multiday) && !is.null(r$load_sum_d3)){
  knitr::kable(data.frame("Dose"=r$load_df$DoseLabel, "AUC median"=round(r$load_sum_d3$median,1), "95% CI"=paste0(round(r$load_sum_d3$p2.5,1),"\u2013",round(r$load_sum_d3$p97.5,1)), "PTA \u226590%"=round(r$load_pta_d3,3), check.names=FALSE), align="c")
}
```

<h2>Maintenance Dose &mdash; AUC 312&ndash;336 h (Day 14)</h2>

```{r plot_maint, fig.height=4.5, fig.width=9}
plot_fun(r$maint_sum, r$maint_pta, r$target_auc, r$maint_df$DoseLabel, r$rec_maint_idx, "AUC 312-336h (mg*h/L)", "Maintenance dose exposure (Day 14)")
```

```{r table_maint}
knitr::kable(data.frame("Dose"=r$maint_df$DoseLabel, "AUC median"=round(r$maint_sum$median,1), "95% CI"=paste0(round(r$maint_sum$p2.5,1),"\\u2013",round(r$maint_sum$p97.5,1)), "PTA"=round(r$maint_pta,3), check.names=FALSE), align="c")
```

</div>

</div>

', file = con)
  cat(report_footer(), "\n", file = con)
}

write_rmd_aposteriori <- function(path, logo_path = NULL) {
  con <- file(path, open = "wt", encoding = "UTF-8"); on.exit(close(con))

  cat('---
title: " "
output:
  html_document:
    self_contained: true
    toc: false
    df_print: kable
params:
  patient: NULL
  res: NULL
  apr: NULL
  expert_rec: NULL
  logo_path: NULL
  logo_html: ""
  mic_label: NULL
  target_auc: NULL
  clinical: NULL
  pharmacologue: NULL
  pk_model_ref: ""
  nsim: 10000
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
r <- params$res; p <- params$patient; apr <- params$apr; clin <- params$clinical
erec <- params$expert_rec
PK_MODEL_REF <- params$pk_model_ref
NSIM <- params$nsim
logo_html_str <- if(nzchar(params$logo_html)) params$logo_html else ""
today_txt <- format(Sys.Date(), "%d/%m/%Y")
mic_label <- params$mic_label; target_auc <- params$target_auc
sq <- function(x) {
  q <- quantile(x, probs=c(0.025,0.5,0.975), na.rm=TRUE)
  c(mean=mean(x,na.rm=TRUE), median=as.numeric(q[2]), p2.5=as.numeric(q[1]), p97.5=as.numeric(q[3]))
}
q_cl <- sq(r$posterior$CL); q_v <- sq(r$posterior$V)
COL_MAP <- "#1565C0"; COL_CM <- "#AD1457"; COL_MEDIAN <- "#37474F"
COL_BAND <- "#90CAF9"; COL_TARGET <- "#B71C1C"; COL_OPT_M <- "#1B5E20"
logo_b64 <- if (!is.null(params$logo_path) && file.exists(params$logo_path)) knitr::image_uri(params$logo_path) else ""
auc_mod <- function(t_start, t_end, dtimes, damounts, CL, V, tinf) {
  k <- CL/V; auc <- 0
  for (i in seq_along(dtimes)) {
    td <- dtimes[i]; D <- damounts[i]; R <- D/tinf
    ds <- t_start - td; de <- t_end - td
    if (de <= 0) next; ds <- max(ds, 0)
    if (ds < tinf) { t1 <- ds; t2 <- min(de, tinf); auc <- auc + (R/CL)*((t2-t1)-(1/k)*(exp(-k*t1)-exp(-k*t2))) }
    if (de > tinf) { t1 <- max(ds, tinf); t2 <- de; auc <- auc + (R/CL)*(1-exp(-k*tinf))*(1/k)*(exp(-k*(t1-tinf))-exp(-k*(t2-tinf))) }
  }
  auc
}
nmtran <- r$nmtran
hist_doses <- if(!is.null(nmtran)) nmtran[nmtran$evid==1,] else NULL
obs_rows <- if(!is.null(nmtran)) nmtran[nmtran$evid==0,] else NULL
OMEGA_MAT_rep <- matrix(c(0.0484,0,0,0.0841), 2, 2)
```

', file = con)
  cat(report_css(), "\n\n", file = con)
  cat(report_header_chunk(), "\n\n", file = con)

  cat('<div class="report-banner">
`r logo_html_str`
<div class="rpt-label">Clinical Pharmacology &mdash; Bayesian Individualisation Report</div>
<h1>Clinical Pharmacology Report</h1>
<div class="rpt-drug">Fluconazole &mdash; Individual Bayesian <em>a posteriori</em> estimation</div>
<div class="rpt-engine">&#9881; Engine: mrgsolve + saemix</div>
<div class="rpt-date">Report date: <strong>`r today_txt`</strong></div>
</div>
<div style="max-width:900px;margin:14px auto;padding:12px 16px;border:1px solid #B71C1C;border-left:5px solid #B71C1C;border-radius:6px;background:#FFF8F8;font-size:12px;color:#444;line-height:1.5;">
<b style="color:#B71C1C;">DISCLAIMER / AVERTISSEMENT</b><br>
This software is provided "AS IS", without warranty of any kind. It is a research prototype intended to support &mdash; never to replace &mdash; the clinical judgment of a qualified physician or pharmacologist. All dose recommendations must be reviewed, validated, and adjusted by a competent healthcare professional taking into account the full clinical context of the patient. The authors and their affiliated institutions decline all responsibility for any adverse outcome, treatment failure, injury, or death that may occur in connection with the use of this software or any decision based on its outputs.
</div>
<div class="page-body">

', file = con)
  cat(report_patient_block(), "\n\n", file = con)

  cat('<div class="section">

<h2>PK/PD Target</h2>
<div class="param-box">
<span class="lab">Target AUC:</span> `r if(!is.null(target_auc)) target_auc else "--"` mg&middot;h/L (AUC/MIC &gt; 100)*<br>
<span class="lab">PK model:</span> <em>`r PK_MODEL_REF`</em><br>

```{r cmi_text, results="asis"}
mic_exact_val <- clin$mic_exact
if (!is.null(mic_exact_val) && !is.na(mic_exact_val) && mic_exact_val > 0) {
  cat(paste0("<span class=\\"lab\\">Target MIC:</span> ", mic_exact_val, " mg/L (CMI exacte)<br>"))
} else {
  cat(paste0("<span class=\\"lab\\">Target MIC:</span> ", if(!is.null(mic_label)) mic_label else "--", "<br>"))
}
```

<span style="font-size:12px;color:#999;font-style:italic;">* Ashbee HR. et al. <em>J Antimicrob Chemother.</em> 2014. doi:10.1093/jac/dkt508.</span>
</div>

<h2>Individual A Posteriori PK Parameters</h2>
<div class="param-box">
<span class="lab">Population CL:</span> `r round(r$CLtyp,3)` L/h &nbsp;|&nbsp;
<span class="lab">MAP CL:</span> `r round(r$map_psi["CL"],3)` L/h<br>
<span class="lab">Population V:</span> `r round(r$Vtyp,1)` L &nbsp;|&nbsp;
<span class="lab">MAP V:</span> `r round(r$map_psi["V"],1)` L<br>
<span class="lab">Posterior CL [95% CI]:</span> `r round(q_cl["median"],3)` [`r round(q_cl["p2.5"],3)`&ndash;`r round(q_cl["p97.5"],3)`] L/h<br>
<span class="lab">Posterior Vd [95% CI]:</span> `r round(q_v["median"],1)` [`r round(q_v["p2.5"],1)`&ndash;`r round(q_v["p97.5"],1)`] L
</div>

<h2>Administered Dosage Regimen</h2>
<div class="param-box">
<span class="lab">Loading dose:</span> `r if(!is.null(apr$dose_first)) paste0(apr$dose_first, " mg") else "--"` le `r if(!is.null(apr$date_first)) apr$date_first else "--"`<br>
<span class="lab">Maintenance dose:</span> `r if(!is.null(apr$dose_maint)) paste0(apr$dose_maint, " mg/24h") else "--"`<br>
<span class="lab">Last administration:</span> `r if(!is.null(apr$date_last)) apr$date_last else "--"`<br>
<span class="lab">Infusion duration:</span> `r apr$tinf` h
</div>

<h2>Measured Concentrations</h2>

```{r obs_table}
if (!is.null(apr$obs_full) && nrow(apr$obs_full) > 0) {
  knitr::kable(data.frame("Date"=apr$obs_full$date, "Time"=sprintf("%.1fh",apr$obs_full$heure), "Concentration (mg/L)"=apr$obs_full$conc, "T post-dose (h)"=round(apr$obs_full$time,1), check.names=FALSE), align="c")
} else if (!is.null(r$obs)) {
  knitr::kable(data.frame("Time post-dose (h)"=r$obs$time, "Concentration (mg/L)"=r$obs$conc, check.names=FALSE), align="c")
}
```

', file = con)

  cat('
', file = con)

  # PK/PD Assessment chunk: written via writeLines to avoid quote-escape pitfalls
  # of the global cat('...') string (apostrophes inside HTML attributes).
  writeLines(c(
    "<h2>PK/PD Assessment per Sample</h2>",
    "",
    "```{r pkpd_assessment, results=\"asis\"}",
    "if (!is.null(hist_doses) && nrow(hist_doses) > 0 && !is.null(obs_rows) && nrow(obs_rows) > 0) {",
    "  d1_pp <- tryCatch(as.Date(apr$date_first, '%d/%m/%Y'), error=function(e) Sys.Date()-7)",
    "  tau_pp <- if(!is.null(apr$tau)) apr$tau else 24",
    "  tinf_pp <- r$tinf; target_pp <- if(!is.null(target_auc)) target_auc else 400",
    "  mic_exact_v <- clin$mic_exact",
    "  mic_is_doc <- !is.null(mic_exact_v) && !is.na(mic_exact_v) && is.finite(mic_exact_v) && mic_exact_v > 0",
    "  mic_is_glab <- !mic_is_doc && isTRUE(clin$glabrata)",
    "  mic_used_pp <- if (mic_is_doc) as.numeric(mic_exact_v) else if (mic_is_glab) 16 else 4",
    "  mic_src_lbl <- if (mic_is_doc) {",
    "    paste0('Documented MIC: <b>', mic_used_pp, ' mg/L</b>')",
    "  } else if (mic_is_glab) {",
    "    'MIC = <b>16 mg/L</b> &mdash; Clinical breakpoint <i>C. glabrata</i>'",
    "  } else {",
    "    'MIC = <b>4 mg/L</b> &mdash; Clinical breakpoint <i>C. albicans</i>, <i>C. parapsilosis</i>, <i>C. tropicalis</i>'",
    "  }",
    "  for (oo in seq_len(nrow(obs_rows))) {",
    "    ah <- obs_rows$time[oo]",
    "    doses_before <- hist_doses$time[hist_doses$time < ah]",
    "    if (length(doses_before) == 0) doses_before <- hist_doses$time[hist_doses$time <= ah]",
    "    t_start <- if(length(doses_before)>0) max(doses_before) else 0",
    "    doses_after_st <- hist_doses$time[hist_doses$time > t_start]",
    "    has_next_st <- length(doses_after_st) > 0 && min(doses_after_st) > ah",
    "    t_end_tau <- if(has_next_st) min(doses_after_st) else t_start + 24",
    "    t_post <- round(ah - t_start, 1)",
    "    auc_map <- auc_mod(t_start, t_end_tau, hist_doses$time, hist_doses$amt, r$map_psi['CL'], r$map_psi['V'], tinf_pp)",
    "    post <- r$posterior; n_s <- min(500, nrow(post)); idx_s <- sample(nrow(post), n_s)",
    "    auc_post <- sapply(idx_s, function(s) auc_mod(t_start, t_end_tau, hist_doses$time, hist_doses$amt, post$CL[s], post$V[s], tinf_pp))",
    "    qa <- quantile(auc_post, c(0.025,0.975), na.rm=TRUE)",
    "    days_since_start <- ah / 24",
    "    at_ss <- days_since_start >= 6",
    "    thr_target <- target_pp; thr_reduce <- 1.5*target_pp; thr_over <- 2.0*target_pp",
    "    tier <- if (auc_map < thr_target) 'not_achieved' else if (auc_map >= thr_over) 'overexposure' else if (auc_map >= thr_reduce && at_ss) 'reduce' else 'achieved'",
    "    bg  <- switch(tier, not_achieved='#FFF3E0', overexposure='#FDECEA', reduce='#E3F2FD', achieved='#E8F5E9')",
    "    brd <- switch(tier, not_achieved='#F57C00', overexposure='#C62828', reduce='#1565C0', achieved='#2E7D32')",
    "    icoTxt <- switch(tier, not_achieved='[!]', overexposure='[OVER]', reduce='[i]', achieved='[OK]')",
    "    auc_mic    <- round(auc_map / mic_used_pp, 0)",
    "    auc_mic_lo <- round(as.numeric(qa[1]) / mic_used_pp, 0)",
    "    auc_mic_hi <- round(as.numeric(qa[2]) / mic_used_pp, 0)",
    "    lbl <- switch(tier,",
    "      not_achieved = paste0(icoTxt, ' PK/PD target NOT ACHIEVED &mdash; Dosage adjustment required'),",
    "      overexposure = paste0(icoTxt, ' OVEREXPOSURE (AUC/MIC &ge; 200) &mdash; Dose reduction recommended'),",
    "      reduce       = paste0(icoTxt, ' Target achieved &mdash; Dose reduction possible (AUC/MIC &ge; 150, steady-state reached)'),",
    "      achieved     = paste0(icoTxt, ' PK/PD target ACHIEVED')",
    "    )",
    "    interval_lbl <- if(has_next_st) paste0('AUC 0&ndash;', round(t_end_tau-t_start,1), 'h [last dose &rarr; next dose]') else 'AUC 0&ndash;24h [last dose &rarr; last dose + 24h]'",
    "    box_html <- paste0(",
    "      \"<div style='background:\", bg, \";border-left:4px solid \", brd, \";border-radius:7px;padding:11px 13px;margin:9px 0;page-break-inside:avoid;box-shadow:0 1px 4px rgba(0,0,0,0.06);'>\",",
    "      \"<div style='font-size:11px;color:#555;font-weight:600;margin-bottom:7px;padding-bottom:5px;border-bottom:1px dashed rgba(0,0,0,0.12);'>Sample \", oo, \" &mdash; T post-dose: \", t_post, \"h &mdash; Observed concentration: \", round(obs_rows$DV[oo],2), \" mg/L</div>\",",
    "      \"<div style='background:#FFFFFF;border-radius:6px;padding:10px 12px;margin-bottom:7px;border:1.5px solid \", brd, \";box-shadow:inset 0 0 0 1px rgba(0,0,0,0.03);'>\",",
    "      \"<div style='font-size:11px;font-weight:800;color:\", brd, \";letter-spacing:1px;text-transform:uppercase;margin-bottom:4px;'>Primary PK/PD Marker</div>\",",
    "      \"<div style='font-size:24px;font-weight:900;color:\", brd, \";line-height:1;margin-bottom:4px;letter-spacing:-0.5px;'>AUC/MIC = \", auc_mic, \"</div>\",",
    "      \"<div style='font-size:12px;font-weight:700;color:#333;margin-bottom:3px;'>Target: AUC/MIC &ge; 100</div>\",",
    "      \"<div style='font-size:11px;color:#555;margin-bottom:5px;'>\", mic_src_lbl, \"</div>\",",
    "      \"<div style='font-size:13px;font-weight:800;color:#FFFFFF;background:\", brd, \";margin-top:5px;padding:5px 7px;border-radius:4px;text-align:center;'>\", lbl, \"</div>\",",
    "      \"</div>\",",
    "      \"<div style='font-size:11px;color:#444;margin-bottom:2px;'><b>AUC MAP = \", round(auc_map,1), \" mg&middot;h/L</b> &nbsp;|&nbsp; 95% CI: [\", round(qa[1],1), \"&ndash;\", round(qa[2],1), \"] mg&middot;h/L &nbsp;|&nbsp; AUC/MIC 95% CI: [\", auc_mic_lo, \"&ndash;\", auc_mic_hi, \"]</div>\",",
    "      \"<div style='font-size:10px;color:#666;font-style:italic;'>\", interval_lbl, \"</div>\",",
    "      \"</div>\"",
    "    )",
    "    cat(box_html, '\\n', sep='')",
    "  }",
    "}",
    "```",
    ""
  ), con)

  cat('
<h2>Individual Interval AUC (Analytical Integral)</h2>

```{r auc_interval}
if (!is.null(hist_doses) && nrow(hist_doses) > 0 && !is.null(obs_rows) && nrow(obs_rows) > 0) {
  d1 <- tryCatch(as.Date(apr$date_first, "%d/%m/%Y"), error=function(e) Sys.Date()-7)
  tau_h <- if(!is.null(apr$tau)) apr$tau else 24
  tinf_r <- r$tinf; target_v <- if(!is.null(target_auc)) target_auc else 400
  auc_rows <- lapply(seq_len(nrow(obs_rows)), function(oo) {
    ah <- obs_rows$time[oo]  # absolute time of observation
    # Last dose STRICTLY before observation (not coincident)
    doses_before <- hist_doses$time[hist_doses$time < ah]
    if (length(doses_before) == 0) doses_before <- hist_doses$time[hist_doses$time <= ah]
    t_start <- if(length(doses_before)>0) max(doses_before) else 0
    # AUC definition:
    #  - If a next dose exists AND obs before it: AUCtau [t_last → t_next] (complete interval)
    #  - If no next dose after obs: AUC [t_last → t_obs] (interval truncated at sample time)
    doses_after <- hist_doses$time[hist_doses$time > t_start]
    has_next_dose <- length(doses_after) > 0 && min(doses_after) > ah
    t_end <- if(has_next_dose) min(doses_after) else t_start + 24
    t_post <- round(ah - t_start, 2)
    date_admin  <- format(d1 + t_start/24, "%d/%m/%Y")
    date_obs_lbl <- format(d1 + ah/24, "%d/%m/%Y")
    # MAP AUC
    auc_mode <- auc_mod(t_start, t_end, hist_doses$time, hist_doses$amt, r$map_psi["CL"], r$map_psi["V"], tinf_r)
    # Posterior distribution
    post <- r$posterior; n_s <- min(500, nrow(post)); idx_s <- sample(nrow(post), n_s)
    auc_post <- sapply(idx_s, function(s) auc_mod(t_start, t_end, hist_doses$time, hist_doses$amt, post$CL[s], post$V[s], tinf_r))
    qa <- quantile(auc_post, c(0.025,0.5,0.975))
    statut <- if(auc_mode >= target_v) "\u2713 Achieved" else "\u2717 Not achieved"
    data.frame(
      "Last dose"=date_admin,
      "End of interval"=format(d1 + t_end/24, "%d/%m/%Y"),
      "T post-dose (h)"=t_post,
      "C obs (mg/L)"=obs_rows$DV[oo],
      "AUC MAP (mg\u00b7h/L)"=round(auc_mode,1),
      "95% CI"=paste0(round(as.numeric(qa[1]),1),"\u2013",round(as.numeric(qa[3]),1)),
      "Target"=target_v,
      "AUC ≥ target?"=statut,
      check.names=FALSE, stringsAsFactors=FALSE)
  })
  auc_df <- do.call(rbind, auc_rows)
  knitr::kable(auc_df, align="c",
    caption="AUCtau MAP — AUC[last dose &#8594; next dose] or AUC 0&#8211;24h if no next dose")
}
```

', file = con)

  cat('
<h2>Concentration-Time Profile</h2>

```{r profile_nmtran, fig.height=5, fig.width=10}
if (!is.null(hist_doses) && nrow(hist_doses) > 0) {
  tryCatch({
    CL_map <- r$map_psi["CL"]; V_map <- r$map_psi["V"]; tinf_r <- r$tinf; post <- r$posterior
    t_last <- max(hist_doses$time); obs_max_t <- if(nrow(obs_rows)>0) max(obs_rows$time) else t_last
    tau_h <- 24; t_end <- max(t_last + tau_h + 12, obs_max_t + 12)
    tgrid <- seq(0, t_end, length.out=500); k_map <- CL_map/V_map
    map_conc <- sapply(tgrid, function(tt) {
      c_tot <- 0
      for(i in seq_len(nrow(hist_doses))) { dt <- tt - hist_doses$time[i]; if(dt>0) c_tot <- c_tot + (hist_doses$amt[i]/tinf_r/CL_map)*ifelse(dt<=tinf_r,1-exp(-k_map*dt),(1-exp(-k_map*tinf_r))*exp(-k_map*(dt-tinf_r))) }
      max(c_tot, 1e-12)
    })
    n_sub <- min(150, nrow(post)); idx_sub <- sample(nrow(post), n_sub)
    post_mat <- matrix(NA, n_sub, length(tgrid))
    for(s in seq_len(n_sub)) {
      cl_s <- post$CL[idx_sub[s]]; v_s <- post$V[idx_sub[s]]; k_s <- cl_s/v_s
      for(j in seq_along(tgrid)) {
        c_tot <- 0
        for(i in seq_len(nrow(hist_doses))) { dt <- tgrid[j]-hist_doses$time[i]; if(dt>0) c_tot <- c_tot + (hist_doses$amt[i]/tinf_r/cl_s)*ifelse(dt<=tinf_r,1-exp(-k_s*dt),(1-exp(-k_s*tinf_r))*exp(-k_s*(dt-tinf_r))) }
        post_mat[s,j] <- max(c_tot, 1e-12)
      }
    }
    p50 <- apply(post_mat,2,median); p025 <- apply(post_mat,2,quantile,0.025); p975 <- apply(post_mat,2,quantile,0.975)
    has_ci <- any(is.finite(p50))
    obs_t <- if(nrow(obs_rows)>0) obs_rows$time/24 else c()
    obs_c <- if(nrow(obs_rows)>0) obs_rows$DV else c()
    tgrid_d <- tgrid/24
    all_vals <- c(if(has_ci) p975, obs_c, map_conc)
    ylim_max <- if(length(all_vals)>0) max(all_vals, na.rm=TRUE)*1.2 else 50
    d1_date <- tryCatch(as.Date(apr$date_first, "%d/%m/%Y"), error=function(e) Sys.Date()-7)
    day_ticks <- seq(0, max(tgrid_d), by=1); date_labels <- format(d1_date + day_ticks, "%d/%m")
    par(mar=c(5,4.5,3,10), bg="white", xpd=FALSE)
    plot(tgrid_d, map_conc, type="n", xlab="", ylab="Fluconazole concentration (mg/L)", main="A posteriori concentration-time profile", ylim=c(0,ylim_max), las=1, xaxt="n")
    axis(1, at=day_ticks, labels=date_labels, cex.axis=0.8); mtext("Date", side=1, line=2.5, cex=0.85)
    grid(col="#ECEFF1", lty=1)
    if(has_ci) polygon(c(tgrid_d,rev(tgrid_d)), c(p025,rev(p975)), col=adjustcolor(COL_BAND,.35), border=NA)
    lines(tgrid_d, map_conc, lwd=3, col=COL_MAP)
    if(length(obs_t)>0) points(obs_t, obs_c, pch=21, cex=2.2, bg="#FF6F00", col="white", lwd=2)
    for(i in seq_len(nrow(hist_doses))) {
      x_d <- hist_doses$time[i]/24
      text(x_d, ylim_max*0.97, paste0(hist_doses$amt[i],"mg"), cex=0.55, col="#00539F", srt=90, adj=c(1,0.5))
      segments(x_d, 0, x_d, ylim_max*0.06, col="#00539F", lwd=0.8, lty=3)
    }
    par(xpd=TRUE)
    legend("right", inset=c(-0.18,0), bty="o", box.col="#E0E0E0", bg="white", cex=0.78,
      legend=c("MAP curve","95% CI","Observations"),
      lty=c(1,NA,NA), lwd=c(3,NA,NA), pch=c(NA,15,21), pt.cex=c(NA,2,1.8),
      pt.bg=c(NA,adjustcolor(COL_BAND,.5),"#FF6F00"),
      col=c(COL_MAP,adjustcolor(COL_BAND,.6),"#FF6F00"))
    par(xpd=FALSE)
  }, error=function(e) { plot.new(); text(0.5,0.5,paste("Profile error:",e$message)) })
}
```

<h2>PK Parameter Distributions (prior vs posterior)</h2>

```{r distrib_cl_v, fig.height=3.5, fig.width=10}
tryCatch({
  L <- t(chol(OMEGA_MAT_rep)); z <- matrix(rnorm(5000*2), ncol=2); eta_p <- z %*% t(L)
  par(mfrow=c(1,2), mar=c(4.5,4.5,3,1.5), bg="white")
  pv_cl <- r$CLtyp*exp(eta_p[,1]); dv_cl <- r$posterior$CL; mv_cl <- r$map_psi["CL"]
  d1c <- density(pv_cl, na.rm=TRUE); d2c <- density(dv_cl, na.rm=TRUE)
  ylm <- max(d1c$y, d2c$y)*1.18
  plot(d1c, main="CL distribution (prior vs posterior)", xlab="CL (L/h)", lwd=2, col="grey70", ylim=c(0,ylm), font.main=2, las=1); grid(col="#ECEFF1",lty=1)
  polygon(c(d2c$x,rev(d2c$x)), c(d2c$y,rep(0,length(d2c$y))), col=adjustcolor(COL_CM,.22), border=NA)
  lines(d2c, lwd=2.5, col=COL_CM); abline(v=mv_cl, col=COL_MAP, lty=1, lwd=3)
  legend("topright", bty="n", legend=c("A priori","A posteriori","MAP"), col=c("grey70",COL_CM,COL_MAP), lty=c(1,1,1), lwd=c(2,2.5,3))
  pv_v <- r$Vtyp*exp(eta_p[,2]); dv_v <- r$posterior$V; mv_v <- r$map_psi["V"]
  d1v <- density(pv_v, na.rm=TRUE); d2v <- density(dv_v, na.rm=TRUE)
  ylm2 <- max(d1v$y, d2v$y)*1.18
  plot(d1v, main="Distribution of Vd values (prior vs posterior)", xlab="V (L)", lwd=2, col="grey70", ylim=c(0,ylm2), font.main=2, las=1); grid(col="#ECEFF1",lty=1)
  polygon(c(d2v$x,rev(d2v$x)), c(d2v$y,rep(0,length(d2v$y))), col=adjustcolor(COL_CM,.22), border=NA)
  lines(d2v, lwd=2.5, col=COL_CM); abline(v=mv_v, col=COL_MAP, lty=1, lwd=3)
  legend("topright", bty="n", legend=c("A priori","A posteriori","MAP"), col=c("grey70",COL_CM,COL_MAP), lty=c(1,1,1), lwd=c(2,2.5,3))
}, error=function(e) { plot.new(); text(0.5,0.5,paste("Distribution error:",e$message)) })
```

```{r pk_table}
knitr::kable(data.frame(
  "Parameter"=c("CL (L/h)","V (L)"),
  "Mean"=round(c(q_cl["mean"],q_v["mean"]),3),
  "Median"=round(c(q_cl["median"],q_v["median"]),3),
  "P2.5"=round(c(q_cl["p2.5"],q_v["p2.5"]),3),
  "P97.5"=round(c(q_cl["p97.5"],q_v["p97.5"]),3),
  check.names=FALSE), align="c")
```

', file = con)

  # Adaptation posologique
  cat('

<h2>Dosage Adjustment</h2>

```{r adaptation, results="asis"}
if (!is.null(erec) && !is.null(erec$schedule)) {
  sched <- erec$schedule; dose_ss <- erec$dose_ss; tau <- erec$tau; tgt <- erec$target
  sched_txt <- paste(sapply(seq_along(sched), function(i) {
    if(sched[i]==0) paste0("Day ",i," : PAUSE") else paste0("Day ",i," : ",sched[i]," mg")
  }), collapse=" &rarr; ")
  has_pause <- any(sched==0)
  cat(paste0("<div class=\\"rec-ok\\"><b>Optimised dosage regimen:</b><br>",
    sched_txt,
    "<br>Maintenance dose: <b>", dose_ss, " mg / ", tau, "h</b>",
    if(has_pause) "<br><em>&#9888; Pause th&eacute;rapeutique propos&eacute;e pour &eacute;viter le surdosage.</em>" else "",
    "</div>"))
  if (!is.null(hist_doses) && nrow(hist_doses) > 0) {
    t_last_hist <- max(hist_doses$time); t_new_start <- t_last_hist + tau
    n_days <- length(sched); n_new <- n_days * ceiling(24/tau)
    new_times <- t_new_start + (seq_len(n_new)-1) * tau
    new_amts <- sapply(new_times, function(t) {
      day_idx <- max(1, min(floor((t-t_new_start)/24)+1, length(sched))); sched[day_idx]
    })
    all_dt <- c(hist_doses$time, new_times); all_da <- c(hist_doses$amt, new_amts)
    tinf_r <- r$tinf
    d1_date <- tryCatch(as.Date(apr$date_first, "%d/%m/%Y"), error=function(e) Sys.Date()-7)
    cat("<table class=\\"table table-sm\\"><thead><tr><th>Day</th><th>Dates</th><th>Dose (mg)</th><th>AUC (mg&middot;h/L)</th><th>Target</th><th>Status</th></tr></thead><tbody>")
    for(d in seq_len(n_days)) {
      ts <- t_new_start + (d-1)*24; te <- t_new_start + d*24
      auc_d <- auc_mod(ts, te, all_dt, all_da, r$map_psi["CL"], r$map_psi["V"], tinf_r)
      dd <- if(d <= length(sched)) sched[d] else dose_ss
      statut <- if(auc_d >= tgt) "&#10003;" else "&#10007;"
      cat(paste0("<tr><td>J",d,"</td><td>",format(d1_date+ts/24,"%d/%m")," &ndash; ",format(d1_date+te/24,"%d/%m"),"</td><td>",if(dd==0)"PAUSE"else dd,"</td><td>",round(auc_d,1),"</td><td>",tgt,"</td><td>",statut,"</td></tr>"))
    }
    cat("</tbody></table>")
  }
} else {
  cat("<em>Dosage adjustment not computed. Use the Find optimal regimen button in the Dosage adjustment tab.</em>")
}
```

```{r adapt_plot, fig.height=5.5, fig.width=10}
if (!is.null(erec) && !is.null(erec$schedule) && !is.null(hist_doses) && nrow(hist_doses) > 0) {
  tryCatch({
    sched <- erec$schedule; tau <- erec$tau; tinf_r <- r$tinf
    CL_map <- r$map_psi["CL"]; V_map <- r$map_psi["V"]; post <- r$posterior
    t_last_hist <- max(hist_doses$time)
    t_new_start <- t_last_hist + tau
    n_days <- length(sched); n_new <- n_days * ceiling(24/tau)
    new_times <- t_new_start + (seq_len(n_new)-1) * tau
    new_amts <- sapply(new_times, function(t) {
      day_idx <- max(1, min(floor((t-t_new_start)/24)+1, length(sched))); sched[day_idx]
    })
    all_times <- c(hist_doses$time, new_times)
    all_amts <- c(hist_doses$amt, new_amts)
    t_end <- max(new_times) + tau
    tgrid <- seq(0, t_end, length.out=600)
    # MAP curve
    k_map <- CL_map/V_map
    map_conc <- sapply(tgrid, function(tt) {
      c_tot <- 0
      for(i in seq_along(all_times)) { dt <- tt - all_times[i]; if(dt>0) c_tot <- c_tot + (all_amts[i]/tinf_r/CL_map)*ifelse(dt<=tinf_r,1-exp(-k_map*dt),(1-exp(-k_map*tinf_r))*exp(-k_map*(dt-tinf_r))) }
      max(c_tot, 1e-12)
    })
    # Posterior CI
    n_sub <- min(150, nrow(post)); idx_sub <- sample(nrow(post), n_sub)
    post_mat <- matrix(NA, n_sub, length(tgrid))
    for(s in seq_len(n_sub)) {
      cl_s <- post$CL[idx_sub[s]]; v_s <- post$V[idx_sub[s]]; k_s <- cl_s/v_s
      for(j in seq_along(tgrid)) {
        c_tot <- 0
        for(i in seq_along(all_times)) { dt <- tgrid[j]-all_times[i]; if(dt>0) c_tot <- c_tot + (all_amts[i]/tinf_r/cl_s)*ifelse(dt<=tinf_r,1-exp(-k_s*dt),(1-exp(-k_s*tinf_r))*exp(-k_s*(dt-tinf_r))) }
        post_mat[s,j] <- max(c_tot, 1e-12)
      }
    }
    p50 <- apply(post_mat,2,median); p025 <- apply(post_mat,2,quantile,0.025); p975 <- apply(post_mat,2,quantile,0.975)
    tgrid_d <- tgrid/24
    obs_t <- if(nrow(obs_rows)>0) obs_rows$time/24 else c()
    obs_c <- if(nrow(obs_rows)>0) obs_rows$DV else c()
    ylim_max <- max(c(p975, obs_c, map_conc), na.rm=TRUE)*1.15
    d1_date <- tryCatch(as.Date(apr$date_first, "%d/%m/%Y"), error=function(e) Sys.Date()-7)
    day_ticks <- seq(0, max(tgrid_d), by=1); date_labels <- format(d1_date + day_ticks, "%d/%m")
    par(mar=c(5,4.5,3,10), bg="white", xpd=FALSE)
    plot(tgrid_d, map_conc, type="n", xlab="", ylab="Fluconazole concentration (mg/L)", main="Dosage adjustment simulation", ylim=c(0,ylim_max), las=1, xaxt="n")
    axis(1, at=day_ticks, labels=date_labels, cex.axis=0.8); mtext("Date", side=1, line=2.5, cex=0.85)
    grid(col="#ECEFF1", lty=1)
    rect(t_new_start/24, 0, t_end/24, ylim_max, col=adjustcolor("#E8F5E9",0.25), border=NA)
    polygon(c(tgrid_d,rev(tgrid_d)), c(p025,rev(p975)), col=adjustcolor(COL_BAND,.35), border=NA)
    lines(tgrid_d, map_conc, lwd=3, col=COL_MAP)
    abline(v=t_new_start/24, col="#2E7D32", lty=2, lwd=2)
    if(length(obs_t)>0) points(obs_t, obs_c, pch=21, cex=1.8, bg="#FF6F00", col="white", lwd=1.5)
    for(i in seq_along(all_times)) {
      x_d <- all_times[i]/24
      col_d <- if(all_times[i] < t_new_start) "#00539F" else "#2E7D32"
      text(x_d, ylim_max*0.97, paste0(all_amts[i],"mg"), cex=0.5, col=col_d, srt=90, adj=c(1,0.5))
      segments(x_d, 0, x_d, ylim_max*0.06, col=col_d, lwd=0.8, lty=3)
    }
    par(xpd=TRUE)
    legend("right", inset=c(-0.18,0), bty="o", box.col="#E0E0E0", bg="white", cex=0.78,
      legend=c("MAP curve","95% CI","Observations","New regimen"),
      lty=c(1,NA,NA,2), lwd=c(3,NA,NA,2), pch=c(NA,15,21,NA), pt.cex=c(NA,2,1.4,NA),
      pt.bg=c(NA,adjustcolor(COL_BAND,.5),"#FF6F00",NA),
      col=c(COL_MAP,adjustcolor(COL_BAND,.6),"#FF6F00","#2E7D32"))
    par(xpd=FALSE)
  }, error=function(e) { plot.new(); text(0.5,0.5,paste("Adjustment plot error:",e$message)) })
}
```

</div>

</div>

', file = con)
  cat(report_footer(), "\n", file = con)
}


# ==============================================================================
# PATIENT MODAL
# ==============================================================================
patient_modal_ui <- function(patient, title="Patient identification") {
  modalDialog(title=title, size="l",
    fluidRow(column(6, textInput("pat_last_name","Last name",value=patient$last_name)), column(6, textInput("pat_first_name","First name",value=patient$first_name))),
    fluidRow(column(4, dateInput("pat_dob","Date of birth",value=patient$dob,format="yyyy-mm-dd")), column(4, textInput("pat_ipp","IPP",value=patient$ipp)), column(4, radioButtons("pat_sexe","Sex",choices=c("Male","Female"),selected=patient$sexe,inline=TRUE))),
    fluidRow(column(4, textInput("pat_service","Department",value=patient$service)), column(4, textInput("pat_uf","Unity number",value=patient$uf)), column(4, textInput("pat_hopital","Hospital",value=patient$hopital))),
    textInput("pat_adresse","Address (optional)",value=patient$adresse),
    footer=tagList(modalButton("Cancel"), actionButton("confirm_patient","Confirm",class="btn-success")), easyClose=FALSE)
}

# ==============================================================================
# DISCLAIMER UI HELPER
# ==============================================================================
# Reusable disclaimer box, displayed at the top of the A priori and A posteriori
# Shiny tabs as a constant reminder to the clinician.
disclaimer_box <- function() {
  div(
    style = paste0(
      "max-width:100%;margin:0 0 14px 0;padding:10px 14px;",
      "border:1px solid #B71C1C;border-left:5px solid #B71C1C;",
      "border-radius:6px;background:#FFF8F8;",
      "font-size:12px;color:#444;line-height:1.5;"
    ),
    HTML(paste0(
      "<b style='color:#B71C1C;'>DISCLAIMER / AVERTISSEMENT</b><br>",
      "This software is provided \"AS IS\", without warranty of any kind. ",
      "It is a research prototype intended to support &mdash; never to replace &mdash; ",
      "the clinical judgment of a qualified physician or pharmacologist. ",
      "All dose recommendations must be reviewed, validated, and adjusted by a ",
      "competent healthcare professional taking into account the full clinical ",
      "context of the patient. The authors and their affiliated institutions ",
      "decline all responsibility for any adverse outcome, treatment failure, ",
      "injury, or death that may occur in connection with the use of this ",
      "software or any decision based on its outputs."
    ))
  )
}

# ==============================================================================
# UI
# ==============================================================================
ui <- fluidPage(
  theme = bs_theme(version=5, bootswatch="flatly"),
  tags$head(
    tags$link(rel="preconnect", href="https://fonts.googleapis.com"),
    tags$link(rel="preconnect", href="https://fonts.gstatic.com", crossorigin=""),
    tags$link(href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800;900&display=swap", rel="stylesheet"),
    tags$meta(name="viewport", content="width=device-width, initial-scale=1, maximum-scale=1")
  ),
  tags$style(HTML("
    /* ---------- Modern foundation ---------- */
    html, body { font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: #1a1d29; background: #f4f6fa; -webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; }
    .container-fluid { max-width: 1600px; margin: 0 auto; padding: 0 20px; }

    /* ---------- Cards ---------- */
    .card { border: 1px solid #e3e8f0 !important; border-radius: 12px !important; box-shadow: 0 1px 3px rgba(15,23,42,0.04), 0 1px 2px rgba(15,23,42,0.03); transition: box-shadow 0.2s ease, transform 0.2s ease; margin-bottom: 16px; background: white; }
    .card:hover { box-shadow: 0 4px 12px rgba(15,23,42,0.06), 0 2px 4px rgba(15,23,42,0.04); }
    .card-header { font-weight: 700 !important; font-size: 14px !important; letter-spacing: 0.2px; text-transform: uppercase; color: #475569; background: #f8fafc !important; border-bottom: 1px solid #e3e8f0 !important; padding: 12px 18px !important; border-radius: 12px 12px 0 0 !important; }
    .card-body { padding: 16px 18px !important; }

    /* ---------- Forms ---------- */
    .compact-form .form-group, .compact-form .shiny-input-container { margin-bottom: 8px !important; }
    .compact-form label, .compact-form .form-label { font-size: 12.5px !important; margin-bottom: 3px !important; }
    .compact-form .form-control, .compact-form .form-select, .compact-form input[type='number'] { font-size: 13.5px !important; padding: 6px 10px !important; height: 34px !important; }
    .compact-form .selectize-input { padding: 5px 10px !important; min-height: 34px !important; font-size: 13.5px !important; }
    .compact-form .checkbox label { font-size: 13px !important; line-height: 1.4 !important; white-space: nowrap !important; }
    .compact-form .form-check-label { font-size: 12.5px !important; white-space: nowrap !important; }
    .selectize-input { background: white !important; }
    .selectize-input > input { background: transparent !important; box-shadow: none !important; opacity: 0 !important; position: absolute !important; left: -9999px !important; width: 1px !important; }
    .selectize-input.has-items > input { display: none !important; }
    .selectize-control.single .selectize-input::after { border-color: #64748b transparent transparent transparent !important; }
    label, .form-label { font-weight: 600 !important; font-size: 13px !important; color: #334155; margin-bottom: 5px !important; }
    .form-control, .form-select, input[type='number'], input[type='text'] { border: 1.5px solid #e2e8f0 !important; border-radius: 8px !important; font-size: 14px !important; padding: 8px 12px !important; transition: all 0.15s ease; }
    .form-control:focus, .form-select:focus, input[type='number']:focus { border-color: #00539B !important; box-shadow: 0 0 0 3px rgba(0,83,155,0.12) !important; outline: none !important; }

    /* ---------- Tabs ---------- */
    .nav-tabs { border-bottom: 2px solid #e3e8f0 !important; margin-bottom: 16px; }
    .nav-tabs .nav-link { border: none !important; border-bottom: 3px solid transparent !important; color: #64748b !important; font-weight: 600 !important; font-size: 14px !important; padding: 10px 20px !important; transition: all 0.15s ease; background: transparent !important; }
    .nav-tabs .nav-link:hover { color: #00539B !important; border-bottom-color: #cbd5e1 !important; }
    .nav-tabs .nav-link.active { color: #00539B !important; background: transparent !important; border-bottom-color: #00539B !important; }

    /* ---------- Buttons ---------- */
    .btn { font-weight: 600 !important; border-radius: 8px !important; transition: all 0.15s ease; letter-spacing: 0.2px; }
    .btn-primary { background: linear-gradient(135deg, #00539B 0%, #0073cc 100%) !important; border: none !important; box-shadow: 0 1px 2px rgba(0,83,155,0.2); }
    .btn-primary:hover { transform: translateY(-1px); box-shadow: 0 4px 12px rgba(0,83,155,0.3) !important; filter: brightness(1.05); }
    .btn-primary:active { transform: translateY(0); }
    .btn-outline-danger:hover { transform: translateY(-1px); }

    /* ---------- Tables ---------- */
    table { width: 100% !important; border-collapse: separate !important; border-spacing: 0; font-size: 13px; }
    th, td { white-space: nowrap; text-align: center; vertical-align: middle; padding: 9px 12px !important; }
    th { font-weight: 700 !important; background: #f1f5f9 !important; color: #475569; font-size: 12px; letter-spacing: 0.3px; text-transform: uppercase; border-bottom: 2px solid #cbd5e1 !important; }
    tbody tr { transition: background 0.1s ease; }
    tbody tr:hover { background: #f8fafc; }
    tbody td { border-bottom: 1px solid #f1f5f9 !important; }
    tbody tr:last-child td { border-bottom: none !important; }

    /* ---------- Existing FlucoShiny rules (preserved) ---------- */
    #main_tabs > .nav { display: none !important; }
    .card-header { font-weight: 700 !important; }
    ul.refdash { list-style: none; padding-left: 0; margin-bottom: 0; }
    ul.refdash li { margin-bottom: 8px; padding-left: 14px; position: relative; line-height: 1.5; }
    ul.refdash li::before { content: ''; position: absolute; left: 0; top: 8px; width: 5px; height: 5px; background: #00539B; border-radius: 50%; }
    .extrapWarn { color: #C62828; font-weight: 800; }
    .tightGap { margin-top: -8px; }
    .tightGapBtn { margin-top: -10px; margin-bottom: -6px; }
    .btnWide { width: 100% !important; }

    /* ---------- Recommendation boxes ---------- */
    .rec-box { border-radius: 10px; padding: 16px 20px; margin-bottom: 14px; border-left: 5px solid; font-size: 14px; box-shadow: 0 1px 3px rgba(15,23,42,0.05); }
    .rec-above { background: linear-gradient(135deg, #f0fdf4 0%, #e8f5e9 100%); border-color: #2E7D32; }
    .rec-below { background: linear-gradient(135deg, #fef2f2 0%, #fdecea 100%); border-color: #C62828; }
    .rec-warn  { background: linear-gradient(135deg, #fffbeb 0%, #fff3e0 100%); border-color: #F57C00; }
    .rec-box .rec-status { font-size: 13px; font-weight: 800; text-transform: uppercase; letter-spacing: 0.4px; margin-bottom: 6px; }
    .rec-box .rec-dose { font-size: 22px; font-weight: 800; margin: 8px 0; display: block; line-height: 1.1; }
    .rec-box .rec-dose.red { color: #C62828; }
    .rec-box .rec-dose.green { color: #2E7D32; }
    .rec-box .rec-dose.orange { color: #E65100; }
    .rec-box .rec-detail { font-size: 13px; color: #475569; line-height: 1.5; }

    /* ---------- Patient badge ---------- */
    .pat-badge { background: linear-gradient(135deg, #eff6ff 0%, #E3F2FD 100%); border: 1px solid #93c5fd; border-radius: 8px; padding: 8px 14px; font-size: 13px; margin-bottom: 10px; box-shadow: 0 1px 2px rgba(59,130,246,0.08); }

    /* ---------- Modal & notifications ---------- */
    .modal-dialog { display: flex !important; align-items: center !important; min-height: calc(100vh - 60px) !important; margin: 30px auto !important; }
    .modal-dialog.modal-lg { max-width: 800px !important; }
    .modal-content { border-radius: 12px !important; border: none !important; box-shadow: 0 20px 60px rgba(15,23,42,0.25); width: 100% !important; }
    .modal-header { border-bottom: 1px solid #e3e8f0 !important; padding: 16px 20px !important; }
    .modal-title { font-weight: 700 !important; }
    .shiny-notification { border-radius: 10px !important; box-shadow: 0 4px 12px rgba(15,23,42,0.15) !important; }

    /* ---------- Subtle scrollbar ---------- */
    ::-webkit-scrollbar { width: 10px; height: 10px; }
    ::-webkit-scrollbar-track { background: #f1f5f9; }
    ::-webkit-scrollbar-thumb { background: #cbd5e1; border-radius: 6px; }
    ::-webkit-scrollbar-thumb:hover { background: #94a3b8; }
  ")),
  div(style="background:linear-gradient(135deg,#003a75 0%,#00539B 60%,#0073cc 100%);color:white;padding:22px 28px;margin-bottom:0;box-shadow:0 4px 12px rgba(0,83,155,0.15);position:relative;overflow:hidden;",
    div(style="position:absolute;top:0;left:0;right:0;height:3px;background:linear-gradient(90deg,#FFD700 0%,#FF8C00 30%,#FF6F00 60%,#FFD700 100%);"),
    fluidRow(
      column(12, div(
        h2(style="margin:0 0 4px 0;font-size:26px;font-weight:900;letter-spacing:-0.5px;","FlucoShiny"),
        p(style="margin:0 0 2px 0;font-size:14px;opacity:.92;line-height:1.5;",
          HTML("Fluconazole dosage individualisation tool &mdash; <em>a priori</em> &amp; Bayesian <em>a posteriori</em> &mdash; Adult ICU")),
        p(style="margin:0;font-size:13px;opacity:.80;line-height:1.4;",
          HTML(paste0("Based on the PK model of Bendjilali-Sabiani et al., digestive ICU, 2026 — n = ",
                      format(NSIM, big.mark="’")," simulations")))
      ))
    )
  ),
  div(style="height:4px;background:linear-gradient(to right,#FFD700,#00539B);margin-bottom:18px;"),
  fluidRow(
    column(4,
      card(card_header(HTML("\U0001F464 \u2002Patient characteristics")),
        div(class="compact-form",
          fluidRow(
            column(6, numericInput("wt","Body weight (kg)",70,40,200,step=5)),
            column(6, numericInput("ckd",HTML("CKD-EPI (mL/min/1.73m\u00b2)"),90,0,250,step=5))
          ),
          fluidRow(
            column(6, selectInput("cirrh","Cirrhosis",choices=c("No","Yes"),selected="No")),
            column(6, numericInput("tinf","Infusion (h)",value=0.5,min=0.05,max=24,step=0.05))
          ),
          fluidRow(
            column(12, numericInput("mic_exact","Exact MIC if available (mg/L)",value=NA,min=0))
          ),
          div(style="margin-top:4px;padding-top:8px;border-top:1px solid #e3e8f0;",
            div(style="white-space:nowrap;overflow:hidden;font-size:13px;",
              checkboxInput("glabrata",HTML("Suspected/proven <i>Candida glabrata</i> infection"))
            )
          )
        )
      ), br(),
      card(card_header(HTML("\U0001F3AF \u2002PK/PD target")), div(style="font-size:13.5px;line-height:1.55;color:#334155;",uiOutput("target_text"))), br(),
      card(card_header(HTML("\U0001F4DA \u2002References")),
        div(style="font-size:12.5px;line-height:1.55;color:#475569;",
          tags$ul(class="refdash",style="margin-bottom:0;",
            tags$li(tags$b("Ashbee HR. et al."), " Therapeutic drug monitoring (TDM) of antifungal agents: guidelines from the British Society for Medical Mycology. ",
                    tags$i("J Antimicrob Chemother."), " 2014. ",
                    tags$span(style="color:#94a3b8;font-size:11.5px;","doi:10.1093/jac/dkt508")),
            tags$li(tags$b("Gómez-López A."), " Antifungal therapeutic drug monitoring: focus on drugs without a clear recommendation. ",
                    tags$i("Clin Microbiol Infect."), " 2020. ",
                    tags$span(style="color:#94a3b8;font-size:11.5px;","doi:10.1016/j.cmi.2020.05.037")),
            tags$li(tags$b("EUCAST"), " v6.0, 2025. ",
                    tags$a(href="http://www.eucast.org",target="_blank",style="color:#00539B;text-decoration:none;","www.eucast.org"))
          )
        )
      )
    ),
    column(8, tabsetPanel(id="main_tabs",
      tabPanel("Home", br(),
        div(style="text-align:center;padding:30px 20px;",
          tags$img(src="flucoshiny.png",style="height:270px;max-width:780px;object-fit:contain;margin-bottom:16px;",onerror="this.style.display='none'"),
          h3(style="color:#00539F;font-weight:900;font-size:28px;margin-bottom:8px;","FlucoShiny"),
          p(style="color:#666;font-size:15px;margin-bottom:30px;","Dosage decision support - Fluconazole IV, adult ICU population"),
          fluidRow(
            column(6,
              actionButton("goto_apriori","A PRIORI",class="btn-primary",style="width:100%;height:120px;font-size:22px;font-weight:800;border-radius:12px;background:#00539F;border:none;",icon=icon("calculator")),
              p(style="margin-top:10px;color:#888;font-size:13px;",paste0("A priori simulation — n = ",format(NSIM,big.mark=" ")," IIV draws"))
            ),
            column(6,
              actionButton("goto_aposteriori","A POSTERIORI",class="btn-success",style="width:100%;height:120px;font-size:22px;font-weight:800;border-radius:12px;background:#2E7D32;border:none;",icon=icon("user-md")),
              p(style="margin-top:10px;color:#888;font-size:13px;","Individual Bayesian estimation (saemix conditional posterior)")
            )
          )
        )
      ),
      tabPanel("A priori", br(),
        fluidRow(
          column(6, actionButton("back_home_ap","Back",class="btn-outline-secondary btn-sm",icon=icon("arrow-left"),style="margin-bottom:10px;")),
          column(6, actionButton("reset_ap","Reset",class="btn-outline-danger btn-sm",icon=icon("refresh"),style="margin-bottom:10px;float:right;"))
        ),
        actionButton("run","Run a priori PK/PD evaluation",class="btn-primary",style="width:100%;height:84px;font-size:19px;font-weight:800;border-radius:12px;background:#00539F;border:none;margin-bottom:12px;"),
        uiOutput("ap_stale_warning"),
        uiOutput("card_pk"), uiOutput("alert_global"), uiOutput("box_loading"), uiOutput("box_maintenance"), br(),
        tabsetPanel(id="apriori_subtabs",
          tabPanel("Plots", br(), uiOutput("card_results")),
          tabPanel("PTA by dose", br(),
            card(card_header("Loading dose \u2013 PTA vs dose (AUC 0-24h)"), plotOutput("pta_plot_load",height=350)),
            br(),
            card(card_header("Maintenance dose \u2013 PTA vs dose (AUC steady-state)"), plotOutput("pta_plot_maint",height=350))),
          tabPanel("Repeated dose", br(),
            uiOutput("repeated_dose_tab_content"))
        ), br(), uiOutput("report_ui_apriori"), br(), disclaimer_box()),
      tabPanel("A posteriori", br(),
        fluidRow(
          column(6, actionButton("back_home_post","Back",class="btn-outline-secondary btn-sm",icon=icon("arrow-left"),style="margin-bottom:10px;")),
          column(6, actionButton("reset_post","Reset",class="btn-outline-danger btn-sm",icon=icon("refresh"),style="margin-bottom:10px;float:right;"))
        ),
        uiOutput("patient_badge_ap"),
        card(card_header("A posteriori Bayesian estimation"),
          tags$b("Administrations"),
          fluidRow(
            column(4, dateInput("ap_date_first","Date of 1st administration",value=Sys.Date()-7,format="yyyy-mm-dd")),
            column(4, numericInput("ap_heure_first","Administration time",value=8,min=0,max=23.99,step=0.5)),
            column(4, numericInput("ap_dose_first","Loading dose (mg)",value=800,min=0,step=25))
          ),
          fluidRow(
            column(4, dateInput("ap_date_maint","Date of 1st maintenance dose",value=Sys.Date()-6,format="yyyy-mm-dd")),
            column(4, numericInput("ap_heure_maint","Administration time",value=8,min=0,max=23.99,step=0.5)),
            column(4, numericInput("ap_dose_maint","Maintenance dose (mg)",value=400,min=0,step=25))
          ),
          checkboxInput("ap_regular","Regular maintenance dose every 24h",value=TRUE),
          conditionalPanel(condition="input.ap_regular==true",
            fluidRow(
              column(4, dateInput("ap_date_last","Last maintenance dose administration",value=Sys.Date()-1,format="yyyy-mm-dd"))
            )
          ),
          conditionalPanel(condition="input.ap_regular==false",
            numericInput("ap_n_doses","Number of administrations",value=3,min=1,max=30,step=1),
            uiOutput("ap_dose_inputs")
          ),
          tags$hr(),
          tags$b("Plasma drug levels"),
          numericInput("ap_n_conc","Number of fluconazole measurements",value=1,min=1,max=10,step=1),
          uiOutput("ap_conc_inputs"),
          uiOutput("ap_ss_warning")
        ), br(),
        card(card_header("Import Monolix file (optional)"),
          fluidRow(
            column(6, fileInput("monolix_file","Excel file (.xlsx)",accept=c(".xlsx",".xls"))),
            column(6, uiOutput("monolix_patient_select"))
          ),
          p(style="font-size:12px;color:#888;","Format : colonnes ID, DateTime, Time, AMT, Infusion duration, OBS, BW, CKD, cirrhosis"),
          uiOutput("monolix_status")
        ), br(),
        actionButton("run_ap","Run a posteriori Bayesian estimation",class="btn-success",style="width:100%;height:84px;font-size:19px;font-weight:800;border-radius:12px;background:#2E7D32;border:none;"),
        br(), br(), uiOutput("ap_summary"), br(),
        tabsetPanel(id="ap_subtabs",
          tabPanel("Profil & distributions", br(), card(plotOutput("ap_plot_profile",height=420)), br(),
            fluidRow(column(6,card(card_header("CL distribution"),plotOutput("ap_plot_cl",height=280))), column(6,card(card_header("Distribution of Vd values"),plotOutput("ap_plot_v",height=280)))), br(),
            card(card_header("PK parameters"),tableOutput("ap_param_table")),
            br(), card(card_header("Individual interval AUC by administration"),tableOutput("ap_auc_interval_table"))),
          tabPanel("Dosage adjustment", br(),
            card(card_header("Dosage modification simulation"),
              fluidRow(
                column(3, numericInput("expert_new_dose","Dose initiale (mg)",value=400,min=0,step=25)),
                column(3, numericInput("expert_new_tau","Interval (h)",value=24,min=6,max=48,step=6)),
                column(3, numericInput("expert_n_days","Simulation days",value=7,min=2,max=30,step=1)),
                column(3, br(), actionButton("fill_uniform","Fill all days",class="btn-outline-secondary btn-sm",style="width:100%;margin-top:4px;"))
              ),
              tags$b("Simulated dose per day"),
              uiOutput("expert_day_doses"),
              br(),
              fluidRow(
                column(6, actionButton("run_expert_sim","Simuler cette posologie",class="btn-warning",style="margin-bottom:10px;width:100%;")),
                column(6, actionButton("find_optimal_dose","Find optimal regimen",class="btn-success",style="margin-bottom:10px;width:100%;background:#2E7D32;border:none;"))
              ),
              uiOutput("expert_recommendation"),
              plotOutput("expert_pk_plot",height=450),
              br(), tableOutput("expert_auc_table")
            ))
        ), br(), uiOutput("report_ui_aposteriori"), br(), disclaimer_box()
      )
    ))
  )
)

# ==============================================================================
# SERVER
# ==============================================================================
server <- function(input, output, session) {
  local({
    candidates <- c("flucoshiny.png","Flucoshiny.png","FlucoShiny.png","FLUCOSHINY.png")
    src <- candidates[file.exists(candidates)]
    if (length(src) > 0 && !file.exists("www/flucoshiny.png")) {
      if (!dir.exists("www")) dir.create("www")
      file.copy(src[1], "www/flucoshiny.png", overwrite=TRUE)
    }
  })

  # Pharmacologist modal at startup
  showModal(modalDialog(title="Welcome to FlucoShiny",
    p("Please enter your name for report validation."),
    textInput("pharma_name","Pharmacologist name",value="",placeholder="Dr First Last"),
    footer=actionButton("confirm_pharma","Confirm",class="btn-primary"), easyClose=FALSE))
  observeEvent(input$confirm_pharma, {
    rv$pharmacologue <- trimws(input$pharma_name)
    removeModal()
  }, ignoreInit=TRUE, once=TRUE)

  patient <- reactiveValues(last_name="",first_name="",dob=Sys.Date()-365*40,service="",ipp="",sexe="Male",uf="",hopital="",adresse="",filled=FALSE)
  rv <- reactiveValues(results=NULL,ap_results=NULL,dose_rec=NULL,ap_inputs=NULL,pending=NULL,pharmacologue="",ap_error_msg=NULL,expert_sim_trigger=NULL,ap_unlocked=FALSE,expert_rec=NULL)
  get_target_auc <- reactive({ if(isTRUE(input$glabrata)) return(1600); if(!is.na(input$mic_exact)&&is.finite(input$mic_exact)&&input$mic_exact>0) return(100*input$mic_exact); 400 })
  get_mic_label <- reactive({ if(isTRUE(input$glabrata)) return(HTML("CMI = 16 mg/L (clinical breakpoint <span style=\"white-space:nowrap\"><i>Candida glabrata</i></span>)")); if(!is.na(input$mic_exact)&&is.finite(input$mic_exact)&&input$mic_exact>0) return(paste0("CMI = ",input$mic_exact," mg/L")); HTML("CMI = 4 mg/L (clinical breakpoint <span style=\"white-space:nowrap\"><i>Candida albicans</i></span>, <span style=\"white-space:nowrap\"><i>Candida parapsilosis</i></span>, <span style=\"white-space:nowrap\"><i>Candida tropicalis</i></span>)") })

  output$target_text <- renderUI({
    if(!is.na(input$mic_exact)&&is.finite(input$mic_exact)&&input$mic_exact>0) HTML(paste0("PK/PD target: AUC ",100*input$mic_exact," mg\u00b7h/L<br>AUC/MIC >100 with MIC ",input$mic_exact," mg/L"))
    else if(isTRUE(input$glabrata)) HTML("PK/PD target: AUC 1600 mg\u00b7h/L<br>AUC/MIC &gt;100 \u2014 MIC 16 mg/L (clinical breakpoint for <span style=\"white-space:nowrap\"><i>Candida glabrata</i></span>)")
    else HTML("PK/PD target: AUC 400 mg\u00b7h/L<br>AUC/MIC &gt;100 \u2014 MIC 4 mg/L (clinical breakpoint for <span style=\"white-space:nowrap\"><i>Candida albicans</i></span>, <span style=\"white-space:nowrap\"><i>Candida parapsilosis</i></span>, <span style=\"white-space:nowrap\"><i>Candida tropicalis</i></span>)")
  })

  get_clinical <- reactive(list(wt=input$wt, ckd=input$ckd, cirrh=ifelse(input$cirrh=="Yes",1,0), dob=if(!is.na(patient$dob))format(patient$dob,"%d/%m/%Y")else NA, mic_exact=input$mic_exact, glabrata=isTRUE(input$glabrata)))

  dose_choice_text <- function(label,wt){if(grepl("mg/kg/day",label,fixed=TRUE)){mgkg<-suppressWarnings(as.numeric(sub(" .*","",label)));if(is.na(mgkg))return(label);paste0(label," (",round(mgkg*wt)," mg/day for ",wt," kg)")}else label}
  dose_status_line <- function(case){if(identical(case,"optimal"))"<b>Optimal</b> (PTA \u2265 90%)" else if(identical(case,"suboptimal"))"<b>Suboptimal</b> (PTA \u2265 75%) \u2013 <b>Risk of underexposure</b>" else "<b>Failure</b> (PTA < 75%)"}

  patient_badge_html <- reactive({if(!patient$filled)return(NULL);div(class="pat-badge",HTML(paste0("<b>",patient$last_name," ",patient$first_name,"</b> &mdash; ",if(!is.na(patient$dob))format(patient$dob,"%d/%m/%Y")else "?"," &mdash; IPP : ",patient$ipp," &mdash; ",patient$service)))})
  output$patient_badge_ap <- renderUI(patient_badge_html())
  show_patient_modal <- function() showModal(patient_modal_ui(patient))
  observeEvent(input$goto_apriori, updateTabsetPanel(session, "main_tabs", "A priori"), ignoreInit=TRUE)
  observeEvent(input$goto_aposteriori, {
    showModal(modalDialog(title="A Posteriori Module Access",
      passwordInput("ap_password","Password",placeholder="Enter password"),
      footer=tagList(modalButton("Cancel"), actionButton("confirm_ap_pwd","Confirm",class="btn-primary")), easyClose=FALSE))
  }, ignoreInit=TRUE)
  observeEvent(input$confirm_ap_pwd, {
    if (identical(input$ap_password, AP_PASSWORD)) { removeModal(); rv$ap_unlocked <- TRUE; updateTabsetPanel(session, "main_tabs", "A posteriori") }
    else showNotification("Incorrect password", type="error", duration=3)
  }, ignoreInit=TRUE)

  # Block bypass: if tab changes to A posteriori without unlock, redirect
  observeEvent(input$main_tabs, {
    if (identical(input$main_tabs, "A posteriori") && !isTRUE(rv$ap_unlocked)) {
      updateTabsetPanel(session, "main_tabs", "Home")
      showNotification("Password-protected access", type="warning")
    }
  }, ignoreInit=TRUE)
  observeEvent(input$back_home_ap, updateTabsetPanel(session, "main_tabs", "Home"), ignoreInit=TRUE)
  observeEvent(input$back_home_post, updateTabsetPanel(session, "main_tabs", "Home"), ignoreInit=TRUE)

  # Reset buttons
  observeEvent(input$reset_ap, {
    rv$results <- NULL
    patient$filled <- FALSE; patient$last_name <- ""; patient$first_name <- ""
    patient$ipp <- ""; patient$service <- ""; patient$sexe <- "Male"
    patient$uf <- ""; patient$hopital <- ""; patient$adresse <- ""; patient$dob <- Sys.Date()-365*40
    updateNumericInput(session, "wt", value=70); updateNumericInput(session, "ckd", value=90)
    updateSelectInput(session, "cirrh", selected="No"); updateNumericInput(session, "tinf", value=0.5)
    updateCheckboxInput(session, "glabrata", value=FALSE); updateNumericInput(session, "mic_exact", value=NA)
    updateTabsetPanel(session, "main_tabs", "Home")
    showNotification("Everything reset", type="message")
  }, ignoreInit=TRUE)

  observeEvent(input$reset_post, {
    rv$ap_results <- NULL; rv$dose_rec <- NULL; rv$ap_inputs <- NULL; rv$ap_error_msg <- NULL
    rv$expert_sim_trigger <- NULL; rv$expert_rec <- NULL
    rv_mono$data <- NULL; rv_mono$patients <- NULL
    rv_fields$doses <- list(); rv_fields$concs <- list()
    patient$filled <- FALSE; patient$last_name <- ""; patient$first_name <- ""
    patient$ipp <- ""; patient$service <- ""; patient$sexe <- "Male"
    patient$uf <- ""; patient$hopital <- ""; patient$adresse <- ""; patient$dob <- Sys.Date()-365*40
    updateNumericInput(session, "wt", value=70); updateNumericInput(session, "ckd", value=90)
    updateSelectInput(session, "cirrh", selected="No"); updateNumericInput(session, "tinf", value=0.5)
    updateCheckboxInput(session, "glabrata", value=FALSE); updateNumericInput(session, "mic_exact", value=NA)
    updateCheckboxInput(session, "ap_regular", value=TRUE)
    updateDateInput(session, "ap_date_first", value=Sys.Date()-7); updateNumericInput(session, "ap_heure_first", value=8)
    updateNumericInput(session, "ap_dose_first", value=800)
    updateDateInput(session, "ap_date_maint", value=Sys.Date()-6); updateNumericInput(session, "ap_heure_maint", value=8)
    updateNumericInput(session, "ap_dose_maint", value=400)
    updateDateInput(session, "ap_date_last", value=Sys.Date()-1)
    updateNumericInput(session, "ap_n_conc", value=1); updateNumericInput(session, "ap_n_doses", value=3)
    updateNumericInput(session, "expert_new_dose", value=400); updateNumericInput(session, "expert_new_tau", value=24)
    updateNumericInput(session, "expert_n_days", value=7)
    showNotification("Everything reset", type="message")
  }, ignoreInit=TRUE)

  # ===================== MONOLIX FILE UPLOAD =====================
  rv_mono <- reactiveValues(data=NULL, patients=NULL)

  observeEvent(input$monolix_file, {
    if(is.null(input$monolix_file))return()
    tryCatch({
      if (!requireNamespace("readxl", quietly=TRUE)) { install.packages("readxl", repos="https://cloud.r-project.org", quiet=TRUE) }
      df <- readxl::read_excel(input$monolix_file$datapath)
      names(df) <- trimws(names(df))
      # Standardize column names
      nm <- tolower(names(df))
      id_col <- which(nm=="id")[1]
      dt_col <- which(nm %in% c("datetime","date","dt"))[1]
      if(is.na(dt_col)) dt_col <- 2  # fallback: column 2 is always DateTime
      amt_col <- which(nm=="amt")[1]
      obs_col <- which(nm=="obs")[1]; bw_col <- which(nm %in% c("bw","wt","weight"))[1]
      ckd_col <- which(nm %in% c("ckd","ckd-epi","ckdepi"))[1]; cirrh_col <- which(nm %in% c("cirrhosis","cirrh","cirrhosis"))[1]
      tinf_col <- which(nm %in% c("infusion duration","tinf","inf_dur"))[1]
      if (is.na(id_col)||is.na(amt_col)||is.na(obs_col)) { showNotification("Columns ID, AMT, OBS required", type="error"); return() }
      df$ID_val <- suppressWarnings(as.numeric(as.character(df[[id_col]])))
      # DateTime: readxl returns POSIXct in UTC - format in UTC to preserve literal values, reparse as local
      dt_raw <- df[[dt_col]]
      if(inherits(dt_raw, "POSIXct")) {
        df$DT_val <- as.POSIXct(format(dt_raw, "%Y-%m-%d %H:%M:%S", tz="UTC"), format="%Y-%m-%d %H:%M:%S")
      } else if(inherits(dt_raw, "character") || inherits(dt_raw, "factor")) {
        dt_str <- as.character(dt_raw)
        df$DT_val <- as.POSIXct(dt_str, format="%m/%d/%Y %H:%M:%S")
        if(all(is.na(df$DT_val))) df$DT_val <- as.POSIXct(dt_str, format="%d/%m/%Y %H:%M:%S")
        if(all(is.na(df$DT_val))) df$DT_val <- as.POSIXct(dt_str, format="%Y-%m-%d %H:%M:%S")
        if(all(is.na(df$DT_val))) df$DT_val <- as.POSIXct(dt_str)
      } else {
        # numeric (Excel serial) or other - try conversion
        df$DT_val <- as.POSIXct(as.numeric(dt_raw)*86400, origin="1899-12-30")
      }
      df$AMT_val <- suppressWarnings(as.numeric(df[[amt_col]])); df$AMT_val[is.na(df$AMT_val)] <- 0
      df$OBS_val <- suppressWarnings(as.numeric(df[[obs_col]])); df$OBS_val[is.na(df$OBS_val)] <- NA
      df$BW_val <- if(!is.na(bw_col)) suppressWarnings(as.numeric(df[[bw_col]])) else NA
      df$CKD_val <- if(!is.na(ckd_col)) suppressWarnings(as.numeric(df[[ckd_col]])) else NA
      df$CIRRH_val <- if(!is.na(cirrh_col)) suppressWarnings(as.numeric(df[[cirrh_col]])) else 0
      df$TINF_val <- if(!is.na(tinf_col)) suppressWarnings(as.numeric(df[[tinf_col]])) else 0.5
      df$TINF_val[is.na(df$TINF_val)] <- 0.5
      rv_mono$data <- df; rv_mono$patients <- sort(unique(df$ID_val[!is.na(df$ID_val)]))
      showNotification(paste0("File loaded: ",length(rv_mono$patients)," patients"), type="message")
    }, error=function(e) showNotification(paste0("Read error: ",e$message), type="error"))
  })

  output$monolix_patient_select <- renderUI({
    if (is.null(rv_mono$patients)) return(NULL)
    tagList(
      selectInput("mono_patient_id","Patient ID",choices=rv_mono$patients),
      actionButton("load_mono_patient","Load this patient",class="btn-warning btn-sm",icon=icon("download"))
    )
  })

  output$monolix_status <- renderUI({
    if (is.null(rv_mono$data) || is.null(input$mono_patient_id)) return(NULL)
    pid <- as.numeric(input$mono_patient_id)
    dd <- rv_mono$data[rv_mono$data$ID_val==pid & !is.na(rv_mono$data$ID_val),]
    n_doses <- sum(dd$AMT_val>0, na.rm=TRUE); n_obs <- sum(!is.na(dd$OBS_val), na.rm=TRUE)
    div(style="font-size:13px;color:#444;padding:4px 0;", HTML(paste0("Patient ", pid, " : <b>", n_doses, " administrations</b>, <b>", n_obs, " observations</b>")))
  })

  observeEvent(input$load_mono_patient, {
    if (is.null(rv_mono$data) || is.null(input$mono_patient_id)) return(NULL)
    pid <- as.numeric(input$mono_patient_id)
    dd <- rv_mono$data[rv_mono$data$ID_val==pid & !is.na(rv_mono$data$ID_val),]
    dd <- dd[order(dd$DT_val),]
    doses <- dd[!is.na(dd$AMT_val) & dd$AMT_val > 0, ]
    obs <- dd[!is.na(dd$OBS_val) & dd$OBS_val > 0, ]
    if (nrow(doses)==0 || nrow(obs)==0) { showNotification("No doses or observations for this patient", type="error"); return() }
    first_dose <- doses[1, ]; last_dose <- doses[nrow(doses), ]
    first_amt <- first_dose$AMT_val
    maint_amt <- if(nrow(doses)>1) doses$AMT_val[nrow(doses)] else first_amt
    bw <- median(dd$BW_val, na.rm=TRUE); ckd <- median(dd$CKD_val, na.rm=TRUE)
    cirrh <- max(dd$CIRRH_val, na.rm=TRUE); tinf <- median(dd$TINF_val, na.rm=TRUE)
    dose_dts <- as.numeric(doses$DT_val)
    dose_diffs_h <- diff(dose_dts) / 3600
    maint_doses <- if(nrow(doses)>1) doses$AMT_val[2:nrow(doses)] else c()
    doses_constant <- length(maint_doses)==0 || length(unique(maint_doses))==1
    is_regular <- length(dose_diffs_h) > 0 && all(abs(dose_diffs_h - 24) < 4) && doses_constant

    # Fill basic fields immediately
    updateNumericInput(session, "wt", value=round(bw,1))
    updateNumericInput(session, "ckd", value=round(ckd,1))  # global CKD (median)
    updateSelectInput(session, "cirrh", selected=if(cirrh>0)"Yes" else "No")
    updateNumericInput(session, "tinf", value=round(tinf,2))
    # Store per-obs CKD for later pre-fill (after ap_n_conc is set)
    rv_mono$ckd_per_obs <- if(nrow(obs)>0 && !is.na(ckd_col)) {
      sapply(seq_len(nrow(obs)), function(ii) {
        t_obs <- obs$DT_val[ii]
        # Find closest CKD measurement in full data at or before obs time
        dd_sorted <- dd[order(dd$DT_val),]
        before <- dd_sorted[dd_sorted$DT_val <= t_obs & is.finite(dd_sorted$CKD_val),]
        if(nrow(before)>0) tail(before$CKD_val,1) else ckd
      })
    } else rep(ckd, nrow(obs))
    # Loading dose
    updateDateInput(session, "ap_date_first", value=as.Date(format(first_dose$DT_val, "%Y-%m-%d")))
    h_first <- as.numeric(format(first_dose$DT_val, "%H")) + as.numeric(format(first_dose$DT_val, "%M"))/60
    updateNumericInput(session, "ap_heure_first", value=round(h_first,1))
    updateNumericInput(session, "ap_dose_first", value=first_amt)
    # Entretien
    if (nrow(doses) > 1) {
      second_dose <- doses[2, ]
      updateDateInput(session, "ap_date_maint", value=as.Date(format(second_dose$DT_val, "%Y-%m-%d")))
      h_maint <- as.numeric(format(second_dose$DT_val, "%H")) + as.numeric(format(second_dose$DT_val, "%M"))/60
      updateNumericInput(session, "ap_heure_maint", value=round(h_maint,1))
    }
    updateNumericInput(session, "ap_dose_maint", value=round(maint_amt))
    # Last administration
    updateDateInput(session, "ap_date_last", value=as.Date(format(last_dose$DT_val, "%Y-%m-%d")))

    if (is_regular) {
      updateCheckboxInput(session, "ap_regular", value=TRUE)
    } else {
      updateCheckboxInput(session, "ap_regular", value=FALSE)
      # Store individual doses in rv_fields for renderUI
      dose_vals <- list()
      for(j in seq_len(nrow(doses))) {
        dose_vals[[as.character(j)]] <- list(
          d=as.Date(format(doses$DT_val[j], "%Y-%m-%d")),
          h=round(as.numeric(format(doses$DT_val[j],"%H")) + as.numeric(format(doses$DT_val[j],"%M"))/60, 1),
          a=doses$AMT_val[j])
      }
      rv_fields$doses <- dose_vals
      updateNumericInput(session, "ap_n_doses", value=nrow(doses))
    }
    # Store concentrations in rv_fields for renderUI
    conc_vals <- list()
    for(j in seq_len(nrow(obs))) {
      ckd_j <- if(!is.null(rv_mono$ckd_per_obs) && j<=length(rv_mono$ckd_per_obs)) rv_mono$ckd_per_obs[j] else ckd
      conc_vals[[as.character(j)]] <- list(
        d=as.Date(format(obs$DT_val[j], "%Y-%m-%d")),
        h=round(as.numeric(format(obs$DT_val[j],"%H")) + as.numeric(format(obs$DT_val[j],"%M"))/60, 1),
        c=obs$OBS_val[j],
        ckd=round(ckd_j, 1))
    }
    rv_fields$concs <- conc_vals
    updateNumericInput(session, "ap_n_conc", value=nrow(obs))

    showNotification(paste0("Patient ",pid," loaded: ",nrow(doses)," doses (",first_amt,"\u2192",maint_amt," mg), ",nrow(obs)," obs. ",
      "1st admin: ",format(first_dose$DT_val,"%d/%m/%Y %H:%M"),", last: ",format(last_dose$DT_val,"%d/%m/%Y %H:%M")),
      type="message", duration=10)
  })

  observeEvent(input$run,{rv$pending<-"apriori";show_patient_modal()},ignoreInit=TRUE)
  observeEvent(input$run_ap,{rv$pending<-"aposteriori";show_patient_modal()},ignoreInit=TRUE)
  observeEvent(input$confirm_patient,{
    if(!isTRUE(nzchar(trimws(input$pat_last_name)))||!isTRUE(nzchar(trimws(input$pat_first_name)))||is.null(input$pat_dob)||!isTRUE(nzchar(trimws(input$pat_service)))||!isTRUE(nzchar(trimws(input$pat_ipp)))){showNotification("Please fill in all required fields",type="warning");return()}
    patient$last_name<-trimws(input$pat_last_name);patient$first_name<-trimws(input$pat_first_name);patient$dob<-as.Date(input$pat_dob);patient$service<-trimws(input$pat_service);patient$ipp<-trimws(input$pat_ipp)
    patient$sexe<-input$pat_sexe;patient$uf<-trimws(input$pat_uf);patient$hopital<-trimws(input$pat_hopital);patient$adresse<-trimws(input$pat_adresse)
    patient$filled<-TRUE;removeModal()
    if(identical(rv$pending,"apriori"))run_apriori();if(identical(rv$pending,"aposteriori"))run_aposteriori();rv$pending<-NULL
  },ignoreInit=TRUE)

  # ===================== A PRIORI =====================
  run_apriori <- function(){
    age_years<-calc_age_years(patient$dob)
    if(!is.na(age_years)&&age_years<AGE_ADULT){
      showModal(modalDialog(title="Model not applicable",
        HTML(paste0("<span class='extrapWarn'><b>Model restricted to adults (\u2265 ",AGE_ADULT,
          " years). Patient aged ",round(age_years,1)," years.</b><br><br>",
          "A priori simulation is not available for paediatric patients.</span>")),
        footer=modalButton("OK"), easyClose=TRUE))
      return()  # BLOCK simulation
    }
    cirrh<-ifelse(input$cirrh=="Yes",1,0);target_auc<-get_target_auc()
    withProgress(message="Monte Carlo simulation...",value=0.2,{
      set.seed(123);pars<-simulate_params(input$wt,input$ckd,cirrh);setProgress(0.4)
      load_df<-loading_regimens(input$wt);maint_df<-maintenance_regimens(input$wt)

      # ── Loading dose: AUC per day, vectorised correctly over 200k patients ──
      # auc_loading_day_scalar: scalar CL/V version for mapply
      auc_load_scalar <- function(dose, CL, V, tinf_v, day_idx) {
        tau <- 24
        dt <- seq(0, by=tau, length.out=day_idx)
        da <- rep(dose, day_idx)
        ts <- (day_idx-1)*tau; te <- day_idx*tau
        k <- CL/V; auc <- 0
        for(i in seq_along(dt)){
          td <- dt[i]; D <- da[i]; R <- D/tinf_v
          ds0 <- ts-td; de0 <- te-td
          if(de0<=0) next
          ds0 <- max(ds0,0)
          if(ds0<tinf_v){ t1<-ds0; t2<-min(de0,tinf_v); auc<-auc+(R/CL)*((t2-t1)-(1/k)*(exp(-k*t1)-exp(-k*t2))) }
          if(de0>tinf_v){ t1<-max(ds0,tinf_v); t2<-de0; auc<-auc+(R/CL)*(1-exp(-k*tinf_v))*(1/k)*(exp(-k*(t1-tinf_v))-exp(-k*(t2-tinf_v))) }
        }
        auc
      }

      compute_load_auc <- function(dose, CL_vec, V_vec, tinf_v, day_idx) {
        mapply(function(cl,v) auc_load_scalar(dose, cl, v, tinf_v, day_idx), CL_vec, V_vec)
      }

      # Use full NSIM for J1 (fast: vectorised scalar), subsample for J2/J3 (3 nested doses)
      N_MULTI <- min(20000L, length(pars$CL))  # 20k still very precise for PTA
      set.seed(42)
      idx_multi <- sample(length(pars$CL), N_MULTI)
      CL_multi <- pars$CL[idx_multi]; V_multi <- pars$V[idx_multi]

      # J1: single dose 0-24h — use fast vectorised analytical formula
      load_auc_d1 <- sapply(load_df$Dose, function(d)
        auc_0_24_inf(d, pars$CL, pars$V, input$tinf, tend=24))
      if(!is.matrix(load_auc_d1)) load_auc_d1 <- matrix(load_auc_d1, nrow=length(pars$CL))
      load_sum_d1 <- summ_auc(load_auc_d1)
      load_pta_d1 <- colMeans(load_auc_d1 >= target_auc)

      # ── If any dose fails PTA≥90% on J1, compute J2 and J3 ──
      need_multiday <- any(load_pta_d1 < PTA_OPTIMAL)
      load_auc_d2 <- load_auc_d3 <- NULL
      load_sum_d2 <- load_sum_d3 <- NULL
      load_pta_d2 <- load_pta_d3 <- NULL
      if(need_multiday){
        setProgress(0.6)
        load_auc_d2 <- sapply(load_df$Dose, function(d)
          compute_load_auc(d, CL_multi, V_multi, input$tinf, 2L))
        if(!is.matrix(load_auc_d2)) load_auc_d2 <- matrix(load_auc_d2, nrow=N_MULTI)
        load_sum_d2 <- summ_auc(load_auc_d2)
        load_pta_d2 <- colMeans(load_auc_d2 >= target_auc)

        load_auc_d3 <- sapply(load_df$Dose, function(d)
          compute_load_auc(d, CL_multi, V_multi, input$tinf, 3L))
        if(!is.matrix(load_auc_d3)) load_auc_d3 <- matrix(load_auc_d3, nrow=N_MULTI)
        load_sum_d3 <- summ_auc(load_auc_d3)
        load_pta_d3 <- colMeans(load_auc_d3 >= target_auc)
      }
      setProgress(0.8)

      # Maintenance
      maint_auc <- sapply(maint_df$Dose, function(d) d/pars$CL)
      if(!is.matrix(maint_auc)) maint_auc <- matrix(maint_auc, ncol=nrow(maint_df))
      maint_sum <- summ_auc(maint_auc); maint_pta <- colMeans(maint_auc >= target_auc)

      decide <- function(pta,dv){if(any(pta>=PTA_OPTIMAL)){ok<-which(pta>=PTA_OPTIMAL);list(case="optimal",idx=ok[which.min(dv[ok])])}else if(any(pta>=PTA_ACCEPTABLE)){ok<-which(pta>=PTA_ACCEPTABLE);list(case="suboptimal",idx=ok[which.min(dv[ok])])}else list(case="failure",idx=NA_integer_)}

      # Best loading recommendation: first day that achieves PTA≥90%, or best available
      ld_d1 <- decide(load_pta_d1, load_df$Dose)
      ld <- ld_d1
      best_load_day <- 1L
      if(ld_d1$case != "optimal" && need_multiday){
        ld_d2 <- decide(load_pta_d2, load_df$Dose)
        if(ld_d2$case == "optimal"){ ld <- ld_d2; best_load_day <- 2L }
        else {
          ld_d3 <- decide(load_pta_d3, load_df$Dose)
          if(ld_d3$case == "optimal"){ ld <- ld_d3; best_load_day <- 3L }
          else if(ld_d2$case == "suboptimal"){ ld <- ld_d2; best_load_day <- 2L }
          else if(ld_d3$case == "suboptimal"){ ld <- ld_d3; best_load_day <- 3L }
        }
      }
      md <- decide(maint_pta, maint_df$Dose)
      setProgress(1)
      rv$results <- list(
        target_auc=target_auc, mic_label=get_mic_label(),
        CLtyp=pars$CLtyp, Vtyp=pars$Vtyp,
        load_df=load_df, maint_df=maint_df,
        # J1
        load_sum_d1=load_sum_d1, load_pta_d1=load_pta_d1,
        # J2/J3 (NULL if not needed)
        load_sum_d2=load_sum_d2, load_pta_d2=load_pta_d2,
        load_sum_d3=load_sum_d3, load_pta_d3=load_pta_d3,
        need_multiday=need_multiday, best_load_day=best_load_day,
        # Keep legacy names pointing to J1 for backward compat (report etc.)
        load_sum=load_sum_d1, load_pta=load_pta_d1,
        maint_sum=maint_sum, maint_pta=maint_pta,
        load_case=ld$case, maint_case=md$case,
        rec_load_idx=ld$idx, rec_maint_idx=md$idx,
        wt=input$wt, ckd=input$ckd, tinf=input$tinf,
        cirrh=cirrh, mic_exact=input$mic_exact, glabrata=isTRUE(input$glabrata),
        wt_extrap=(input$wt<WT_LOWER)||(input$wt>WT_UPPER),
        ckd_extrap=(input$ckd<CKD_LOWER)||(input$ckd>CKD_UPPER),
        age_years=age_years,
        global_failure=identical(ld$case,"failure")&&identical(md$case,"failure"))
    })
  }

  output$card_pk<-renderUI({if(is.null(rv$results))return(NULL);r<-rv$results;ww<-if(isTRUE(r$wt_extrap))HTML(paste0("<span class='extrapWarn'><b>Extrapolated: weight outside ",WT_LOWER,"\u2013",WT_UPPER," kg</b></span>"))else HTML("");cw<-if(isTRUE(r$ckd_extrap))HTML(paste0("<span class='extrapWarn'><b>Extrapolated: CKD-EPI outside ",CKD_LOWER,"\u2013",CKD_UPPER," mL/min/1.73m\u00b2</b></span>"))else HTML("");card(style="background-color:#E3F2FD;",card_header("PK summary"),tagList(p(HTML(paste0("<b>Target AUC:</b> ",r$target_auc," mg\u00b7h/L | <b>Infusion:</b> ",r$tinf," h | <b>CL:</b> ",round(r$CLtyp,3)," L/h | <b>V:</b> ",round(r$Vtyp,1)," L"))),if(isTRUE(r$wt_extrap))p(ww),if(isTRUE(r$ckd_extrap))p(cw)))})

  # ── Stale-inputs warning ─────────────────────────────────────────────────────
  output$ap_stale_warning <- renderUI({
    if(is.null(rv$results)) return(NULL)
    r <- rv$results

    # Detect any change vs the snapshot stored at last run
    cur_cirrh  <- ifelse(input$cirrh == "Yes", 1, 0)
    cur_mic    <- input$mic_exact
    cur_glab   <- isTRUE(input$glabrata)
    snap_mic   <- if(!is.null(r$mic_exact)) r$mic_exact else NA_real_
    snap_glab  <- if(!is.null(r$glabrata)) isTRUE(r$glabrata) else FALSE
    snap_cirrh <- if(!is.null(r$cirrh))    r$cirrh else 0

    wt_changed    <- !isTRUE(all.equal(input$wt,    r$wt,    tolerance = 1e-6))
    ckd_changed   <- !isTRUE(all.equal(input$ckd,   r$ckd,   tolerance = 1e-6))
    tinf_changed  <- !isTRUE(all.equal(input$tinf,  r$tinf,  tolerance = 1e-6))
    cirrh_changed <- !isTRUE(all.equal(cur_cirrh,   snap_cirrh, tolerance = 1e-6))
    mic_changed   <- !isTRUE(all.equal(
      if(is.na(cur_mic)||!is.finite(cur_mic)) NA_real_ else cur_mic,
      if(is.na(snap_mic)||!is.finite(snap_mic)) NA_real_ else snap_mic,
      tolerance = 1e-6))
    glab_changed  <- !identical(cur_glab, snap_glab)

    if(!any(wt_changed, ckd_changed, tinf_changed, cirrh_changed, mic_changed, glab_changed))
      return(NULL)

    # Build human-readable summary of the STORED (valid) values
    mic_disp <- if(isTRUE(snap_glab)) "MIC 16 mg/L (C. glabrata)" else
                if(!is.na(snap_mic) && is.finite(snap_mic) && snap_mic > 0)
                  paste0("MIC ", snap_mic, " mg/L") else
                "MIC 4 mg/L (C. albicans / C. tropicalis)"

    div(
      style = paste0(
        "margin:0 0 10px 0;padding:12px 16px;",
        "background:linear-gradient(135deg,#FFF8E1,#FFF3CD);",
        "border:1.5px solid #F9A825;border-left:5px solid #E65100;",
        "border-radius:10px;font-size:13px;line-height:1.65;",
        "box-shadow:0 1px 4px rgba(230,81,0,0.10);"
      ),
      HTML(paste0(
        "<b style='color:#E65100;font-size:14px;'>",
        "&#9888;&nbsp;Inputs changed \u2014 results below are NOT updated",
        "</b><br>",
        "<span style='color:#555555;'>",
        "Displayed results are valid for: ",
        "<b>weight ", r$wt, " kg</b>, ",
        "<b>CKD-EPI ", r$ckd, " mL/min</b>, ",
        "<b>infusion ", r$tinf, " h</b>, ",
        "<b>", mic_disp, "</b>, ",
        "cirrhosis: <b>", if(snap_cirrh == 1) "Yes" else "No", "</b>.",
        "<br>Click <em>Run a priori PK/PD evaluation</em> again to update with the new values.",
        "</span>"
      ))
    )
  })

  output$alert_global<-renderUI({if(is.null(rv$results))return(NULL);r<-rv$results;if(isTRUE(r$global_failure))card(style="background:#B71C1C;color:white;",card_header("PK/PD Failure"),p("Loading dose: no dose achieves PTA \u2265 75%."),p("Maintenance dose: no dose achieves PTA \u2265 75%."),strong("Consider an antifungal alternative."))})
  output$box_loading<-renderUI({
    if(is.null(rv$results))return(NULL);r<-rv$results
    if(isTRUE(r$global_failure))return(NULL)

    # Labels for each day
    day_labels   <- c("Day 1 (0\u201324h)","Day 2 (24\u201348h)","Day 3 (48\u201372h)")
    day_ranges   <- c("0\u201324h","0\u201348h","0\u201372h")  # cumulative range for N days
    best_lbl     <- day_labels[min(r$best_load_day,3L)]
    n_days       <- r$best_load_day

    # Failure across all 3 days
    if(r$load_case=="failure") return(card(
      style="background-color:#FDECEA;border-left:5px solid #B71C1C;",
      card_header(HTML("&#10007; Loading dose \u2013 No adequate dose identified")),
      p(HTML(paste0(
        "<b style='color:#B71C1C;'>No loading dose achieves PTA \u2265 75% for the PK/PD target ",
        "(AUC \u2265 ",r$target_auc," mg\u00b7h/L) over 3 consecutive days (0\u201372h).</b>"
      ))),
      p("Expert pharmacological opinion is required. Consider an alternative antifungal.")))

    txt <- dose_choice_text(r$load_df$DoseLabel[r$rec_load_idx],r$wt)

    # Case: repeated loading dose needed (best PTA on Day 2 or Day 3)
    if(n_days > 1) {
      cum_range <- day_ranges[min(n_days,3L)]
      subopt_note <- if(r$load_case=="suboptimal") HTML(paste0(
        "<div style='margin-top:8px;padding:8px 12px;border-radius:6px;background:#FFF8E1;border:1px solid #F9A825;font-size:12.5px;color:#555;'>",
        "<b style='color:#E65100;'>&#9888; Suboptimal PTA:</b> Even after ",n_days,
        " days of repeated loading, PTA remains between 75\u201390%. ",
        "Consider a higher dose or an alternative antifungal agent.</div>")) else NULL
      return(card(
        style="background-color:#FFF3CD;border-left:5px solid #E65100;",
        card_header(HTML("&#9888; Loading dose \u2013 Repeated loading dose required")),
        tagList(
          div(style="color:#B71C1C;font-weight:700;font-size:13.5px;margin-bottom:10px;",
            HTML(paste0(
              "&#9888;&nbsp;No loading dose achieves PTA \u2265 90% for the PK/PD target ",
              "(AUC \u2265 ",r$target_auc," mg\u00b7h/L) with a <u>single administration</u> (Day 1 alone)."
            ))),
          div(style="font-size:14px;margin-bottom:8px;",
            HTML(paste0(
              "&#8594;&nbsp;<b>Recommended loading dose: ",
              "<span style='font-size:16px;color:#E65100;'>",txt,"</span>",
              " repeated for <span style='color:#E65100;font-size:16px;font-weight:900;'>",
              n_days," consecutive days</span>",
              " (",cum_range,").</b>"
            ))),
          div(style="font-size:13px;color:#444;",
            HTML(paste0(
              "&#10003;&nbsp;PTA \u2265 90% is first achieved on <b>",best_lbl,"</b>, ",
              "after cumulated dosing over ",n_days," day",if(n_days>1)"s" else "","."
            ))),
          subopt_note
        )
      ))
    }

    # Case: single loading dose sufficient (Day 1)
    if(r$load_case=="optimal")
      card(style="background-color:#E8F5E9;border-left:5px solid #2E7D32;",
        card_header(HTML("&#10003; Loading dose \u2013 Optimal")),
        tagList(p(HTML(dose_status_line("optimal"))),
                p(HTML(paste0("<b>Loading dose</b>: ",txt)))))
    else
      card(style="background-color:#FFF3E0;border-left:5px solid #F57C00;",
        card_header("Loading dose \u2013 Suboptimal"),
        tagList(p(HTML(dose_status_line("suboptimal"))),
                p(HTML(paste0("<b>Loading dose</b>: ",txt)))))
  })
  output$box_maintenance<-renderUI({if(is.null(rv$results))return(NULL);r<-rv$results;if(isTRUE(r$global_failure))return(NULL);if(r$maint_case=="failure")return(card(style="background-color:#FDECEA;",card_header("Maintenance dose \u2013 failure"),p(HTML(dose_status_line("failure"))),p("No dose achieves PTA \u2265 75%."),tags$em("Consider alternative.")));txt<-dose_choice_text(r$maint_df$DoseLabel[r$rec_maint_idx],r$wt);if(r$maint_case=="optimal")card(style="background-color:#E8F5E9;",card_header("Maintenance dose"),tagList(p(HTML(dose_status_line("optimal"))),p(HTML(paste0("<b>Maintenance dose</b>: ",txt)))))else card(style="background-color:#FFF3E0;",card_header("Maintenance dose"),tagList(p(HTML(dose_status_line("suboptimal"))),p(HTML(paste0("<b>Maintenance dose</b>: ",txt)))))})
  output$card_results<-renderUI({
    if(is.null(rv$results))return(NULL); r<-rv$results
    load_card <- if(isTRUE(r$need_multiday)){
      tabsetPanel(
        tabPanel("Day 1 (0\u201324h)",  plotOutput("plot_load_d1",height=280), tableOutput("table_load_d1")),
        tabPanel("Day 2 (24\u201348h)", plotOutput("plot_load_d2",height=280), tableOutput("table_load_d2")),
        tabPanel("Day 3 (48\u201372h)", plotOutput("plot_load_d3",height=280), tableOutput("table_load_d3"))
      )
    } else {
      tagList(plotOutput("plot_load_d1",height=280), tableOutput("table_load_d1"))
    }
    tagList(
      card(card_header(if(isTRUE(r$need_multiday))
        paste0("Loading dose \u2013 Multi-day AUC (PTA <90% on Day 1 \u2192 Days 1\u20133 simulated)")
        else "Loading dose \u2013 AUC 0\u201324 h (mg\u00b7h/L)"),
        load_card),
      br(),
      card(card_header("Maintenance dose \u2013 AUC 312\u2013336 h (Day 14)"),
           plotOutput("plot_maint",height=280), tableOutput("table_maint"))
    )
  })

  # ── Loading plot/table helpers ──
  .make_load_plot <- function(r, sum_d, pta_d, day_label) {
    if(is.null(sum_d)) return(NULL)
    plot_fun(sum_d, pta_d, r$target_auc, r$load_df$DoseLabel, r$rec_load_idx,
             paste0("AUC ", day_label, " (mg\u00b7h/L)"),
             paste0("Loading dose exposure \u2013 ", day_label))
  }
  .make_load_table <- function(r, sum_d, pta_d, day_label) {
    if(is.null(sum_d)) return(NULL)
    data.frame("Loading dose"=r$load_df$DoseLabel,
               "AUC median"=round(sum_d$median,1),
               "95% CI"=paste0(round(sum_d$p2.5,1),"\u2013",round(sum_d$p97.5,1)),
               "PTA \u226590%"=round(pta_d,3),
               "Day"=day_label,
               check.names=FALSE)
  }

  output$plot_load_d1 <- renderPlot({if(is.null(rv$results))return(NULL);r<-rv$results;.make_load_plot(r,r$load_sum_d1,r$load_pta_d1,"0\u201324h")})
  output$plot_load_d2 <- renderPlot({if(is.null(rv$results)||!isTRUE(rv$results$need_multiday))return(NULL);r<-rv$results;.make_load_plot(r,r$load_sum_d2,r$load_pta_d2,"24\u201348h")})
  output$plot_load_d3 <- renderPlot({if(is.null(rv$results)||!isTRUE(rv$results$need_multiday))return(NULL);r<-rv$results;.make_load_plot(r,r$load_sum_d3,r$load_pta_d3,"48\u201372h")})
  output$table_load_d1 <- renderTable({if(is.null(rv$results))return(NULL);r<-rv$results;.make_load_table(r,r$load_sum_d1,r$load_pta_d1,"0\u201324h")},align="c")
  output$table_load_d2 <- renderTable({if(is.null(rv$results)||!isTRUE(rv$results$need_multiday))return(NULL);r<-rv$results;.make_load_table(r,r$load_sum_d2,r$load_pta_d2,"24\u201348h")},align="c")
  output$table_load_d3 <- renderTable({if(is.null(rv$results)||!isTRUE(rv$results$need_multiday))return(NULL);r<-rv$results;.make_load_table(r,r$load_sum_d3,r$load_pta_d3,"48\u201372h")},align="c")
  output$table_maint<-renderTable({if(is.null(rv$results))return(NULL);r<-rv$results;data.frame("Maintenance dose"=r$maint_df$DoseLabel,"Median AUC (mg.h/L)"=round(r$maint_sum$median,1),"95% CI"=paste0(round(r$maint_sum$p2.5,1),"\u2013",round(r$maint_sum$p97.5,1)),"PTA"=round(r$maint_pta,3),check.names=FALSE)},align="c")
  output$plot_maint<-renderPlot({if(is.null(rv$results))return(NULL);r<-rv$results;plot_fun(r$maint_sum,r$maint_pta,r$target_auc,r$maint_df$DoseLabel,r$rec_maint_idx,"AUC 312-336 h (mg*h/L)","Maintenance dose exposure (Day 14)")})
  output$report_ui_apriori<-renderUI({if(is.null(rv$results))return(NULL);downloadButton("dl_apriori","Download a priori report (HTML)",class="btn-primary",style="width:100%;font-size:16px;padding:12px;")})

  output$pta_plot_load<-renderPlot({
    if(is.null(rv$results))return(NULL);r<-rv$results
    ord<-order(r$load_df$Dose,decreasing=TRUE)
    doses<-r$load_df$Dose[ord]
    pta1<-r$load_pta_d1[ord]*100
    par(mar=c(5,5,3,5),bg="white")
    plot(doses,pta1,type="n",xlab="Dose (mg)",ylab="PTA (%)",xlim=rev(range(doses)),ylim=c(0,105),las=1,
         main=if(isTRUE(r$need_multiday))"Loading dose - PTA by day (AUC 0-24h / 24-48h / 48-72h)" else "Loading dose - PTA (AUC 0-24h)")
    grid(col="#ECEFF1",lty=1)
    abline(h=90,col="#2E7D32",lty=2,lwd=2);abline(h=75,col="#F57C00",lty=2,lwd=1.5)
    polygon(c(doses,rev(doses)),c(rep(0,length(doses)),rev(pta1)),col=grDevices::adjustcolor("#1565C0",.10),border=NA)
    lines(doses,pta1,lwd=2.5,col="#1565C0"); points(doses,pta1,pch=19,col=ifelse(pta1>=90,"#2E7D32",ifelse(pta1>=75,"#F57C00","#C62828")),cex=1.2)
    if(isTRUE(r$need_multiday) && !is.null(r$load_pta_d2)){
      pta2<-r$load_pta_d2[ord]*100; pta3<-r$load_pta_d3[ord]*100
      lines(doses,pta2,lwd=2.5,col="#7B1FA2",lty=1); points(doses,pta2,pch=17,col=ifelse(pta2>=90,"#2E7D32",ifelse(pta2>=75,"#F57C00","#7B1FA2")),cex=1.2)
      lines(doses,pta3,lwd=2.5,col="#E65100",lty=1); points(doses,pta3,pch=15,col=ifelse(pta3>=90,"#2E7D32",ifelse(pta3>=75,"#F57C00","#E65100")),cex=1.2)
      legend("right",inset=c(-0.12,0),bty="o",box.col="#E0E0E0",bg="white",cex=0.78,
        legend=c("Day 1 (0-24h)","Day 2 (24-48h)","Day 3 (48-72h)","PTA \u226590%","PTA \u226575%"),
        lty=c(1,1,1,2,2),lwd=c(2.5,2.5,2.5,2,1.5),col=c("#1565C0","#7B1FA2","#E65100","#2E7D32","#F57C00"))
    } else {
      legend("right",bty="n",cex=0.8,legend=c("PTA \u2265 90% (optimal)","PTA \u2265 75% (acceptable)"),lty=c(2,2),lwd=c(2,1.5),col=c("#2E7D32","#F57C00"))
    }
  })

  output$pta_plot_maint<-renderPlot({
    if(is.null(rv$results))return(NULL);r<-rv$results
    ord<-order(r$maint_df$Dose,decreasing=TRUE)
    doses<-r$maint_df$Dose[ord]; pta<-r$maint_pta[ord]*100
    par(mar=c(5,5,3,2),bg="white")
    plot(doses,pta,type="n",xlab="Dose (mg)",ylab="PTA (%)",xlim=rev(range(doses)),ylim=c(0,105),las=1,main="Maintenance dose - PTA (AUC steady-state)",xaxt="n")
    axis(1, at=c(50,100,200,400,600,800,1000,1200), cex.axis=0.8)
    grid(col="#ECEFF1",lty=1)
    abline(h=90,col="#2E7D32",lty=2,lwd=2);abline(h=75,col="#F57C00",lty=2,lwd=1.5)
    polygon(c(doses,rev(doses)),c(rep(0,length(doses)),rev(pta)),col=grDevices::adjustcolor("#1565C0",.12),border=NA)
    lines(doses,pta,lwd=3,col="#1565C0")
    points(doses,pta,pch=19,col=ifelse(pta>=90,"#2E7D32",ifelse(pta>=75,"#F57C00","#C62828")),cex=1.3)
    text(doses,pta+4,paste0(round(pta,0),"%"),cex=0.65,col="#333333")
    legend("right",bty="n",cex=0.8,legend=c("PTA >= 90% (optimal)","PTA >= 75% (acceptable)"),lty=c(2,2),lwd=c(2,1.5),col=c("#2E7D32","#F57C00"))
  })

  # ── Repeated dose tab: PTA vs time for all loading doses ──────────────────
  output$repeated_dose_tab_content <- renderUI({
    if(is.null(rv$results)) return(
      div(style="padding:40px;text-align:center;color:#94a3b8;",
        icon("clock", style="font-size:2rem;margin-bottom:12px;display:block;"),
        p(style="font-size:15px;","Run the a priori simulation first."))
    )
    r <- rv$results
    if(!isTRUE(r$need_multiday)) return(
      div(
        style=paste0("padding:20px;border-radius:10px;background:linear-gradient(135deg,#f0fdf4,#e8f5e9);",
                     "border:1px solid #2E7D32;border-left:5px solid #2E7D32;margin-top:4px;"),
        HTML(paste0(
          "<b style='color:#2E7D32;font-size:15px;'>&#10003; All loading doses achieve PTA \u2265 90% on Day 1</b><br>",
          "<span style='color:#475569;font-size:13px;'>",
          "A repeated loading dose is not required. This tab is only active when at least one ",
          "loading dose fails to achieve PTA \u2265 90% on Day 1 (AUC 0\u201324h).",
          "</span>"
        ))
      )
    )
    tagList(
      div(style="margin-bottom:12px;",
        HTML(paste0(
          "<div style='padding:14px 18px;border-radius:10px;background:linear-gradient(135deg,#FFF3CD,#FFF8E1);",
          "border:1px solid #E65100;border-left:6px solid #B71C1C;font-size:13.5px;line-height:1.7;'>",
          "<b style='color:#B71C1C;font-size:15px;'>",
          "&#9888;&nbsp;Repeated loading dose required &mdash; ",r$best_load_day," consecutive day",
          if(r$best_load_day>1)"s" else "","</b><br>",
          "<span style='color:#555;'>",
          "No loading dose achieves PTA &ge;&nbsp;90% for the PK/PD target ",
          "(<b>AUC &ge;&nbsp;",r$target_auc,"&nbsp;mg&middot;h/L</b>) after a single administration (Day 1 alone).",
          "</span><br>",
          "A loading dose repeated over <b>",r$best_load_day,"&nbsp;consecutive days</b> is required.<br>",
          "<span style='color:#2E7D32;font-weight:700;'>",
          "&#10003;&nbsp;The plot and table below show PTA (%) for each loading dose on Day 1, Day 2 and Day 3.",
          "</span>",
          "</div>"
        ))
      ),
      card(
        card_header("Loading dose \u2013 PTA (%) vs Day for all doses"),
        plotOutput("pta_repeated_dose_plot", height=420)
      ),
      br(),
      card(
        card_header("Summary table \u2013 PTA (%) by dose and day"),
        tableOutput("pta_repeated_dose_table")
      )
    )
  })

  output$pta_repeated_dose_plot <- renderPlot({
    if(is.null(rv$results) || !isTRUE(rv$results$need_multiday)) return(NULL)
    r <- rv$results
    doses   <- r$load_df$Dose
    labels  <- r$load_df$DoseLabel
    pta_d1  <- r$load_pta_d1 * 100
    pta_d2  <- if(!is.null(r$load_pta_d2)) r$load_pta_d2 * 100 else rep(NA, length(doses))
    pta_d3  <- if(!is.null(r$load_pta_d3)) r$load_pta_d3 * 100 else rep(NA, length(doses))

    days    <- c(1, 2, 3)
    day_labels <- c("Day 1 (0\u201324h)", "Day 2 (24\u201348h)", "Day 3 (48\u201372h)")
    n_doses <- length(doses)

    # Qualitative palette — one distinct colour per loading dose (up to 8)
    pal_all <- c("#1565C0", "#2E7D32", "#AD1457", "#E65100",
                 "#6A1B9A", "#00695C", "#F9A825", "#37474F")
    pal <- pal_all[seq_len(min(n_doses, length(pal_all)))]
    if(n_doses > length(pal_all))
      pal <- c(pal, colorRampPalette(c("#90CAF9","#0D47A1"))(n_doses - length(pal_all)))

    par(mar = c(5, 5, 4, 14.5), bg = "white", xpd = FALSE)
    plot(days, rep(NA, 3), type = "n",
         xlim = c(0.7, 3.3), ylim = c(0, 105),
         xaxt = "n", yaxt = "n",
         xlab = "", ylab = "PTA (%)",
         main = paste0(
           "PTA vs Day \u2014 Loading dose candidates\n",
           "AUC target: \u2265 ", r$target_auc, " mg\u00b7h/L \u2014 ",
           r$best_load_day, " consecutive day",
           if(r$best_load_day > 1) "s" else "", " required"),
         font.main = 2, cex.main = 1.0, las = 1)
    axis(1, at = days, labels = day_labels, cex.axis = 0.92, font = 2)
    axis(2, at = seq(0, 100, by = 10), las = 1)
    mtext("Day", side = 1, line = 3.2, cex = 0.88)
    grid(col = "#ECEFF1", lty = 1)

    # Shaded zone 75–90%
    rect(0.7, 75, 3.3, 90, col = grDevices::adjustcolor("#FFF9C4", 0.45), border = NA)
    rect(0.7, 90, 3.3, 105, col = grDevices::adjustcolor("#E8F5E9", 0.35), border = NA)

    # Reference lines
    abline(h = 90, col = "#2E7D32", lty = 2, lwd = 2.2)
    abline(h = 75, col = "#F57C00", lty = 2, lwd = 1.8)
    text(3.28, 91.5, "90%", col = "#2E7D32", cex = 0.80, font = 2, xpd = TRUE)
    text(3.28, 76.5, "75%", col = "#F57C00", cex = 0.80, font = 2, xpd = TRUE)

    # Draw one line + points per dose
    for(i in seq_len(n_doses)) {
      pta_vec <- c(pta_d1[i], pta_d2[i], pta_d3[i])
      valid   <- !is.na(pta_vec)
      if(sum(valid) < 1) next
      lines(days[valid], pta_vec[valid], lwd = 2.5, col = pal[i])
      pt_fill <- ifelse(pta_vec[valid] >= 90, "#2E7D32",
                        ifelse(pta_vec[valid] >= 75, "#F57C00", "#C62828"))
      points(days[valid], pta_vec[valid],
             pch = 21, cex = 2.0,
             bg  = pt_fill,
             col = pal[i], lwd = 2.0)
    }

    # Legend — doses (line colour) + threshold key
    par(xpd = TRUE)
    legend("right", inset = c(-0.31, 0.1),
           bty = "o", box.col = "#E0E0E0", bg = "white",
           cex  = 0.74, title = expression(bold("Loading dose")),
           legend = rev(labels),
           lty    = rep(1, n_doses),
           lwd    = rep(2.5, n_doses),
           pch    = rep(21, n_doses),
           pt.cex = rep(1.6, n_doses),
           pt.bg  = rev(pal),
           col    = rev(pal))
    legend("right", inset = c(-0.31, -0.45),
           bty = "o", box.col = "#E0E0E0", bg = "white",
           cex = 0.74, title = expression(bold("Point = PTA level")),
           legend = c("PTA \u2265 90% (optimal)", "PTA \u2265 75% (acceptable)", "PTA < 75% (failure)"),
           pch    = c(21, 21, 21), pt.cex = c(1.6, 1.6, 1.6),
           pt.bg  = c("#2E7D32", "#F57C00", "#C62828"),
           col    = c("#555555", "#555555", "#555555"),
           lty    = NA, lwd    = NA)
    par(xpd = FALSE)
  })

  output$pta_repeated_dose_table <- renderTable({
    if(is.null(rv$results) || !isTRUE(rv$results$need_multiday)) return(NULL)
    r <- rv$results
    doses  <- r$load_df$Dose
    labels <- r$load_df$DoseLabel
    pta_d1 <- round(r$load_pta_d1 * 100, 1)
    pta_d2 <- if(!is.null(r$load_pta_d2)) round(r$load_pta_d2 * 100, 1) else rep(NA, length(doses))
    pta_d3 <- if(!is.null(r$load_pta_d3)) round(r$load_pta_d3 * 100, 1) else rep(NA, length(doses))

    # Best day to achieve >= 90%
    best_day <- sapply(seq_along(doses), function(i) {
      if(!is.na(pta_d1[i]) && pta_d1[i] >= 90) return("Day 1")
      if(!is.na(pta_d2[i]) && pta_d2[i] >= 90) return("Day 2")
      if(!is.na(pta_d3[i]) && pta_d3[i] >= 90) return("Day 3")
      if(!is.na(pta_d2[i]) && pta_d2[i] >= 75) return("Day 2 (suboptimal)")
      if(!is.na(pta_d3[i]) && pta_d3[i] >= 75) return("Day 3 (suboptimal)")
      "Failure"
    })

    data.frame(
      "Loading dose"      = labels,
      "PTA J1 (%)"        = pta_d1,
      "PTA J2 (%)"        = pta_d2,
      "PTA J3 (%)"        = pta_d3,
      "First day \u226590%" = best_day,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }, align = "c", striped = TRUE)

  output$dl_apriori<-downloadHandler(
    filename=function(){paste0("FlucoShiny_Apriori_",gsub("[^0-9A-Za-z_-]","",patient$ipp),"_",format(Sys.Date(),"%Y%m%d"),".html")},
    content=function(file){
      tryCatch({
        if(is.null(rv$results)||!isTRUE(patient$filled)){file.create(file);return()}
        logo_path<-Filter(file.exists,c("flucoshiny.png","Flucoshiny.png","FlucoShiny.png","www/flucoshiny.png","www/Flucoshiny.png","www/FlucoShiny.png"));logo_abs<-if(length(logo_path))normalizePath(logo_path[1])else NULL
        logo_html_val<-if(!is.null(logo_abs)&&file.exists(logo_abs)){b64<-knitr::image_uri(logo_abs);paste0('<img src="',b64,'" style="position:absolute;top:16px;right:22px;height:140px;max-width:400px;object-fit:contain;opacity:.92;" alt="logo">')}else ""
        tmpdir<-tempfile("fluco_ap_");dir.create(tmpdir,recursive=TRUE)
        rmd_path<-file.path(tmpdir,"report_apriori.Rmd");write_rmd_apriori(rmd_path,logo_abs)
        rmarkdown::render(rmd_path,output_file=file,params=list(patient=list(last_name=patient$last_name,first_name=patient$first_name,dob=format(patient$dob,"%d/%m/%Y"),service=patient$service,ipp=patient$ipp,sexe=patient$sexe,uf=patient$uf,hopital=patient$hopital,adresse=patient$adresse),res=rv$results,logo_path=logo_abs,logo_html=logo_html_val,nsim=NSIM,clinical=get_clinical(),pharmacologue=rv$pharmacologue,pk_model_ref=PK_MODEL_REF),envir=new.env(parent=globalenv()),quiet=FALSE)
      }, error=function(e){
        message("Report apriori error: ",e$message)
        writeLines(paste0("<html><body><h2>Report generation error</h2><p>",e$message,"</p></body></html>"),file)
      })
    })

  # ===================== A POSTERIORI =====================
  run_aposteriori<-function(){
    if(is.null(input$wt)||is.null(input$ckd)||is.null(input$tinf)||is.null(input$ap_dose_maint)||is.null(input$ap_n_conc)||is.null(input$ap_date_first))return()
    # Block if age < 18
    age_years<-calc_age_years(patient$dob)
    if(!is.na(age_years)&&age_years<AGE_ADULT){showModal(modalDialog(title="Model not applicable",HTML(paste0("<span class='extrapWarn'><b>Model restricted to adults (\u2265 ",AGE_ADULT," years). Patient aged ",round(age_years,1)," years.</b></span>")),footer=modalButton("OK"),easyClose=TRUE));return()}
    # Extrapolation warnings
    wt_warn<-(input$wt<WT_LOWER)||(input$wt>WT_UPPER); ckd_warn<-(input$ckd<CKD_LOWER)||(input$ckd>CKD_UPPER)
    if(wt_warn||ckd_warn){
      warn_msg<-c()
      if(wt_warn) warn_msg<-c(warn_msg,paste0("Weight outside model range (",WT_LOWER,"\u2013",WT_UPPER," kg)"))
      if(ckd_warn) warn_msg<-c(warn_msg,paste0("CKD-EPI outside model range (",CKD_LOWER,"\u2013",CKD_UPPER," mL/min/1.73m\u00b2)"))
      showModal(modalDialog(title="Extrapolation",HTML(paste0("<span class='extrapWarn'>",paste(warn_msg,collapse="<br>"),"</span><br><br>Results to be interpreted with caution.")),footer=actionButton("confirm_extrap_ap","Continue anyway",class="btn-warning"),easyClose=TRUE))
      observeEvent(input$confirm_extrap_ap,{removeModal();run_aposteriori_core()}, ignoreInit=TRUE, once=TRUE)
      return()
    }
    run_aposteriori_core()
  }
  run_aposteriori_core<-function(){
   tryCatch({
    # Hard validation - block before any mrgsolve call
    wt_val <- input$wt; ckd_val <- input$ckd; tinf_val <- input$tinf
    if(is.null(wt_val)||!is.finite(wt_val)||wt_val<=0) { showNotification("Invalid weight",type="error"); return() }
    if(is.null(ckd_val)||!is.finite(ckd_val)||ckd_val<=0) { showNotification("Invalid CKD-EPI",type="error"); return() }
    if(is.null(tinf_val)||!is.finite(tinf_val)||tinf_val<=0) { showNotification("Invalid infusion duration",type="error"); return() }
    cirrh<-ifelse(input$cirrh=="Yes",1,0)
    is_regular <- isTRUE(input$ap_regular)

    if (is_regular) {
      tau_h <- 24  # checkbox = toutes les 24h
      total_daily_dose <- input$ap_dose_maint * (24/tau_h); n_admin <- round(24/tau_h)
      heure_first <- if(!is.null(input$ap_heure_first)) input$ap_heure_first else 8
      heure_maint <- if(!is.null(input$ap_heure_maint)) input$ap_heure_maint else heure_first
      # Validate date of last administration
      if(is.null(input$ap_date_last) || is.na(input$ap_date_last)) {
        showModal(modalDialog(title="Missing date",
          HTML("<b>Please enter the date of the last maintenance dose administration.</b><br><br>This field is required, including after Monolix file import."),
          footer=modalButton("OK"), easyClose=TRUE))
        return()
      }
      if(input$ap_date_last < input$ap_date_maint) {
        showModal(modalDialog(title="Date error",
          "The last administration date must be on or after the 1st maintenance dose date.",
          footer=modalButton("OK"), easyClose=TRUE))
        return()
      }
      obs_full <- build_obs_from_dates(input$ap_n_conc, input, input$ap_date_first, heure_first=heure_first, tau=tau_h)
      obs <- data.frame(time=obs_full$time, conc=obs_full$conc)
      errs <- c(if(nrow(obs)==0)"At least one valid concentration.",if(any(obs$conc<=0))"Concentrations > 0.",if(input$ap_dose_maint<=0)"Dose > 0.",if(input$tinf<=0)"Infusion duration > 0.")
      if(length(errs)){showModal(modalDialog(title="Error",tags$ul(lapply(errs,tags$li)),footer=modalButton("OK"),easyClose=TRUE));return()}
      rv$ap_inputs <- list(total_daily_dose=total_daily_dose,n_admin=n_admin,tinf=input$tinf,tau=tau_h,
        date_first=format(input$ap_date_first,"%d/%m/%Y"),heure_first=heure_first,
        dose_first=input$ap_dose_first,
        date_maint=format(input$ap_date_maint,"%d/%m/%Y"),heure_maint=heure_maint,
        date_last=format(input$ap_date_last,"%d/%m/%Y"),
        dose_maint=input$ap_dose_maint,obs_full=obs_full)
      target_auc<-get_target_auc();mic_val<-if(!is.na(input$mic_exact)&&is.finite(input$mic_exact)&&input$mic_exact>0)input$mic_exact else if(isTRUE(input$glabrata))16 else 4
      withProgress(message="Bayesian estimation (saemix)...",value=0,{
        setProgress(0.1,detail="MAP estimation...")
        res<-tryCatch(estimate_aposteriori(obs_full,input$wt,input$ckd,cirrh,total_daily_dose,n_admin,input$tinf,
                                            date_first=format(input$ap_date_first,"%d/%m/%Y"),heure_first=heure_first,
                                            date_maint=format(input$ap_date_maint,"%d/%m/%Y"),heure_maint=heure_maint,
                                            date_last=format(input$ap_date_last,"%d/%m/%Y"),tau=tau_h,
                                            dose_first=input$ap_dose_first,dose_maint=input$ap_dose_maint),
                      error=function(e){message("AP error: ",e$message);rv$ap_error_msg<-e$message;NULL})
        setProgress(1);rv$ap_results<-res
      })
    } else {
      # Individual dose mode
      n_doses <- input$ap_n_doses; if(is.null(n_doses)||n_doses<1) return()
      # Get first dose date/time for obs computation
      d1 <- tryCatch(as.Date(input$ddose_1), error=function(e) NA)
      h1 <- as.numeric(input$hdose_1)
      if(is.na(d1)||is.na(h1)){showNotification("Date/time of 1st admin. required",type="error");return()}
      obs_full <- build_obs_from_dates(input$ap_n_conc, input, d1, heure_first=h1, tau=24)
      obs <- data.frame(time=obs_full$time, conc=obs_full$conc)
      if(nrow(obs)==0){showNotification("At least one valid concentration required",type="error");return()}
      nmtran <- build_nmtran_individual(n_doses, input, obs_full, input$wt, input$ckd, cirrh, input$tinf)
      if(is.null(nmtran)){showNotification("Check administrations",type="error");return()}
      # Get dose info for ap_inputs
      dose_rows_nm <- nmtran[nmtran$evid==1,]
      total_daily_dose <- if(nrow(dose_rows_nm)>1) dose_rows_nm$amt[nrow(dose_rows_nm)] else dose_rows_nm$amt[1]
      rv$ap_inputs<-list(total_daily_dose=total_daily_dose,n_admin=1,tinf=input$tinf,tau=24,date_first=format(d1,"%d/%m/%Y"),heure_first=h1,dose_first=dose_rows_nm$amt[1],date_maint=format(d1,"%d/%m/%Y"),heure_maint=h1,date_last=format(d1+max(dose_rows_nm$time)/24,"%d/%m/%Y"),dose_maint=total_daily_dose,obs_full=obs_full)
      target_auc<-get_target_auc();mic_val<-if(!is.na(input$mic_exact)&&is.finite(input$mic_exact)&&input$mic_exact>0)input$mic_exact else if(isTRUE(input$glabrata))16 else 4
      withProgress(message="Bayesian estimation (saemix)...",value=0,{
        setProgress(0.1,detail="MAP estimation...")
        # Direct saemix call with pre-built NM-TRAN
        res <- tryCatch({
          wt_c <- max(WT_LOWER,min(input$wt,WT_UPPER)); ckd_c <- max(CKD_LOWER,min(input$ckd,CKD_UPPER)); tinf_c <- max(0.05,input$tinf)
          CLtyp <- predict_CL_typ(ckd_c, cirrh); Vtyp <- predict_V_typ(wt_c)
          message("=== NM-TRAN individual ==="); print(nmtran)
          est <- fluco_saemix_estimate(nmtran)
          CL_map <- tryCatch(fluco_get_CL(est),error=function(e)CLtyp)
          V_map <- tryCatch(fluco_get_V(est),error=function(e)Vtyp)
          etas <- tryCatch(fluco_get_etas(est),error=function(e)c(ETA1=0,ETA2=0))
          eta1 <- as.numeric(etas[1]); eta2 <- as.numeric(etas[2])
          if(!is.finite(eta1)) eta1<-0; if(!is.finite(eta2)) eta2<-0
          CL_cm_i <- if(!is.null(est$cond_mean)) as.numeric(est$cond_mean["CL"]) else CL_map
          V_cm_i  <- if(!is.null(est$cond_mean)) as.numeric(est$cond_mean["V"])  else V_map
          if(!is.finite(CL_cm_i) || CL_cm_i<=0) CL_cm_i <- CL_map
          if(!is.finite(V_cm_i)  || V_cm_i <=0) V_cm_i  <- V_map
          # ---- Posterior sampling — APPROCHE LEVOSHINY (saemix nsamp=1000) ----
          cond_s <- fluco_get_cond_samples(est)
          if(!is.null(cond_s) && length(cond_s$CL)>0 && length(cond_s$V)>0) {
            set.seed(12345)
            idx <- sample.int(length(cond_s$CL), size=NPOST_SAMP, replace=TRUE)
            cl_v <- pmax(as.numeric(cond_s$CL[idx]),0.001)
            v_v  <- pmax(as.numeric(cond_s$V[idx]), 0.1)
            post <- data.frame(etaCL=log(cl_v/CLtyp), etaV=log(v_v/Vtyp), CL=cl_v, V=v_v)
          } else {
            samp_eta <- rmvnorm_chol(NPOST_SAMP, c(eta1,eta2), OMEGA_MAT)
            post <- data.frame(etaCL=samp_eta[,1], etaV=samp_eta[,2])
            post$CL <- pmax(CLtyp*exp(post$etaCL),0.001); post$V <- pmax(Vtyp*exp(post$etaV),0.1)
          }
          dose <- total_daily_dose; tau_v <- 24
          dose_times_auc <- seq(0, by=tau_v, length.out=14)
          dose_amounts_auc <- rep(dose, length(dose_times_auc))
          post$AUC24 <- mapply(function(cl, v) {
            auc_model_interval((length(dose_times_auc)-1)*tau_v, length(dose_times_auc)*tau_v,
                               dose_times_auc, dose_amounts_auc, cl, v, tinf_c)
          }, post$CL, post$V)
          tgrid <- seq(0,tau_v,length.out=POST_PROFILE_GRID)
          prof <- tryCatch({
            mod_patient <- mod_fluco %>% param(ETA1=eta1,ETA2=eta2,CKD=ckd_c,WT=wt_c,CIRRH=cirrh) %>% zero_re()
            sim_prof <- mod_patient %>% ev(amt=dose,rate=dose/tinf_c,ii=tau_v,addl=13,cmt=1) %>% mrgsim(end=14*tau_v,delta=tau_v/POST_PROFILE_GRID) %>% as.data.frame()
            ss_start <- 13*tau_v; sim_ss <- sim_prof[sim_prof$time>=ss_start&sim_prof$time<=(ss_start+tau_v),]
            dup <- !duplicated(sim_ss$time-ss_start)
            pm <- approx((sim_ss$time-ss_start)[dup],sim_ss$CP[dup],xout=tgrid,rule=2)$y
            if(any(is.na(pm))) pm[is.na(pm)]<-0
            data.frame(time=tgrid,map=pm,condmean=pm,p2.5=pm*0.8,p50=pm,p97.5=pm*1.2)
          }, error=function(e) data.frame(time=tgrid,map=rep(0,length(tgrid)),condmean=rep(0,length(tgrid)),p2.5=rep(0,length(tgrid)),p50=rep(0,length(tgrid)),p97.5=rep(0,length(tgrid))))
          mod_patient <- tryCatch(mod_fluco %>% param(ETA1=eta1,ETA2=eta2,CKD=ckd_c,WT=wt_c,CIRRH=cirrh) %>% zero_re(), error=function(e) NULL)
          obs_compat <- data.frame(time=obs_full$time, conc=obs_full$conc)
          list(CLtyp=CLtyp,Vtyp=Vtyp,dose=dose,tau=tau_v,tinf=tinf_c,obs=obs_compat,
               map=c(etaCL=eta1,etaV=eta2),map_psi=c(CL=CL_map,V=V_map),cond_mean_psi=c(CL=CL_cm_i,V=V_cm_i),
               posterior=post,profile=prof,convergence=0,objective=0,nmtran=nmtran,saemix_est=est,mod_patient=mod_patient)
        }, error=function(e){message("AP error: ",e$message);rv$ap_error_msg<-e$message;NULL})
        setProgress(1);rv$ap_results<-res
      })
    }
    if(is.null(rv$ap_results)){err_txt<-if(!is.null(rv$ap_error_msg))rv$ap_error_msg else "Erreur inconnue";showModal(modalDialog(title="Estimation failed",HTML(paste0("Please check your data.<br><br><small style='color:#999;'>Detail: ",err_txt,"</small>")),footer=modalButton("OK"),easyClose=TRUE))}
   }, error=function(e){
     message("FATAL run_aposteriori_core: ", e$message)
     showModal(modalDialog(title="Unexpected error",
       HTML(paste0("<b>Estimation failed.</b><br><br>Detail: <code>",gsub("<","&lt;",e$message),"</code><br><br>Please check weight, CKD-EPI and entered data.")),
       footer=modalButton("OK"),easyClose=TRUE))
   })
  }

  # Stored values for dynamic inputs (survives re-render)
  rv_fields <- reactiveValues(doses=list(), concs=list())

  # Dynamic dose administration inputs with X button
  output$ap_dose_inputs<-renderUI({
    n<-input$ap_n_doses;if(is.null(n)||!is.finite(n)||n<1)return(NULL)
    tagList(lapply(seq_len(n),function(i){
      # Read from stored values or from existing input or defaults
      stored <- isolate(rv_fields$doses[[as.character(i)]])
      cur_d <- if(!is.null(stored)) stored$d else isolate(input[[paste0("ddose_",i)]])
      cur_h <- if(!is.null(stored)) stored$h else isolate(input[[paste0("hdose_",i)]])
      cur_a <- if(!is.null(stored)) stored$a else isolate(input[[paste0("adose_",i)]])
      dd <- if(!is.null(cur_d)) cur_d else Sys.Date()-7+i-1
      hh <- if(!is.null(cur_h)&&is.finite(cur_h)) cur_h else 8
      aa <- if(!is.null(cur_a)&&is.finite(cur_a)) cur_a else if(i==1) 800 else 400
      tagList(fluidRow(
        column(4, dateInput(paste0("ddose_",i),paste0("Admin ",i),value=dd,format="yyyy-mm-dd")),
        column(3, numericInput(paste0("hdose_",i),"Administration time",value=hh,min=0,max=23.99,step=0.5)),
        column(4, numericInput(paste0("adose_",i),"Dose (mg)",value=aa,min=0,step=25)),
        column(1, tags$label(HTML("&nbsp;"),class="form-label d-block"), actionButton(paste0("del_dose_",i),"\u2716",class="btn-outline-danger btn-sm",style="padding:2px 8px;font-size:11px;"))
      ))
    }))
  })

  # Delete dose handler
  observe({
    n <- input$ap_n_doses; if(is.null(n)||n<1) return()
    lapply(seq_len(n), function(i){
      observeEvent(input[[paste0("del_dose_",i)]], {
        cur_n <- input$ap_n_doses; if(is.null(cur_n)||cur_n<=1) return()
        # Save all current values except deleted one
        vals <- list()
        new_idx <- 1
        for(j in seq_len(cur_n)){
          if(j==i) next
          vals[[as.character(new_idx)]] <- list(
            d=input[[paste0("ddose_",j)]],
            h=as.numeric(input[[paste0("hdose_",j)]]),
            a=as.numeric(input[[paste0("adose_",j)]])
          )
          new_idx <- new_idx + 1
        }
        rv_fields$doses <- vals
        updateNumericInput(session, "ap_n_doses", value=cur_n-1)
      }, ignoreInit=TRUE, once=TRUE)
    })
  })

  # Build NM-TRAN from individual doses (non-regular mode)
  build_nmtran_individual <- function(n_doses, input_vals, obs_full, wt, ckd, cirrh, tinf) {
    wt <- max(WT_LOWER, min(wt, WT_UPPER)); ckd <- max(CKD_LOWER, min(ckd, CKD_UPPER)); tinf <- max(0.05, tinf)
    dose_list <- lapply(seq_len(n_doses), function(i) {
      d <- tryCatch(as.Date(input_vals[[paste0("ddose_",i)]]), error=function(e) NA)
      h <- as.numeric(input_vals[[paste0("hdose_",i)]])
      amt <- as.numeric(input_vals[[paste0("adose_",i)]])
      if (is.na(d)||is.na(h)||is.na(amt)||amt<=0) return(NULL)
      list(dt=as.POSIXct(paste(d,sprintf("%02d:%02d:00",floor(h),round((h%%1)*60))),format="%Y-%m-%d %H:%M:%S"), amt=amt)
    })
    dose_list <- Filter(Negate(is.null), dose_list)
    if (length(dose_list)==0) return(NULL)
    dts <- sapply(dose_list, function(x) x$dt)
    ord <- order(dts); dose_list <- dose_list[ord]
    t0 <- dose_list[[1]]$dt
    dose_rows <- do.call(rbind, lapply(dose_list, function(x) {
      data.frame(ID=1, time=as.numeric(difftime(x$dt,t0,units="hours")), evid=1, cmt=1, amt=x$amt, rate=x$amt/tinf, mdv=1, DV=0, CKD=ckd, WT=wt, CIRRH=cirrh, stringsAsFactors=FALSE)
    }))
    if (is.null(obs_full)||nrow(obs_full)==0) return(NULL)
    obs_rows <- data.frame(ID=1, time=obs_full$abs_h, evid=0, cmt=1, amt=0, rate=0, mdv=0, DV=obs_full$conc, CKD=ckd, WT=wt, CIRRH=cirrh, stringsAsFactors=FALSE)
    nmtran <- rbind(dose_rows, obs_rows)
    nmtran <- nmtran[order(nmtran$time, -nmtran$evid),]
    rownames(nmtran) <- NULL
    nmtran
  }

  output$ap_conc_inputs<-renderUI({n<-input$ap_n_conc;if(is.null(n)||!is.finite(n)||n<1)return(NULL);tagList(lapply(seq_len(n),function(i){
    stored <- isolate(rv_fields$concs[[as.character(i)]])
    cur_d <- if(!is.null(stored)) stored$d else isolate(input[[paste0("dconc_",i)]])
    cur_h <- if(!is.null(stored)) stored$h else isolate(input[[paste0("hconc_",i)]])
    cur_c <- if(!is.null(stored)) stored$c else isolate(input[[paste0("cconc_",i)]])
    dd <- if(!is.null(cur_d)) cur_d else Sys.Date()
    hh <- if(!is.null(cur_h)&&is.finite(cur_h)) cur_h else if(i==1) 8 else NA
    cc <- if(!is.null(cur_c)&&is.finite(cur_c)) cur_c else if(i==1) 10 else NA
    cur_ckd <- if(!is.null(stored$ckd)) stored$ckd else isolate(input[[paste0("ckdconc_",i)]])
    ckd_val <- if(!is.null(cur_ckd) && !is.na(suppressWarnings(as.numeric(cur_ckd)))) as.numeric(cur_ckd) else isolate(input$ckd)
    fluidRow(
      column(3, dateInput(paste0("dconc_",i),paste0("Sample ",i," - Date"),value=dd,format="yyyy-mm-dd")),
      column(2, numericInput(paste0("hconc_",i),"Time (h)",value=hh,min=0,max=23.99,step=0.5)),
      column(3, numericInput(paste0("cconc_",i),"Concentration (mg/L)",value=cc,min=0.001,step=0.1)),
      column(3, numericInput(paste0("ckdconc_",i),"CKD-EPI at sampling (mL/min/1.73m²)",value=ckd_val,min=0,max=300,step=1)),
      column(1, tags$label(HTML("&nbsp;"),class="form-label d-block"), actionButton(paste0("del_conc_",i),"\u2716",class="btn-outline-danger btn-sm",style="padding:2px 8px;font-size:11px;"))
    )
  }))})

  # Delete concentration handler
  observe({
    n <- input$ap_n_conc; if(is.null(n)||n<1) return()
    lapply(seq_len(n), function(i){
      observeEvent(input[[paste0("del_conc_",i)]], {
        cur_n <- input$ap_n_conc; if(is.null(cur_n)||cur_n<=1) return()
        vals <- list(); new_idx <- 1
        for(j in seq_len(cur_n)){
          if(j==i) next
          vals[[as.character(new_idx)]] <- list(
            d=input[[paste0("dconc_",j)]],
            h=as.numeric(input[[paste0("hconc_",j)]]),
            c=as.numeric(input[[paste0("cconc_",j)]])
          )
          new_idx <- new_idx + 1
        }
        rv_fields$concs <- vals
        updateNumericInput(session, "ap_n_conc", value=cur_n-1)
      }, ignoreInit=TRUE, once=TRUE)
    })
  })

  output$ap_ss_warning<-renderUI({
    if(is.null(input$ap_date_first)||is.null(input$ap_date_last))return(NULL)
    days_treatment<-as.numeric(difftime(input$ap_date_last,input$ap_date_first,units="days"))
    if(days_treatment<5) div(style="background:#FFF3E0;border-left:4px solid #F57C00;padding:8px 12px;border-radius:4px;margin-bottom:8px;",HTML(paste0("<b>Warning:</b> ",round(days_treatment,1)," days of treatment. Steady-state not guaranteed (fluconazole t1/2 ~30h, SS after ~6 days).")))
  })
  output$ap_summary<-renderUI({tryCatch({if(is.null(rv$ap_results))return(NULL);r<-rv$ap_results;q_cl<-summ_quant(r$posterior$CL);q_v<-summ_quant(r$posterior$V)
    wt_ext <- (input$wt<WT_LOWER)||(input$wt>WT_UPPER); ckd_ext <- (input$ckd<CKD_LOWER)||(input$ckd>CKD_UPPER)
    extrap_ui <- if(wt_ext||ckd_ext) div(style="background:#FFF3E0;border-left:4px solid #F57C00;padding:8px 12px;border-radius:4px;margin-bottom:8px;", HTML(paste0(
      "<b>\u26a0 Extrapolation:</b> ",
      if(wt_ext) paste0("Weight outside model range (",WT_LOWER,"\u2013",WT_UPPER," kg). ") else "",
      if(ckd_ext) paste0("CKD-EPI outside model range (",CKD_LOWER,"\u2013",CKD_UPPER," mL/min/1.73m\u00b2). ") else "",
      "Results to be interpreted with caution."))) else NULL
    tagList(extrap_ui,
      # ── AUC MAP: primary result, highlighted prominently ──
      uiOutput("ap_auc_summary_box"),
      # ── PK parameters card ──
      card(style="background-color:#E8F0FE;",card_header("Individual a posteriori PK parameters"),
        p(HTML(paste0("<b>Population CL:</b> ",round(r$CLtyp,3)," L/h | <b>MAP CL:</b> ",round(r$map_psi["CL"],3)," L/h"))),
        p(HTML(paste0("<b>Population V:</b> ",round(r$Vtyp,1)," L | <b>MAP V:</b> ",round(r$map_psi["V"],1)," L"))),
        p(HTML(paste0("<b>Posterior CL [95% CI]:</b> ",round(q_cl["median"],3)," [",round(q_cl["p2.5"],3),"-",round(q_cl["p97.5"],3),"] L/h"))),
        p(HTML(paste0("<b>Posterior Vd [95% CI]:</b> ",round(q_v["median"],1)," [",round(q_v["p2.5"],1),"-",round(q_v["p97.5"],1),"] L")))
      )
    )
  },error=function(e)NULL)})
  output$ap_param_table<-renderTable({tryCatch({if(is.null(rv$ap_results))return(NULL);r<-rv$ap_results;q_cl<-summ_quant(r$posterior$CL);q_v<-summ_quant(r$posterior$V);data.frame(Parameter=c("CL (L/h)","V (L)"),Mean=round(c(q_cl["mean"],q_v["mean"]),3),Median=round(c(q_cl["median"],q_v["median"]),3),"P2.5"=round(c(q_cl["p2.5"],q_v["p2.5"]),3),"P97.5"=round(c(q_cl["p97.5"],q_v["p97.5"]),3),check.names=FALSE)},error=function(e)NULL)},align="c")

  # Prominent AUC MAP summary box - shown above PK params
  output$ap_auc_summary_box <- renderUI({
    tryCatch({
      if(is.null(rv$ap_results)||is.null(rv$ap_inputs)) return(NULL)
      r <- rv$ap_results; apr <- rv$ap_inputs
      obs_f <- apr$obs_full
      if(is.null(obs_f)||nrow(obs_f)==0) return(NULL)
      target_auc <- get_target_auc()
      tau_h <- if(!is.null(apr$tau)) apr$tau else 24
      d1 <- as.Date(apr$date_first, "%d/%m/%Y")
      dose_ld <- if(!is.null(apr$dose_first)&&apr$dose_first>0) apr$dose_first else apr$total_daily_dose
      dose_md <- if(!is.null(apr$dose_maint)&&apr$dose_maint>0) apr$dose_maint else apr$total_daily_dose
      days_total <- as.numeric(difftime(as.Date(apr$date_last,"%d/%m/%Y"),d1,units="days"))
      dose_times <- 0; dose_amounts <- dose_ld
      if(days_total>0){mt<-seq(tau_h,days_total*24,by=tau_h);dose_times<-c(dose_times,mt);dose_amounts<-c(dose_amounts,rep(dose_md,length(mt)))}

      # ---- MIC used (same logic as get_target_auc) ---------------------------
      # Fluconazole PK/PD target: AUC/MIC >= 100 (constant)
      # MIC = mic_exact if entered, else 16 (glabrata), else 4 mg/L
      # (C. albicans / C. parapsilosis / C. tropicalis - EUCAST clinical breakpoint)
      mic_used <- if(!is.na(input$mic_exact) && is.finite(input$mic_exact) && input$mic_exact > 0) {
        as.numeric(input$mic_exact)
      } else if(isTRUE(input$glabrata)) {
        16
      } else {
        4
      }
      mic_is_documented <- !is.na(input$mic_exact) && is.finite(input$mic_exact) && input$mic_exact > 0
      mic_is_glabrata   <- isTRUE(input$glabrata) && !mic_is_documented

      mic_src_lbl <- if(mic_is_documented) {
        paste0("Documented MIC: <b>", mic_used, " mg/L</b>")
      } else if(mic_is_glabrata) {
        paste0("MIC = <b>16 mg/L</b> &mdash; <i>Clinical breakpoint <span style='font-style:italic'>C. glabrata</span></i>")
      } else {
        paste0("MIC = <b>4 mg/L</b> &mdash; <i>Clinical breakpoint <span style='font-style:italic'>C. albicans</span>, <span style='font-style:italic'>C. parapsilosis</span>, <span style='font-style:italic'>C. tropicalis</span></i>")
      }
      target_range_lbl <- "AUC/MIC \u2265 100"
      ratio_target_min <- 100  # PK/PD target

      auc_boxes <- lapply(seq_len(nrow(obs_f)), function(oo){
        ah <- obs_f$abs_h[oo]
        doses_before <- dose_times[dose_times < ah]
        if(length(doses_before)==0) doses_before <- dose_times[dose_times<=ah]
        t_start <- if(length(doses_before)>0) max(doses_before) else 0
        # AUC definition:
        #  - Next dose after obs → AUCtau [t_last → t_next] (complete interval)
        #  - No next dose after obs → AUC [t_last → t_obs] (truncated at sample time)
        doses_after_st <- dose_times[dose_times > t_start]
        has_next_st <- length(doses_after_st) > 0 && min(doses_after_st) > ah
        t_end_tau <- if(has_next_st) min(doses_after_st) else t_start + 24
        t_post <- round(ah - t_start, 1)
        auc_map <- auc_model_interval(t_start, t_end_tau, dose_times, dose_amounts, r$map_psi["CL"], r$map_psi["V"], r$tinf)
        # 95% CI from posterior
        post <- r$posterior; n_s <- min(500,nrow(post)); idx_s <- sample(nrow(post),n_s)
        auc_post <- sapply(idx_s, function(s) auc_model_interval(t_start, t_end_tau, dose_times, dose_amounts, post$CL[s], post$V[s], r$tinf))
        qa <- quantile(auc_post, c(0.025,0.975))
        # Days since first dose (for "dose reduction possible" threshold)
        days_since_start <- ah / 24
        at_ss <- days_since_start >= 6  # steady-state reached after ~6 days (5\u00d7t\u00bd)

        # Dynamic thresholds based on MIC (all derived from target_auc = 100\u00d7MIC)
        thr_target    <- target_auc          # 100\u00d7MIC \u2014 minimum target
        thr_reduce    <- 1.5 * target_auc   # 150\u00d7MIC \u2014 reduction possible if SS reached
        thr_over      <- 2.0 * target_auc   # 200\u00d7MIC \u2014 overexposure

        # 4-tier classification
        tier <- if (auc_map < thr_target) {
          "not_achieved"
        } else if (auc_map >= thr_over) {
          "overexposure"
        } else if (auc_map >= thr_reduce && at_ss) {
          "reduce"
        } else {
          "achieved"
        }

        bg  <- switch(tier, not_achieved="#FFF3E0", overexposure="#FDECEA", reduce="#E3F2FD", achieved="#E8F5E9")
        brd <- switch(tier, not_achieved="#F57C00", overexposure="#C62828", reduce="#1565C0", achieved="#2E7D32")
        ico <- switch(tier, not_achieved="\u26a0\ufe0f", overexposure="\ud83d\udd34", reduce="\ud83d\udd35", achieved="\u2705")

        # AUC/MIC ratio (point estimate on AUC MAP + 95% CI from posterior)
        auc_mic    <- round(auc_map / mic_used, 0)
        auc_mic_lo <- round(as.numeric(qa[1]) / mic_used, 0)
        auc_mic_hi <- round(as.numeric(qa[2]) / mic_used, 0)

        lbl <- switch(tier,
          not_achieved = paste0(ico," PK/PD target NOT ACHIEVED \u2014 Dosage adjustment required"),
          overexposure = paste0(ico," OVEREXPOSURE (AUC/MIC \u2265 200) \u2014 Dose reduction recommended"),
          reduce       = paste0(ico," Target achieved \u2014 Dose reduction possible (AUC/MIC \u2265 150, steady-state reached)"),
          achieved     = paste0(ico," PK/PD target ACHIEVED")
        )

        div(style=paste0("background:",bg,";border-left:6px solid ",brd,";border-radius:12px;padding:18px 22px;margin-bottom:14px;"),
          # Sample info
          div(style="font-size:13px;color:#555;margin-bottom:8px;",
              paste0("Sample ", oo, " \u2014 T post-dose: ", t_post, "h \u2014 Observed concentration: ", round(obs_f$conc[oo],2), " mg/L")),

          # === MAIN RESULT: PK/PD ===
          div(style=paste0("background:rgba(255,255,255,0.75);border-radius:10px;padding:14px 16px;margin-bottom:10px;border:2px solid ",brd,";"),
            div(style="font-size:12px;font-weight:700;color:#555;letter-spacing:0.5px;text-transform:uppercase;margin-bottom:6px;",
                "\U0001F3AF Primary PK/PD marker"),
            div(style=paste0("font-size:28px;font-weight:900;color:",brd,";line-height:1.1;margin-bottom:4px;"),
                paste0("AUC/MIC = ", auc_mic)),
            div(style="font-size:14px;font-weight:700;color:#444;margin-bottom:4px;",
                paste0("Target: ", target_range_lbl)),
            div(style="font-size:13px;color:#555;margin-bottom:4px;", HTML(mic_src_lbl)),
            div(style=paste0("font-size:16px;font-weight:800;color:",brd,";margin-top:8px;padding:6px 10px;background:rgba(0,0,0,0.04);border-radius:6px;"), lbl)
          ),

          # AUC MAP (secondary)
          div(style="font-size:13px;color:#444;margin-bottom:2px;",
              HTML(paste0("<b>AUC MAP = ", round(auc_map, 1), " mg\u00b7h/L</b>",
                " &nbsp;|&nbsp; 95% CI: [", round(qa[1],1), "\u2013", round(qa[2],1), "] mg\u00b7h/L",
                " &nbsp;|&nbsp; AUC/MIC 95% CI: [", auc_mic_lo, "\u2013", auc_mic_hi, "]"))),
          div(style="font-size:12px;color:#666;margin-bottom:8px;",
              if(has_next_st) paste0("AUC 0\u2013", round(t_end_tau-t_start,1), "h [last dose \u2192 next dose]")
              else paste0("AUC 0\u201324h [last dose \u2192 last dose + 24h]"))
        )
      })
      tagList(auc_boxes)
    }, error=function(e) NULL)
  })

  output$ap_auc_interval_table<-renderTable({
    tryCatch({
    if(is.null(rv$ap_results)||is.null(rv$ap_inputs))return(NULL)
    r <- rv$ap_results; apr <- rv$ap_inputs
    obs_f <- apr$obs_full
    if (is.null(obs_f) || nrow(obs_f) == 0) return(data.frame("Info"="No observations available", check.names=FALSE))
    tau_h <- if(!is.null(apr$tau)) apr$tau else 24
    d1 <- as.Date(apr$date_first, "%d/%m/%Y")
    dose_ld <- if (!is.null(apr$dose_first) && apr$dose_first > 0) apr$dose_first else apr$total_daily_dose
    dose_md <- if (!is.null(apr$dose_maint) && apr$dose_maint > 0) apr$dose_maint else apr$total_daily_dose
    days_total <- as.numeric(difftime(as.Date(apr$date_last, "%d/%m/%Y"), d1, units="days"))
    dose_times <- 0; dose_amounts <- dose_ld
    if (days_total > 0) { mt <- seq(tau_h, days_total*24, by=tau_h); dose_times <- c(dose_times, mt); dose_amounts <- c(dose_amounts, rep(dose_md, length(mt))) }
    post <- r$posterior
    auc_rows <- lapply(seq_len(nrow(obs_f)), function(oo) {
      ah <- obs_f$abs_h[oo]  # absolute time of observation (h from t0)
      # Last dose time STRICTLY before observation (not coincident)
      doses_before <- dose_times[dose_times < ah]
      if (length(doses_before) == 0) doses_before <- dose_times[dose_times <= ah]
      t_start <- if(length(doses_before)>0) max(doses_before) else 0
      # AUC definition:
      #  - Next dose exists AND obs before it → AUCtau [t_last → t_next] (complete)
      #  - No next dose after obs → AUC [t_last → t_obs] (truncated at sample time)
      doses_after_start <- dose_times[dose_times > t_start]
      has_next <- length(doses_after_start) > 0 && min(doses_after_start) > ah
      t_end <- if(has_next) min(doses_after_start) else t_start + 24
      t_post <- round(ah - t_start, 2)
      date_obs <- format(d1 + ah/24, "%d/%m/%Y %H:%M")
      date_admin <- format(d1 + t_start/24, "%d/%m/%Y")
      # AUC[t_last → t_obs] with MAP CL/V (conditional mode)
      auc_mode <- auc_model_interval(t_start, t_end, dose_times, dose_amounts, r$map_psi["CL"], r$map_psi["V"], r$tinf)
      # Posterior distribution AUC[t_last → t_obs]
      n_s <- min(500, nrow(post)); idx_s <- sample(nrow(post), n_s)
      auc_post <- sapply(idx_s, function(s) auc_model_interval(t_start, t_end, dose_times, dose_amounts, post$CL[s], post$V[s], r$tinf))
      qa <- quantile(auc_post, c(0.025,0.5,0.975))
      data.frame(
        Administration=date_admin,
        EndInterval=format(d1 + t_end/24, "%d/%m/%Y"),
        T_post_dose=t_post,
        C_obs=obs_f$conc[oo],
        AUC_mode=round(auc_mode,1),
        IC95=paste0(round(as.numeric(qa[1]),1),"\u2013",round(as.numeric(qa[3]),1)),
        check.names=FALSE, stringsAsFactors=FALSE)
    })
    auc_df <- do.call(rbind, auc_rows)
    names(auc_df) <- c("Last dose","End of interval","T post-dose (h)","C obs (mg/L)",
                       "AUC MAP (mg\u00b7h/L)","95% CI")
    auc_df
    }, error=function(e){message("AUC table error:",e$message);NULL})
  }, align="c")
  output$ap_plot_profile<-renderPlot({
    tryCatch({
    if(is.null(rv$ap_results)||is.null(rv$ap_inputs))return(NULL)
    r <- rv$ap_results; apr <- rv$ap_inputs; post <- r$posterior
    nmtran <- r$nmtran; if(is.null(nmtran)) return(NULL)
    hist_doses <- nmtran[nmtran$evid==1,]; obs_rows <- nmtran[nmtran$evid==0,]
    if(nrow(hist_doses)==0) return(NULL)
    t_last <- max(hist_doses$time); tau_h <- if(!is.null(apr$tau)) apr$tau else 24
    # Ensure time grid covers all observations too
    obs_max_t <- if(nrow(obs_rows)>0) max(as.numeric(obs_rows$time), na.rm=TRUE) else 0
    t_end <- max(t_last + tau_h + 12, obs_max_t + 12)

    # MAP digital twin (like prof: mod_patient %>% zero_re())
    eta1 <- r$map["etaCL"]; eta2 <- r$map["etaV"]
    if(!is.finite(eta1)) eta1<-0; if(!is.finite(eta2)) eta2<-0
    mod_pat <- tryCatch(mod_fluco %>%
      param(ETA1=eta1, ETA2=eta2, CKD=max(5,hist_doses$CKD[1]), WT=max(20,hist_doses$WT[1]), CIRRH=hist_doses$CIRRH[1]) %>%
      zero_re(), error=function(e){message("mod_pat error:",e$message);NULL})
    if(is.null(mod_pat)) return(NULL)
    ev_df <- data.frame(ID=1, time=hist_doses$time, amt=hist_doses$amt, rate=hist_doses$rate, evid=1, cmt=1)

    # MAP simulation
    sim_map <- tryCatch(mod_pat %>% ev(as.ev(ev_df)) %>% mrgsim(end=t_end, delta=0.5) %>% as.data.frame(), error=function(e){message("Profile sim error:",e$message);NULL})
    if(is.null(sim_map)||nrow(sim_map)==0) return(NULL)
    tgrid_sim <- sim_map$time; map_conc <- sim_map$CP

    # Posterior CI - use common time grid with interpolation for robustness
    tgrid_common <- seq(0, t_end, by=0.5)
    n_sub <- min(150, nrow(post)); idx_sub <- sample(nrow(post), n_sub)
    post_mat <- matrix(NA, n_sub, length(tgrid_common))
    for(s in seq_len(n_sub)){
      mod_s <- tryCatch(mod_fluco %>%
        param(ETA1=post$etaCL[idx_sub[s]], ETA2=post$etaV[idx_sub[s]],
              CKD=max(5,hist_doses$CKD[1]), WT=max(20,hist_doses$WT[1]), CIRRH=hist_doses$CIRRH[1]) %>%
        zero_re(), error=function(e) NULL)
      if(is.null(mod_s)) next
      sim_s <- tryCatch(mod_s %>% ev(as.ev(ev_df)) %>% mrgsim(end=t_end, delta=0.5) %>% as.data.frame(), error=function(e) NULL)
      if(!is.null(sim_s) && nrow(sim_s)>1) {
        dup <- !duplicated(sim_s$time)
        post_mat[s,] <- approx(sim_s$time[dup], sim_s$CP[dup], xout=tgrid_common, rule=2)$y
      }
    }
    keep <- apply(post_mat, 1, function(x) !any(is.na(x))); post_mat <- post_mat[keep,,drop=FALSE]
    has_ci <- nrow(post_mat) >= 2
    if(has_ci){
      p50 <- apply(post_mat,2,median); p025 <- apply(post_mat,2,quantile,0.025); p975 <- apply(post_mat,2,quantile,0.975)
    }

    # Interpolate MAP on common grid (deduplicate to avoid warnings)
    dup_map <- !duplicated(tgrid_sim)
    map_interp <- approx(tgrid_sim[dup_map], map_conc[dup_map], xout=tgrid_common, rule=2)$y
    tgrid_d <- tgrid_common / 24

    # Observations - from NM-TRAN or from ap_inputs$obs_full
    obs_t <- numeric(0); obs_c <- numeric(0)
    if(nrow(obs_rows)>0) {
      obs_t <- as.numeric(obs_rows$time) / 24
      obs_c <- as.numeric(obs_rows$DV)
    } else if(!is.null(apr$obs_full) && nrow(apr$obs_full)>0) {
      obs_t <- apr$obs_full$abs_h / 24
      obs_c <- apr$obs_full$conc
    }
    # Remove any NA
    valid <- is.finite(obs_t) & is.finite(obs_c) & obs_c > 0
    obs_t <- obs_t[valid]; obs_c <- obs_c[valid]
    message("=== Profile plot: ", nrow(hist_doses), " doses, ", length(obs_c), " obs, CI samples=", sum(keep), " ===")

    all_vals <- c(if(has_ci) p975 else NULL, map_interp, obs_c)
    ylim_max <- if(length(all_vals)>0) max(all_vals, na.rm=TRUE) * 1.2 else 50

    d1_date <- tryCatch(as.Date(apr$date_first, "%d/%m/%Y"), error=function(e) Sys.Date()-7)
    day_ticks <- seq(0, max(tgrid_d), by=1); date_labels <- format(d1_date + day_ticks, "%d/%m")

    par(mar=c(5,4.5,3,10), bg="white", xpd=FALSE)
    plot(tgrid_d, map_interp, type="n", xlab="", ylab="Fluconazole concentration (mg/L)", main="", ylim=c(0,ylim_max), las=1, xaxt="n")
    axis(1, at=day_ticks, labels=date_labels, cex.axis=0.8); mtext("Date", side=1, line=2.5, cex=0.85)
    grid(col="#ECEFF1", lty=1)
    if(has_ci) polygon(c(tgrid_d,rev(tgrid_d)), c(p025,rev(p975)), col=grDevices::adjustcolor(COL_BAND,.35), border=NA)
    lines(tgrid_d, map_interp, lwd=3, col=COL_MAP)
    # Observations - always draw last, on top of everything
    if(length(obs_t)>0) {
      points(obs_t, obs_c, pch=21, cex=2.2, bg="#FF6F00", col="white", lwd=2)
    }
    # Dose annotations
    for(i in seq_len(nrow(hist_doses))){
      x_d <- hist_doses$time[i]/24
      text(x_d, ylim_max*0.97, paste0(hist_doses$amt[i],"mg"), cex=0.55, col="#00539F", srt=90, adj=c(1,0.5))
      segments(x_d, 0, x_d, ylim_max*0.06, col="#00539F", lwd=0.8, lty=3)
    }
    par(xpd=TRUE)
    legend("right", inset=c(-0.18,0), bty="o", box.col="#E0E0E0", bg="white", cex=0.78,
      legend=c("MAP curve",if(has_ci)"95% CI" else NULL,"Observations"),
      lty=c(1,if(has_ci)1 else NULL,if(has_ci)NA else NULL,NA), lwd=c(3,if(has_ci)2 else NULL,if(has_ci)NA else NULL,NA),
      pch=c(NA,if(has_ci)NA else NULL,if(has_ci)15 else NULL,21), pt.cex=c(NA,if(has_ci)NA else NULL,if(has_ci)2 else NULL,1.8),
      pt.bg=c(NA,if(has_ci)NA else NULL,if(has_ci)grDevices::adjustcolor(COL_BAND,.5) else NULL,"#FF6F00"),
      col=c(COL_MAP,if(has_ci)grDevices::adjustcolor(COL_BAND,.6) else NULL,"#FF6F00"))
    par(xpd=FALSE)
    }, error=function(e){message("Profile plot error:",e$message);plot.new();text(0.5,0.5,paste("Erreur:",e$message),cex=0.8)})
  })

  plot_density_panel<-function(r,param,xlab,title){
    tryCatch({
    eta_p<-rmvnorm_chol(5000,c(0,0),OMEGA_MAT);pv<-if(param=="CL")r$CLtyp*exp(eta_p[,1])else r$Vtyp*exp(eta_p[,2]);dv<-if(param=="CL")r$posterior$CL else r$posterior$V;mv<-if(param=="CL")r$map_psi["CL"]else r$map_psi["V"];cv<-if(param=="CL")r$cond_mean_psi["CL"]else r$cond_mean_psi["V"]
    d1<-density(pv,na.rm=TRUE);d2<-density(dv,na.rm=TRUE);ylm<-max(d1$y,d2$y)*1.18
    par(mar=c(4.5,4.5,3.5,1.5),bg="white");plot(d1,main=title,xlab=xlab,lwd=2,col="grey70",ylim=c(0,ylm),font.main=2,las=1);grid(col="#ECEFF1",lty=1)
    polygon(c(d2$x,rev(d2$x)),c(d2$y,rep(0,length(d2$y))),col=grDevices::adjustcolor(COL_CM,.22),border=NA);lines(d2,lwd=2.5,col=COL_CM);abline(v=mv,col=COL_MAP,lty=1,lwd=3);abline(v=cv,col=COL_OPT_M,lty=2,lwd=2)
    legend("topright",bty="n",legend=c("A priori","A posteriori","MAP","Moy. cond."),col=c("grey70",COL_CM,COL_MAP,COL_OPT_M),lty=c(1,1,1,2),lwd=c(2,2.5,3,2))
    }, error=function(e){message("Density plot error:",e$message);plot.new();text(0.5,0.5,paste("Erreur:",e$message),cex=0.8)})
  }
  output$ap_plot_cl<-renderPlot({tryCatch({if(is.null(rv$ap_results))return(NULL);plot_density_panel(rv$ap_results,"CL","CL (L/h)","CL distribution (prior vs posterior)")},error=function(e){plot.new();text(0.5,0.5,"Error",cex=0.8)})})
  output$ap_plot_v<-renderPlot({tryCatch({if(is.null(rv$ap_results))return(NULL);plot_density_panel(rv$ap_results,"V","V (L)","Distribution of Vd values (prior vs posterior)")},error=function(e){plot.new();text(0.5,0.5,"Error",cex=0.8)})})

  # ===================== EXPERT PK SIMULATION =====================

  # Dynamic per-day dose inputs (J2, J3, ... since J1 = Dose initiale)
  output$expert_day_doses <- renderUI({
    n <- input$expert_n_days; if(is.null(n)||!is.finite(n)||n<2) return(NULL)
    n_extra <- n - 1
    def_dose <- if(!is.null(input$expert_new_dose)) input$expert_new_dose else 400
    # Use fixed width per input, wrap in scrollable div if many
    inputs <- lapply(seq_len(n_extra), function(d) {
      day_num <- d + 1
      cur <- isolate(input[[paste0("expert_d",day_num)]])
      val <- if(!is.null(cur)&&is.finite(cur)) cur else def_dose
      div(style="display:inline-block;width:90px;margin-right:4px;", numericInput(paste0("expert_d",day_num), paste0("Day ",day_num), value=val, min=0, step=25))
    })
    div(style="overflow-x:auto;white-space:nowrap;padding:4px 0;", inputs)
  })

  # Fill all days with dose initiale
  observeEvent(input$fill_uniform, {
    n <- input$expert_n_days; d <- input$expert_new_dose
    if(is.null(n)||is.null(d)) return()
    for(i in 2:n) updateNumericInput(session, paste0("expert_d",i), value=d)
  }, ignoreInit=TRUE)

  # Read full schedule: J1 from dose initiale, J2+ from per-day inputs
  get_expert_schedule <- function() {
    n <- input$expert_n_days; if(is.null(n)||n<1) return(NULL)
    j1 <- if(!is.null(input$expert_new_dose)) max(0, input$expert_new_dose) else 400
    if(n==1) return(j1)
    rest <- sapply(2:n, function(d) {
      v <- input[[paste0("expert_d",d)]]
      if(is.null(v)||!is.finite(v)) j1 else max(0, v)
    })
    c(j1, rest)
  }

  observeEvent(input$run_expert_sim, {
    if(is.null(rv$ap_results)) { showNotification("Please run the Bayesian estimation first", type="warning"); return() }
    rv$expert_rec <- NULL
    rv$expert_sim_trigger <- runif(1)
  }, ignoreInit=TRUE)

  # Find optimal dosing SCHEME: adaptation + maintenance, including pause (dose=0)
  observeEvent(input$find_optimal_dose, {
    if(is.null(rv$ap_results)) { showNotification("Please run the Bayesian estimation first", type="warning"); return() }
    r <- rv$ap_results; apr <- rv$ap_inputs
    CL_ind <- r$map_psi["CL"]; V_ind <- r$map_psi["V"]; tinf_r <- r$tinf
    nmtran <- r$nmtran; if(is.null(nmtran)) return()
    hist_doses <- nmtran[nmtran$evid==1,]
    t_last_hist <- max(hist_doses$time)
    new_tau <- if(!is.null(input$expert_new_tau)) input$expert_new_tau else 24
    n_days <- if(!is.null(input$expert_n_days)) input$expert_n_days else 7
    target <- tryCatch(get_target_auc(), error=function(e) 400)

    withProgress(message="Searching for optimal dosage regimen...", value=0, {
      candidates <- c(0, seq(25, 2000, by=25))  # 0 = pause
      t_new_start <- t_last_hist + new_tau

      # Step 1: SS maintenance dose
      dose_ss <- ceiling(target * CL_ind / (24/new_tau) / 25) * 25
      dose_ss <- min(dose_ss, 2000)
      setProgress(0.2, detail="Maintenance SS dose...")

      # Step 2: Optimize each day
      schedule <- rep(dose_ss, n_days)
      for(d in seq_len(n_days)) {
        setProgress(0.2 + 0.7*d/n_days, detail=paste0("Optimising Day ",d,"..."))
        found <- FALSE
        for(d_cand in candidates) {
          test_sched <- schedule; test_sched[d] <- d_cand
          n_new <- n_days * ceiling(24/new_tau)
          new_times <- t_new_start + (seq_len(n_new)-1) * new_tau
          new_amts <- sapply(new_times, function(t) {
            day_idx <- max(1, min(floor((t - t_new_start)/24)+1, n_days))
            test_sched[day_idx]
          })
          all_dt <- c(hist_doses$time, new_times)
          all_da <- c(hist_doses$amt, new_amts)
          ts <- t_new_start + (d-1)*24; te <- t_new_start + d*24
          auc_d <- auc_model_interval(ts, te, all_dt, all_da, CL_ind, V_ind, tinf_r)
          if(auc_d >= target && auc_d <= target*1.5) { schedule[d] <- d_cand; found <- TRUE; break }
        }
        if(!found) {
          # If no dose keeps AUC in range, find closest to target
          best_cand <- candidates[1]; best_diff <- Inf
          for(d_cand in candidates) {
            test_sched <- schedule; test_sched[d] <- d_cand
            n_new <- n_days * ceiling(24/new_tau)
            new_times <- t_new_start + (seq_len(n_new)-1) * new_tau
            new_amts <- sapply(new_times, function(t) {
              day_idx <- max(1, min(floor((t - t_new_start)/24)+1, n_days))
              test_sched[day_idx]
            })
            all_dt <- c(hist_doses$time, new_times); all_da <- c(hist_doses$amt, new_amts)
            ts <- t_new_start + (d-1)*24; te <- t_new_start + d*24
            auc_d <- auc_model_interval(ts, te, all_dt, all_da, CL_ind, V_ind, tinf_r)
            diff_t <- abs(auc_d - target)
            if(auc_d >= target && diff_t < best_diff) { best_diff <- diff_t; best_cand <- d_cand }
          }
          schedule[d] <- best_cand
        }
        if(d < n_days && schedule[d] == dose_ss) { for(dd in (d+1):n_days) schedule[dd] <- dose_ss; break }
      }
      setProgress(1, detail="Done")
    })

    # Fill inputs: J1 = dose initiale, J2+ = per-day
    updateNumericInput(session, "expert_new_dose", value=schedule[1])
    for(i in 2:length(schedule)) updateNumericInput(session, paste0("expert_d",i), value=schedule[i])

    first_ss_day <- which(schedule == dose_ss)[1]
    if(is.na(first_ss_day)) first_ss_day <- n_days + 1
    adapt_days <- if(first_ss_day > 1) schedule[1:(first_ss_day-1)] else integer(0)
    has_pause <- any(schedule == 0)

    rv$expert_rec <- list(schedule=schedule, dose_ss=dose_ss, adapt_days=adapt_days,
                          first_ss_day=first_ss_day, tau=new_tau, target=target, n_days=n_days, has_pause=has_pause)
    rv$expert_sim_trigger <- runif(1)
  }, ignoreInit=TRUE)

  output$expert_recommendation <- renderUI({
    if(is.null(rv$expert_rec)) return(NULL)
    rec <- rv$expert_rec
    if(is.null(rec$schedule)) return(NULL)
    sched <- rec$schedule; dose_ss <- rec$dose_ss; target <- rec$target; tau <- rec$tau

    # Build description
    sched_txt <- paste(sapply(seq_along(sched), function(i) {
      if(sched[i]==0) paste0("<b>J",i," : PAUSE</b>")
      else paste0("<b>J",i," : ",sched[i]," mg</b>")
    }), collapse=" \u2192 ")

    has_pause <- any(sched==0)
    all_same <- length(unique(sched))==1

    if(all_same) {
      div(class="rec-box rec-above", style="margin:10px 0;", HTML(paste0(
        "\u2713 <b>Optimal regimen: ", sched[1], " mg / ", tau, "h</b>",
        " \u2014 AUC target \u2265 ", target, " mg\u00b7h/L achieved from Day 1.")))
    } else {
      note <- if(has_pause) "<br><em>\u26a0 Treatment pause proposed to avoid overdose.</em>" else ""
      div(class=paste0("rec-box ",if(has_pause)"rec-warn"else"rec-above"), style="margin:10px 0;", HTML(paste0(
        "\u2713 <b>Optimal dosage regimen:</b><br>",
        sched_txt,
        "<br>Maintenance dose: <b>", dose_ss, " mg / ", tau, "h</b>",
        note)))
    }
  })

  output$expert_pk_plot <- renderPlot({
    trg <- rv$expert_sim_trigger  # reactive dependency on both buttons
    if(is.null(rv$ap_results)) return(NULL)
    if(is.null(trg) && (is.null(input$run_expert_sim)||input$run_expert_sim==0)) return(NULL)
    r <- rv$ap_results; apr <- rv$ap_inputs
    rec <- rv$expert_rec
    new_tau <- if(!is.null(input$expert_new_tau)) input$expert_new_tau else 24
    n_days <- input$expert_n_days
    if(is.null(new_tau)||is.null(n_days)||n_days<1) return(NULL)

    # Read schedule: from rec if just optimized, else from per-day inputs
    schedule <- if(!is.null(rec) && !is.null(rec$schedule)) rec$schedule else get_expert_schedule()
    if(is.null(schedule)||length(schedule)==0) return(NULL)
    if(length(schedule)<n_days) schedule <- c(schedule, rep(schedule[length(schedule)], n_days-length(schedule)))

    CL_map <- r$map_psi["CL"]; V_map <- r$map_psi["V"]; tinf_r <- r$tinf
    post <- r$posterior
    nmtran <- r$nmtran; hist_doses <- nmtran[nmtran$evid==1,]; obs_rows <- nmtran[nmtran$evid==0,]
    t_last_hist <- max(hist_doses$time)

    t_new_start <- t_last_hist + new_tau
    n_new <- n_days * ceiling(24/new_tau)
    new_times <- t_new_start + (seq_len(n_new)-1) * new_tau
    new_amts <- sapply(new_times, function(t) {
      day_idx <- max(1, min(floor((t-t_new_start)/24)+1, length(schedule)))
      schedule[day_idx]
    })

    # All dose schedule
    all_times <- c(hist_doses$time, new_times)
    all_amts <- c(hist_doses$amt, new_amts)
    t_end <- max(new_times) + new_tau
    tgrid <- seq(0, t_end, length.out=600)

    # MAP curve by superposition
    k_map <- CL_map/V_map
    map_conc <- sapply(tgrid, function(tt) {
      c_tot <- 0
      for(i in seq_along(all_times)) { dt <- tt - all_times[i]; if(dt>0) c_tot <- c_tot + (all_amts[i]/tinf_r/CL_map)*ifelse(dt<=tinf_r,1-exp(-k_map*dt),(1-exp(-k_map*tinf_r))*exp(-k_map*(dt-tinf_r))) }
      max(c_tot,1e-12)
    })

    # Posterior CI (subsample)
    n_sub <- min(300, nrow(post)); idx_sub <- sample(nrow(post), n_sub)
    post_mat <- matrix(NA, n_sub, length(tgrid))
    for(s in seq_len(n_sub)) {
      cl_s <- post$CL[idx_sub[s]]; v_s <- post$V[idx_sub[s]]; k_s <- cl_s/v_s
      for(j in seq_along(tgrid)) {
        c_tot <- 0
        for(i in seq_along(all_times)) { dt <- tgrid[j]-all_times[i]; if(dt>0) c_tot <- c_tot + (all_amts[i]/tinf_r/cl_s)*ifelse(dt<=tinf_r,1-exp(-k_s*dt),(1-exp(-k_s*tinf_r))*exp(-k_s*(dt-tinf_r))) }
        post_mat[s,j] <- max(c_tot,1e-12)
      }
    }
    p50 <- apply(post_mat,2,median); p025 <- apply(post_mat,2,quantile,0.025); p975 <- apply(post_mat,2,quantile,0.975)

    tgrid_d <- tgrid/24
    ylim_max <- max(c(p975, obs_rows$DV, map_conc), na.rm=TRUE)*1.15

    # Date labels
    d1_date <- tryCatch(as.Date(apr$date_first, "%d/%m/%Y"), error=function(e) Sys.Date()-7)
    day_ticks <- seq(0, max(tgrid_d), by=1)
    date_labels <- format(d1_date + day_ticks, "%d/%m")

    par(mar=c(5,4.5,3,10), bg="white", xpd=FALSE)
    plot(tgrid_d, map_conc, type="n", xlab="", ylab="Fluconazole concentration (mg/L)", main="", ylim=c(0,ylim_max), las=1, xaxt="n")
    axis(1, at=day_ticks, labels=date_labels, cex.axis=0.8)
    mtext("Date", side=1, line=2.5, cex=0.85)
    grid(col="#ECEFF1", lty=1)
    # New regimen zone
    rect(t_new_start/24, 0, t_end/24, ylim_max, col=grDevices::adjustcolor("#E8F5E9",0.25), border=NA)
    # IC 95% band
    polygon(c(tgrid_d,rev(tgrid_d)), c(p025,rev(p975)), col=grDevices::adjustcolor(COL_BAND,.35), border=NA)
    lines(tgrid_d, map_conc, lwd=3, col=COL_MAP)
    # Transition line
    abline(v=t_new_start/24, col="#2E7D32", lty=2, lwd=2)
    # Observations
    points(obs_rows$time/24, obs_rows$DV, pch=21, cex=1.8, bg="#FF6F00", col="white", lwd=1.5)
    # Dose annotations
    for(i in seq_along(all_times)) {
      x_d <- all_times[i]/24
      col_d <- if(all_times[i] < t_new_start) "#00539F" else "#2E7D32"
      text(x_d, ylim_max*0.97, paste0(all_amts[i],"mg"), cex=0.5, col=col_d, srt=90, adj=c(1,0.5))
      segments(x_d, 0, x_d, ylim_max*0.06, col=col_d, lwd=0.8, lty=3)
    }
    par(xpd=TRUE)
    legend("right", inset=c(-0.18,0), bty="o", box.col="#E0E0E0", bg="white", cex=0.78,
      legend=c("MAP curve","95% CI","Observations","New regimen"),
      lty=c(1,NA,NA,2), lwd=c(3,NA,NA,2), pch=c(NA,15,21,NA), pt.cex=c(NA,2,1.4,NA),
      pt.bg=c(NA,grDevices::adjustcolor(COL_BAND,.5),"#FF6F00",NA),
      col=c(COL_MAP,grDevices::adjustcolor(COL_BAND,.6),"#FF6F00","#2E7D32"))
    par(xpd=FALSE)
  })

  output$expert_auc_table <- renderTable({
    trg <- rv$expert_sim_trigger
    if(is.null(rv$ap_results)) return(NULL)
    if(is.null(trg) && (is.null(input$run_expert_sim)||input$run_expert_sim==0)) return(NULL)
    r <- rv$ap_results; apr <- rv$ap_inputs; rec <- rv$expert_rec
    new_tau <- if(!is.null(input$expert_new_tau)) input$expert_new_tau else 24
    n_days <- input$expert_n_days
    if(is.null(new_tau)||is.null(n_days)||n_days<1) return(NULL)
    CL_ind <- r$map_psi["CL"]; V_ind <- r$map_psi["V"]; tinf_r <- r$tinf
    nmtran <- r$nmtran; if(is.null(nmtran)) return(NULL)
    hist_doses <- nmtran[nmtran$evid==1,]
    t_last_hist <- max(hist_doses$time)
    target <- tryCatch(get_target_auc(), error=function(e) 400)
    d1_date <- tryCatch(as.Date(apr$date_first, "%d/%m/%Y"), error=function(e) Sys.Date()-7)

    # Read schedule from rec or per-day inputs
    schedule <- if(!is.null(rec) && !is.null(rec$schedule)) rec$schedule else get_expert_schedule()
    if(is.null(schedule)||length(schedule)==0) return(NULL)
    if(length(schedule)<n_days) schedule <- c(schedule, rep(schedule[length(schedule)], n_days-length(schedule)))

    t_new_start <- t_last_hist + new_tau
    n_new <- n_days * ceiling(24/new_tau)
    new_times <- t_new_start + (seq_len(n_new)-1) * new_tau
    new_amts <- sapply(new_times, function(t) {
      day_idx <- max(1, min(floor((t-t_new_start)/24)+1, length(schedule)))
      schedule[day_idx]
    })
    all_dt <- c(hist_doses$time, new_times); all_da <- c(hist_doses$amt, new_amts)

    rows <- list()
    for(d in seq_len(n_days)) {
      ts <- t_new_start + (d-1)*24; te <- t_new_start + d*24
      auc_map <- auc_model_interval(ts, te, all_dt, all_da, CL_ind, V_ind, tinf_r)
      date_from <- format(d1_date + ts/24, "%d/%m/%Y"); date_to <- format(d1_date + te/24, "%d/%m/%Y")
      day_dose <- if(d <= length(schedule)) schedule[d] else schedule[length(schedule)]
      statut <- if(auc_map >= target) "\u2713 Achieved" else "\u2717 Not achieved"
      rows[[d]] <- data.frame(
        "Day" = paste0("Day ", d), "Dates" = paste0(date_from, " \u2013 ", date_to),
        "Dose (mg)" = if(day_dose==0) "PAUSE" else day_dose,
        "AUC 0-24h (mg.h/L)" = round(auc_map, 1), "Target" = target, "Status" = statut,
        check.names=FALSE, stringsAsFactors=FALSE)
    }
    auc_df <- do.call(rbind, rows)
    # SS row: use last day's dose as maintenance
    dose_last <- schedule[length(schedule)]
    if(dose_last > 0) {
      auc_ss <- auc_ss_24(dose_last, CL_ind, V_ind, tinf_r, new_tau)
      auc_df <- rbind(auc_df, data.frame(
        "Day" = "SS", "Dates" = "Steady-state",
        "Dose (mg)" = dose_last,
        "AUC 0-24h (mg.h/L)" = round(auc_ss, 1), "Target" = target,
        "Status" = if(auc_ss >= target) "\u2713 Achieved" else "\u2717 Not achieved",
        check.names=FALSE, stringsAsFactors=FALSE))
    }
    auc_df
  }, align="c")
  output$report_ui_aposteriori<-renderUI({if(is.null(rv$ap_results))return(NULL);downloadButton("dl_aposteriori","Download a posteriori report (HTML)",class="btn-primary",style="width:100%;font-size:16px;padding:12px;")})

  output$dl_aposteriori<-downloadHandler(
    filename=function(){paste0("FlucoShiny_Aposteriori_",gsub("[^0-9A-Za-z_-]","",patient$ipp),"_",format(Sys.Date(),"%Y%m%d"),".html")},
    content=function(file){
      tryCatch({
        if(is.null(rv$ap_results)){file.create(file);return()}
        logo_path<-Filter(file.exists,c("flucoshiny.png","Flucoshiny.png","FlucoShiny.png","www/flucoshiny.png","www/Flucoshiny.png","www/FlucoShiny.png"));logo_abs<-if(length(logo_path))normalizePath(logo_path[1])else NULL
        logo_html_val<-if(!is.null(logo_abs)&&file.exists(logo_abs)){b64<-knitr::image_uri(logo_abs);paste0('<img src="',b64,'" style="position:absolute;top:16px;right:22px;height:140px;max-width:400px;object-fit:contain;opacity:.92;" alt="logo">')}else ""
        tmpdir<-tempfile("fluco_post_");dir.create(tmpdir,recursive=TRUE)
        rmd_path<-file.path(tmpdir,"report_aposteriori.Rmd");write_rmd_aposteriori(rmd_path,logo_abs)
        rmarkdown::render(rmd_path,output_file=file,params=list(
          patient=list(last_name=patient$last_name,first_name=patient$first_name,dob=format(patient$dob,"%d/%m/%Y"),service=patient$service,ipp=patient$ipp,sexe=patient$sexe,uf=patient$uf,hopital=patient$hopital,adresse=patient$adresse),
          res=rv$ap_results,apr=rv$ap_inputs,expert_rec=rv$expert_rec,
          logo_path=logo_abs,logo_html=logo_html_val,nsim=NSIM,mic_label=get_mic_label(),target_auc=get_target_auc(),clinical=get_clinical(),pharmacologue=rv$pharmacologue,pk_model_ref=PK_MODEL_REF
        ),envir=new.env(parent=globalenv()),quiet=FALSE)
      }, error=function(e){
        message("Report aposteriori error: ",e$message)
        writeLines(paste0("<html><body><h2>Report generation error</h2><pre>",e$message,"</pre></body></html>"),file)
      })
    })
}

shinyApp(ui, server)
