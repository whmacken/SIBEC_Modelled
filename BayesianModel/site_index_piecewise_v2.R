# Standalone interface for site_index_piecewise_v2.stan. Sourcing never fits a model.
.si_pwv2_file <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
.si_pwv2_dir <- if(is.null(.si_pwv2_file)) getwd() else dirname(normalizePath(.si_pwv2_file))
.si_pwv2_restore <- function(x,template) {
  if(inherits(template,"data.table")) data.table::as.data.table(x) else x
}
.si_pwv2_int <- function(x) is.numeric(x) && length(x)==1L && is.finite(x) && x>=1 && x==floor(x)

# Preserve all current climate priors; append interpretable edatopic priors.
si_priors_piecewise_v2 <- function(climate_priors=read.csv(file.path(.si_pwv2_dir,"moisture_optimum_priors.csv"))) {
  p <- as.data.frame(climate_priors)
  for(v in intersect(c("A","C_median","C_log_sd","w_C_median","w_C_log_sd"),names(p))) p[[v]] <- NULL
  # Convert old prior centres once, never multiply sampled P by sampled cold fraction.
  if(!"P_median" %in% names(p)) p$P_median <- 30
  if(!"F_median" %in% names(p))
    p$F_median <- p$P_median * (if("cold200_median" %in% names(p)) p$cold200_median else .02)
  if(!"F_logit_sd" %in% names(p))
    p$F_logit_sd <- if("cold200_logit_sd" %in% names(p)) p$cold200_logit_sd else .5
  if(!"P_share_logit_sd" %in% names(p)) p$P_share_logit_sd <- .55
  if(!"K_median" %in% names(p)) p$K_median <- 1800
  if(!"K_gap_log_sd" %in% names(p)) p$K_gap_log_sd <- .4
  for(v in intersect(c("cold200_median","cold200_logit_sd","P_log_sd"),names(p))) p[[v]] <- NULL
  defaults <- c(O_M_median=4,O_M_logit_sd=.65,r0_median=.02,r0_logit_sd=.5,
    r8_median=.02,r8_logit_sd=.5,rA_median=.4,rA_log_sd=.4,b_NN_mean=-.5,b_NN_sd=.3)
  for(v in setdiff(names(defaults),names(p))) p[[v]] <- unname(defaults[v])
  .si_pwv2_restore(p,climate_priors)
}
.si_pwv2_priors <- function(priors) {
  p <- as.data.frame(priors)
  med <- paste0(c("K","V50","w_V","D50","p_CMD"),"_median")
  sds <- c("K_gap_log_sd",paste0(c("V50","w_V","D50","p_CMD"),"_log_sd"))
  cols <- c("P_median","P_share_logit_sd","F_median","F_logit_sd",med,sds,"beta_CMD_mean","beta_CMD_sd")
  if(!all(c("species",cols) %in% names(p)) || !nrow(p)) stop("Invalid prior columns; use si_priors_piecewise_v2().")
  p$species <- as.character(p$species)
  if(anyNA(p$species) || any(!nzchar(p$species)) || anyDuplicated(p$species)) stop("Prior species must be unique and nonempty.")
  vcols <- c("V50_median","V50_log_sd","w_V_median","w_V_log_sd")
  if(anyNA(p[vcols])) stop("Fill all four VPD prior columns (V50_median, V50_log_sd, w_V_median, w_V_log_sd). Medians use kPa.")
  if(!all(vapply(p[cols],is.numeric,logical(1))) || any(!is.finite(as.matrix(p[cols])))) stop("Priors must be finite numeric values.")
  if(any(as.matrix(p[setdiff(cols,"beta_CMD_mean")])<=0)) stop("Prior medians and SDs must be positive; beta_CMD_mean may have either sign.")
  if(any(p$F_median>=p$P_median | p$P_median>=50)) stop("Require 0 < F_median < P_median < 50 metres.")
  ec <- c("O_M_median","O_M_logit_sd","r0_median","r0_logit_sd",
    "r8_median","r8_logit_sd","rA_median","rA_log_sd","b_NN_mean","b_NN_sd")
  if(!all(ec %in% names(p)) || !all(vapply(p[ec],is.numeric,logical(1))) ||
     any(!is.finite(as.matrix(p[ec])))) stop("Supply finite edatopic priors; use si_priors_piecewise_v2().")
  if(any(p$O_M_median<=0 | p$O_M_median>=8 | p$r0_median<=0 | p$r0_median>=1 |
         p$r8_median<=0 | p$r8_median>=1 | p$rA_median<=0)) stop("Require 0 < O_M < 8, endpoint ratios in (0,1), and rA > 0.")
  es <- as.matrix(p[c("O_M_logit_sd","r0_logit_sd","r8_logit_sd","rA_log_sd","b_NN_sd")])
  if(any(es<=0)) stop("Edatopic SDs must be positive.")
  em <- cbind(qlogis(p$O_M_median/8),qlogis(p$r0_median),qlogis(p$r8_median),log(p$rA_median),p$b_NN_mean)
  if(any(p$K_median<=200)) stop("K_median must exceed the 200 DD5 floor.")
  medvals <- as.matrix(p[med]); medvals[,1] <- medvals[,1]-200
  list(table=p,ed_centre=em,ed_sd=es,centre=cbind(qlogis(p$F_median/50),log(medvals),p$beta_CMD_mean),
    sd=cbind(p$F_logit_sd,as.matrix(p[sds]),p$beta_CMD_sd))
}
.si_pwv2_data <- function(data,response=FALSE) {
  required <- c("species","DD5","VPD","CMD","nutrient","moisture",if(response) "site_index")
  if(!is.data.frame(data) || !nrow(data) || !all(required %in% names(data)))
    stop("Required columns: ",paste(required,collapse=", "),". VPD must be in kPa.")
  d <- as.data.frame(data)
  d$species <- as.character(d$species)
  if(anyNA(d$species) || any(!nzchar(d$species))) stop("Invalid species.")
  labels <- as.character(d$nutrient)
  n <- suppressWarnings(as.numeric(labels)); k <- match(labels,LETTERS[1:5])
  n[!is.na(k)] <- k[!is.na(k)]
  if(any(!is.finite(n)) || any(n<1 | n>5)) stop("nutrient must be A-E or numeric SNR in [1,5].")
  d$snr_numeric <- n
  for(v in c("DD5","VPD","CMD","moisture",if(response) "site_index"))
    if(!is.numeric(d[[v]]) || any(!is.finite(d[[v]]))) stop("Invalid numeric ",v)
  if(any(d$DD5<0 | d$VPD<0 | d$CMD<0)) stop("Climate variables must be nonnegative.")
  if(any(d$moisture<0 | d$moisture>8)) stop("SMR must be in [0,8].")
  if(response && any(d$site_index<=0)) stop("SI must be positive.")
  d
}
prepare_si_data_piecewise_v2 <- function(data,priors,n_check=500L,
    save_log_lik=FALSE,prior_only=FALSE,seed=41,
    gamma_CMD_mean=0,gamma_CMD_mean_sd=.25,gamma_CMD_tau_sd=.20) {
  d <- .si_pwv2_data(data,TRUE)
  p <- as.data.frame(priors)
  if(!"species" %in% names(p)) stop("Missing species in prior table.")
  p <- p[p$species %in% unique(d$species),,drop=FALSE]
  pp <- .si_pwv2_priors(p); p <- pp$table
  sp <- match(d$species,p$species)
  if(anyNA(sp)) stop("Species missing from prior table: ",paste(unique(d$species[is.na(sp)]),collapse=", "))
  if(!.si_pwv2_int(n_check)) stop("n_check must be a positive integer.")
  for(x in list(save_log_lik,prior_only))
    if(!is.logical(x) || length(x)!=1L || is.na(x)) stop("Flags must be TRUE or FALSE.")
  set.seed(seed); N <- nrow(d); idx <- sort(sample.int(N,min(N,n_check)))
  if(any(!is.finite(c(gamma_CMD_mean,gamma_CMD_mean_sd,gamma_CMD_tau_sd))) ||
     length(gamma_CMD_mean)!=1 || length(gamma_CMD_mean_sd)!=1 ||
     length(gamma_CMD_tau_sd)!=1 || gamma_CMD_mean_sd<=0 || gamma_CMD_tau_sd<=0)
    stop("Supply one finite gamma mean and positive hyperprior SDs.")
  sd <- list(gamma_CMD_mean=gamma_CMD_mean,gamma_CMD_mean_sd=gamma_CMD_mean_sd,
    gamma_CMD_tau_sd=gamma_CMD_tau_sd,N=N,S=nrow(p),y=d$site_index,DD5=d$DD5,VPD=d$VPD,CMD=d$CMD,
    snr=d$snr_numeric,smr=d$moisture,species=sp,P_share_centres=qlogis((p$P_median-p$F_median)/(50-p$F_median)),P_share_sd=p$P_share_logit_sd,
    cold_anchor_DD5=200L, prior_centre=unname(pp$centre),prior_sd=unname(pp$sd),
    ed_prior_centre=unname(pp$ed_centre),ed_prior_sd=unname(pp$ed_sd),
    prior_only=as.integer(prior_only),N_check=length(idx),check_idx=idx,
    N_log_lik=if(save_log_lik) N else 0L)
  list(model="dd5_piecewise_v2_smr_cmd_temperature",stan_data=sd,data=.si_pwv2_restore(d,data),
    priors=.si_pwv2_restore(p,priors),species=p$species,check_idx=idx,
    climate_units=c(DD5="degree-days",VPD="kPa",CMD="mm"))
}
fit_si_piecewise_v2 <- function(prepared,
    stan_file=file.path(.si_pwv2_dir,"site_index_piecewise_v2.stan"),
    output_dir="si_fits_piecewise_v2",seed=2026,chains=4,parallel_chains=chains,
    iter_warmup=1000,iter_sampling=1000,adapt_delta=.95,max_treedepth=12,step_size=.001) {
  if(!identical(prepared$model,"dd5_piecewise_v2_smr_cmd_temperature")) stop("Use prepare_si_data_piecewise_v2().")
  if(!identical(prepared$stan_data$cold_anchor_DD5,200L)) stop("Rebuild prepared data for the DD5=200 anchor.")
  if(is.null(prepared$stan_data$P_share_centres)) stop("Rebuild prepared data for the 50 m bounded-plateau model.")
  if(!requireNamespace("cmdstanr",quietly=TRUE)) stop("cmdstanr is required.")
  S <- prepared$stan_data$S
  init <- list(gamma_CMD_bar=prepared$stan_data$gamma_CMD_mean,
    tau_gamma_CMD=.10,z_gamma_CMD=rep(0,S),delta_P=0,rho_P=.5,z_P=rep(0,S),z_climate=matrix(0,S,7),
    delta_climate=rep(0,6),rho_climate=rep(.5,6),delta_ed=rep(0,5),rho_ed=rep(.5,5),
    z_ed=matrix(0,S,5),c_M_bar=0,tau_c_M=.05,z_c_M=rep(0,S),
    log_sigma_bar=log(.25),tau_log_sigma=.15,z_log_sigma=rep(0,S))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  run_dir <- tempfile("run_",tmpdir=normalizePath(output_dir)); dir.create(run_dir)
  saveRDS(prepared,file.path(run_dir,"prepared.rds"))
  model <- cmdstanr::cmdstan_model(stan_file)
  fit <- model$sample(data=prepared$stan_data,init=replicate(chains,init,simplify=FALSE),
    seed=seed,chains=chains,parallel_chains=parallel_chains,
    iter_warmup=iter_warmup,iter_sampling=iter_sampling,adapt_delta=adapt_delta,
    max_treedepth=max_treedepth,step_size=step_size,output_dir=run_dir,sig_figs=10)
  fit$save_object(file.path(run_dir,"fit.rds"))
  list(fit=fit,prepared=prepared,run_dir=run_dir)
}

# Log multipliers other than DD5 and amplitude, shared by predictions and
# the observation adjustment. getv returns one vector of posterior draws.
.si_pwv2_log_site <- function(getv,s,n,smr,VPD,CMD,DD5) {
  m <- (smr-3)/3; opt <- getv("O_M",s)
  ed <- ifelse(smr<opt,getv("log_r0",s)*((smr-opt)/opt)^2,
    getv("log_r8",s)*((smr-opt)/(8-opt))^2) +
    (getv("b_NN",s)-getv("log_rA",s))*n+getv("b_NN",s)*n^2
  loss <- if(CMD==0) 0 else log(2)*exp(getv("p_CMD",s)*
    (log(CMD)-log(getv("D50",s))-getv("beta_CMD",s)*m+getv("gamma_CMD",s)*(DD5-1500)/1000))
  ed-getv("nutrient_qmax",s)+plogis((getv("V50",s)-VPD)/getv("w_V",s),log.p=TRUE)-loss
}

.si_pwv2_prediction_vars <- c("log_P","nutrient_qmax","cold_floor","K_gap","V50","w_V","D50","p_CMD","beta_CMD","gamma_CMD","c_M","O_M","log_r0","log_r8","log_rA","b_NN","sigma_species")
# One joint posterior scenario per row, one input observation per column.
draw_si_piecewise_v2 <- function(object,newdata,n_draws=100L,
    type=c("mean","predictive"),draw_ids=NULL,seed=NULL) {
  type <- match.arg(type)
  if(!identical(object$prepared$model,"dd5_piecewise_v2_smr_cmd_temperature")) stop("Requires a piecewise_v2 DD5 fit.")
  d <- .si_pwv2_data(newdata); sp <- match(d$species,object$prepared$species)
  if(anyNA(sp)) stop("Unknown species.")
  if(!.si_pwv2_int(n_draws)) stop("n_draws must be positive integer.")
  if(!is.null(seed)) set.seed(seed)
  vars <- .si_pwv2_prediction_vars
  draws <- if(!is.null(object$draws)) as.matrix(object$draws) else
    as.matrix(object$fit$draws(variables=vars,format="matrix"))
  if(!nrow(draws)) stop("No draws available.")
  if(is.null(draw_ids)) draw_ids <- sample.int(nrow(draws),n_draws,replace=n_draws>nrow(draws))
  else {
    if(!is.numeric(draw_ids) || !length(draw_ids) || any(!is.finite(draw_ids)) ||
       any(draw_ids!=floor(draw_ids) | draw_ids<1 | draw_ids>nrow(draws))) stop("Invalid draw_ids.")
    if(!missing(n_draws) && n_draws!=length(draw_ids)) stop("n_draws must equal length(draw_ids).")
  }
  draws <- draws[draw_ids,,drop=FALSE]
  getv <- function(v,s) draws[,paste0(v,"[",s,"]")]
  eta <- matrix(NA_real_,nrow(draws),nrow(d))
  for(i in seq_len(nrow(d))) {
    s <- sp[i]; n <- (d$snr_numeric[i]-3)/2; m <- (d$moisture[i]-3)/3
    ramp <- pmin(1,pmax(0,d$DD5[i]-200)/(getv("K_gap",s)*exp(getv("c_M",s)*m)))
    floor_frac <- getv("cold_floor",s)
    eta[,i] <- getv("log_P",s)+.si_pwv2_log_site(getv,s,n,d$moisture[i],d$VPD[i],d$CMD[i],d$DD5[i])+
      log(floor_frac+(1-floor_frac)*ramp)
  }
  values <- exp(eta)
  if(type=="predictive") {
    sig <- draws[,paste0("sigma_species[",sp,"]"),drop=FALSE]
    values <- exp(eta-.5*sig^2+matrix(rnorm(length(eta)),nrow(eta))*sig)
  }
  if(any(!is.finite(values))) stop("Nonfinite predictions; inspect posterior draws.")
  attr(values,"posterior_draw_ids") <- draw_ids
  attr(values,"log_mean") <- eta
  values
}

# Posterior summaries use the same parameter draw across all conditioning rows.
predict_si_piecewise_v2 <- function(object,newdata,interval=.90,n_draws=1000L,
    type=c("mean","predictive"),seed=73,return_draws=FALSE) {
  type <- match.arg(type)
  if(!is.numeric(interval) || length(interval)!=1L || !is.finite(interval) || interval<=0 || interval>=1)
    stop("interval must be between 0 and 1.")
  y <- draw_si_piecewise_v2(object,newdata,n_draws=n_draws,type=type,seed=seed)
  q <- t(vapply(seq_len(ncol(y)),function(i) quantile(y[,i],
    c((1-interval)/2,.5,(1+interval)/2),names=FALSE),numeric(3)))
  d <- as.data.frame(newdata)
  d$lower <- q[,1]; d$median <- q[,2]; d$upper <- q[,3]
  result <- .si_pwv2_restore(d,newdata)
  if(return_draws) return(list(summary=result,draws=y,
    posterior_draw_ids=attr(y,"posterior_draw_ids"),interval=interval,type=type))
  result
}
# One output row per input row x species. Uses the same posterior scenarios
# throughout; batch_size limits how many prediction columns are held in memory.
predict_si_validation_piecewise_v2 <- function(object,newdata,species=NULL,
    interval=.90,n_draws=1000L,type=c("mean","predictive"),seed=73,batch_size=256L) {
  type <- match.arg(type)
  if(!identical(object$prepared$model,"dd5_piecewise_v2_smr_cmd_temperature")) stop("Requires a piecewise_v2 DD5 fit.")
  fitted <- object$prepared$species
  if(is.null(species)) species <- fitted
  if(!is.character(species) || !length(species) || anyNA(species) ||
     any(!species %in% fitted) || anyDuplicated(species)) stop("Select distinct fitted species codes.")
  if(!is.numeric(interval) || length(interval)!=1L || !is.finite(interval) || interval<=0 || interval>=1)
    stop("interval must be between 0 and 1.")
  if(!.si_pwv2_int(n_draws) || !.si_pwv2_int(batch_size)) stop("Draw and batch sizes must be positive integers.")
  cols <- c("moisture","nutrient","DD5","CMD","VPD")
  if(!is.data.frame(newdata) || !nrow(newdata) || !all(cols %in% names(newdata)))
    stop("Supply a nonempty predictor table with: ",paste(cols,collapse=", "),".")
  d <- as.data.frame(newdata)
  reserved <- c("validation_row","SI_median","SI_lower","SI_upper")
  if(anyDuplicated(names(d)) || any(reserved %in% names(d))) stop("Input columns must be unique; reserved names: ",paste(reserved,collapse=", "),".")
  if("species" %in% names(d)) {
    if("observed_species" %in% names(d)) stop("Rename either species or observed_species to preserve both.")
    names(d)[match("species",names(d))] <- "observed_species"
  }
  probe <- d[cols]; probe$species <- species[1]; .si_pwv2_data(probe)
  m <- if(!is.null(object$draws)) as.matrix(object$draws) else
    as.matrix(object$fit$draws(variables=.si_pwv2_prediction_vars,format="matrix"))
  if(!nrow(m)) stop("No posterior draws available.")
  if(!is.null(seed)) set.seed(seed)
  ids <- sample.int(nrow(m),n_draws,replace=n_draws>nrow(m))
  local <- list(prepared=object$prepared,draws=m[ids,,drop=FALSE])
  idx <- rep(seq_len(nrow(d)),each=length(species))
  out <- d[idx,,drop=FALSE]; rownames(out) <- NULL
  out$validation_row <- idx
  out$species <- rep(species,nrow(d))
  q <- matrix(NA_real_,nrow(out),3)
  probs <- c((1-interval)/2,.5,(1+interval)/2)
  for(start in seq.int(1L,nrow(out),by=batch_size)) {
    ii <- seq.int(start,min(nrow(out),start+batch_size-1))
    y <- draw_si_piecewise_v2(local,out[ii,c("species",cols),drop=FALSE],
      draw_ids=seq_along(ids),type=type)
    q[ii,] <- t(vapply(seq_along(ii),function(j) quantile(y[,j],probs,names=FALSE),numeric(3)))
  }
  out$SI_median <- q[,2]; out$SI_lower <- q[,1]; out$SI_upper <- q[,3]
  .si_pwv2_restore(out,newdata)
}

.si_pwv2_limits <- function(x,name) {
  if(!is.null(x) && (!is.numeric(x) || length(x)!=2L || any(!is.finite(x)) ||
    x[1]<0 || x[2]<=x[1])) stop(name," must contain two increasing nonnegative values.")
}
.si_pwv2_climate <- function(DD5,VPD,CMD) {
  for(x in list(DD5,VPD,CMD)) if(!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<0)
    stop("Supply one finite nonnegative value for each fixed climate variable.")
}

# Align every observed point to the chosen reference site/climate. No residual
# noise is generated. DD5 stays fixed; the full conditional SI ratio adjusts y.
adjust_si_data_piecewise_v2 <- function(object,data=NULL,species=NULL,
    VPD=1,CMD=200,nutrient="C",moisture=3,interval=.90,n_draws=1000L,
    seed=73,draw_ids=NULL,return_draws=FALSE) {
  if(!identical(object$prepared$model,"dd5_piecewise_v2_smr_cmd_temperature")) stop("Requires a piecewise_v2 DD5 fit.")
  .si_pwv2_climate(0,VPD,CMD)
  if(length(nutrient)!=1L || length(moisture)!=1L) stop("Choose one SNR and SMR.")
  if(!is.numeric(interval) || length(interval)!=1L || !is.finite(interval) || interval<=0 || interval>=1) stop("interval must be between 0 and 1.")
  if(!.si_pwv2_int(n_draws)) stop("n_draws must be a positive integer.")
  template <- if(is.null(data)) object$prepared$data else data
  d <- .si_pwv2_data(template,TRUE)
  if(is.null(species)) species <- unique(d$species)
  if(!is.character(species) || !length(species) || anyNA(species) || any(!species %in% object$prepared$species)) stop("Unknown species.")
  rows <- which(d$species %in% species)
  if(!length(rows)) stop("No observed rows for the selected species.")
  d <- d[rows,,drop=FALSE]
  ref <- .si_pwv2_data(data.frame(species=species[1],DD5=0,VPD=VPD,CMD=CMD,nutrient=nutrient,moisture=moisture))
  vars <- c("K_gap","cold_floor","nutrient_qmax","c_M","V50","w_V","D50","p_CMD","beta_CMD","gamma_CMD","O_M","log_r0","log_r8","log_rA","b_NN")
  m <- if(!is.null(object$draws)) as.matrix(object$draws) else as.matrix(object$fit$draws(variables=vars,format="matrix"))
  if(!nrow(m)) stop("No posterior draws.")
  if(is.null(draw_ids)) {
    if(!is.null(seed)) set.seed(seed)
    draw_ids <- sample.int(nrow(m),n_draws,replace=n_draws>nrow(m))
  } else if(!is.numeric(draw_ids) || !length(draw_ids) || any(!is.finite(draw_ids)) ||
    any(draw_ids!=floor(draw_ids) | draw_ids<1 | draw_ids>nrow(m))) stop("Invalid draw_ids.")
  m <- m[draw_ids,,drop=FALSE]
  getv <- function(v,s) m[,paste0(v,"[",s,"]")]
  ref_m <- (moisture-3)/3; ref_n <- (ref$snr_numeric-3)/2
  probs <- c((1-interval)/2,.5,(1+interval)/2)
  xs <- ys <- matrix(NA_real_,nrow(d),3)
  xd <- yd <- if(return_draws) matrix(NA_real_,nrow(m),nrow(d)) else NULL
  for(sp in unique(d$species)) {
    s <- match(sp,object$prepared$species)
    cm <- getv("c_M",s)
    for(i in which(d$species==sp)) {
      orig_m <- (d$moisture[i]-3)/3
      x <- rep(d$DD5[i],length(cm))
      ref_site <- .si_pwv2_log_site(getv,s,ref_n,moisture,VPD,CMD,d$DD5[i])
      ramp_ref <- pmin(1,pmax(0,d$DD5[i]-200)/(getv("K_gap",s)*exp(cm*ref_m)))
      ramp_orig <- pmin(1,pmax(0,d$DD5[i]-200)/(getv("K_gap",s)*exp(cm*orig_m)))
      floor <- getv("cold_floor",s)
      ref_site <- ref_site+log(floor+(1-floor)*ramp_ref)-log(floor+(1-floor)*ramp_orig)
      orig_site <- .si_pwv2_log_site(getv,s,(d$snr_numeric[i]-3)/2,d$moisture[i],d$VPD[i],d$CMD[i],d$DD5[i])
      y <- d$site_index[i]*exp(ref_site-orig_site)
      if(any(!is.finite(x)) || any(!is.finite(y))) stop("Nonfinite adjustment at input row ",rows[i],"; inspect that observation and the reference conditions.")
      xs[i,] <- quantile(x,probs,names=FALSE); ys[i,] <- quantile(y,probs,names=FALSE)
      if(return_draws) { xd[,i] <- x; yd[,i] <- y }
    }
  }
  d$observation_row <- rows
  d$DD5_aligned_lower <- xs[,1]; d$DD5_aligned <- xs[,2]; d$DD5_aligned_upper <- xs[,3]
  d$SI_adjusted_lower <- ys[,1]; d$SI_adjusted <- ys[,2]; d$SI_adjusted_upper <- ys[,3]
  d$reference_VPD <- VPD; d$reference_CMD <- CMD; d$reference_nutrient <- nutrient; d$reference_moisture <- moisture
  result <- .si_pwv2_restore(d,template)
  if(return_draws) list(summary=result,DD5_aligned=xd,SI_adjusted=yd,posterior_draw_ids=draw_ids) else result
}

# Conditional expected-SI curves, not marginal curves or prediction intervals.
plot_si_curves_piecewise_v2 <- function(object,variable=c("DD5","VPD","CMD"),
    species=NULL,DD5=2000,VPD=1,CMD=200,nutrient="C",moisture=3,
    grid=NULL,xlim=NULL,ylim=NULL,interval=.90,n_draws=1000L,
    individual_draws=0L,seed=73,plot=TRUE,observations=FALSE,data=NULL,
    adjustment_intervals=FALSE,point_cex=.5) {
  variable <- match.arg(variable)
  if(observations && variable!="DD5") stop("Aligned observation overlays are defined for DD5 only.")
  .si_pwv2_climate(DD5,VPD,CMD); .si_pwv2_limits(xlim,"xlim"); .si_pwv2_limits(ylim,"ylim")
  if(length(nutrient)!=1L || length(moisture)!=1L) stop("Choose one SNR and SMR.")
  if(is.null(species)) species <- object$prepared$species
  if(!is.character(species) || !length(species) || anyNA(species) || any(!species %in% object$prepared$species)) stop("Unknown species.")
  species <- unique(species)
  if(!is.numeric(individual_draws) || length(individual_draws)!=1L || !is.finite(individual_draws) ||
     individual_draws<0 || individual_draws!=floor(individual_draws)) stop("individual_draws must be a nonnegative integer.")
  if(is.null(grid)) {
    limits <- if(!is.null(xlim)) xlim else switch(variable,DD5=c(0,6000),VPD=c(0,4),CMD=c(0,2000))
    grid <- seq(limits[1],limits[2],length.out=201)
  }
  if(!is.numeric(grid) || length(grid)<2L || any(!is.finite(grid)) || any(grid<0) || length(unique(grid))<2L)
    stop("grid must contain at least two distinct nonnegative finite values.")
  grid <- sort(unique(grid))
  if(!is.null(xlim)) grid <- sort(unique(c(grid,seq(xlim[1],xlim[2],length.out=201))))
  if(is.null(xlim)) xlim <- range(grid)
  d <- expand.grid(x=grid,species=species,stringsAsFactors=FALSE)
  d$DD5 <- DD5; d$VPD <- VPD; d$CMD <- CMD
  d[[variable]] <- d$x; d$nutrient <- nutrient; d$moisture <- moisture
  pr <- predict_si_piecewise_v2(object,d,interval=interval,n_draws=n_draws,seed=seed,return_draws=TRUE)
  p <- as.data.frame(pr$summary)
  obs_y <- NULL
  if(observations) {
    pr$observations <- adjust_si_data_piecewise_v2(object,data=data,species=species,
      VPD=VPD,CMD=CMD,nutrient=nutrient,moisture=moisture,interval=interval,draw_ids=pr$posterior_draw_ids)
    adj <- as.data.frame(pr$observations)
    visible <- adj$DD5_aligned>=xlim[1] & adj$DD5_aligned<=xlim[2]
    obs_y <- if(adjustment_intervals) adj$SI_adjusted_upper[visible] else adj$SI_adjusted[visible]
    if(any(!visible)) message(sum(!visible)," aligned observations outside xlim; retained in $observations.")
  }
  if(is.null(ylim)) ylim <- c(0,max(1,p$upper[p$x>=xlim[1] & p$x<=xlim[2]],obs_y))
  if(plot) {
    old <- par(no.readonly=TRUE); on.exit(par(old))
    nc <- min(2L,length(species)); par(mfrow=c(ceiling(length(species)/nc),nc),mar=c(4,4,4,1))
    for(sp in species) {
      keep <- which(p$species==sp); z <- p[keep,,drop=FALSE]
      plot(z$x,z$median,type="n",xlim=xlim,ylim=ylim,main=sp,
        xlab=if(observations) "DD5" else switch(variable,DD5="DD5",VPD="Summer VPD (kPa)",CMD="CMD (mm)"),
        ylab=if(observations) "Site index (m)" else "Expected SI (m)")
      polygon(c(z$x,rev(z$x)),c(z$lower,rev(z$upper)),col="#0072B230",border=NA)
      if(individual_draws>0) matlines(z$x,t(pr$draws[seq_len(min(individual_draws,nrow(pr$draws))),keep,drop=FALSE]),
        col="#0072B230",lty=1)
      if(observations) {
        w <- adj[adj$species==sp,,drop=FALSE]
        if(adjustment_intervals) {
          segments(w$DD5_aligned_lower,w$SI_adjusted,w$DD5_aligned_upper,w$SI_adjusted,col="#44444425")
          segments(w$DD5_aligned,w$SI_adjusted_lower,w$DD5_aligned,w$SI_adjusted_upper,col="#44444425")
        }
        points(w$DD5_aligned,w$SI_adjusted,pch=16,cex=point_cex,col="#44444450")
      }
      lines(z$x,z$median,col="#0072B2",lwd=2)
      fixed <- c(paste0("DD5 ",DD5),paste0("VPD ",VPD," kPa"),paste0("CMD ",CMD," mm"))
      mtext(paste(c(fixed[c("DD5","VPD","CMD")!=variable],paste("SNR",nutrient),paste("SMR",moisture)),collapse=" | "),cex=.65)
      obs <- object$prepared$data
      # Marginal climate coverage only; rug does not condition on site variables.
      if(observations) rug(w$DD5_aligned,col="#00000030") else if(!is.null(obs)) rug(obs[[variable]][obs[["species"]]==sp],col="#00000030")
      legend("topright",c("Posterior median",paste0(100*interval,"% credible interval"),if(observations) "Adjusted observations"),
        col=c("#0072B2","#0072B230",if(observations) "#444444"),lty=c(1,NA,if(observations) NA),
        pch=c(NA,15,if(observations) 16),bty="n",cex=.65)
    }
  }
  invisible(pr)
}

# Square A-E by SMR 0-8 grid; numbers are median and central credible interval.
plot_si_edatopic_piecewise_v2 <- function(object,species,DD5,VPD,CMD,
    SMR=0:8,interval=.90,n_draws=1000L,digits=1,zlim=NULL,
    text_cex=.65,seed=73,plot=TRUE) {
  .si_pwv2_climate(DD5,VPD,CMD); .si_pwv2_limits(zlim,"zlim")
  if(!is.character(species) || length(species)!=1L || is.na(species) || !species %in% object$prepared$species)
    stop("Choose one fitted species.")
  if(!is.numeric(SMR) || !length(SMR) || any(!is.finite(SMR)) || any(SMR<0 | SMR>8) ||
     any(SMR!=floor(SMR))) stop("For grid cells, SMR must contain integers in [0,8].")
  SMR <- sort(unique(SMR))
  if(!is.numeric(digits) || length(digits)!=1L || !digits %in% 0:3) stop("digits must be 0-3.")
  if(!is.numeric(text_cex) || length(text_cex)!=1L || !is.finite(text_cex) || text_cex<=0) stop("Invalid text_cex.")
  d <- expand.grid(nutrient=LETTERS[1:5],moisture=SMR,stringsAsFactors=FALSE)
  d$species <- species; d$DD5 <- DD5; d$VPD <- VPD; d$CMD <- CMD
  p <- as.data.frame(predict_si_piecewise_v2(object,d,interval=interval,n_draws=n_draws,seed=seed))
  if(is.null(zlim)) zlim <- c(0,max(1,p$median))
  palette <- grDevices::hcl.colors(100,"YlGnBu",rev=TRUE)
  idx <- pmax(1L,pmin(100L,1L+floor(99*(p$median-zlim[1])/diff(zlim))))
  if(plot) {
    old <- par(no.readonly=TRUE); on.exit({layout(1);par(old)})
    layout(matrix(1:2,nrow=1),widths=c(5,1))
    par(mar=c(4,4,5,0))
    plot(NA,xlim=c(.5,5.5),ylim=c(max(SMR)+.5,min(SMR)-.5),asp=1,
      xaxs="i",yaxs="i",axes=FALSE,xlab="Soil nutrient regime",ylab="Soil moisture regime",
      main="")
    title(main=paste(species,"- expected SI (m)"),line=3.5)
    axis(1,at=1:5,labels=LETTERS[1:5]); axis(2,at=SMR,las=1)
    mtext(paste0("DD5 ",DD5," | VPD ",VPD," kPa | CMD ",CMD," mm"),line=2,cex=.7)
    mtext(paste0("Median [",100*interval,"% credible interval]"),line=.6,cex=.7)
    for(i in seq_len(nrow(p))) {
      x <- match(p$nutrient[i],LETTERS[1:5]); y <- p$moisture[i]
      fill <- palette[idx[i]]
      rect(x-.5,y-.5,x+.5,y+.5,col=fill,border="white")
      rgb <- grDevices::col2rgb(fill)/255
      ink <- if(sum(rgb*c(.2126,.7152,.0722))<.5) "white" else "black"
      label <- paste0(formatC(p$median[i],format="f",digits=digits),"\n[",
        formatC(p$lower[i],format="f",digits=digits),", ",formatC(p$upper[i],format="f",digits=digits),"]")
      text(x,y,label,cex=text_cex,col=ink)
    }
    par(mar=c(4,1,5,3))
    plot(NA,xlim=c(0,1),ylim=zlim,xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab="")
    breaks <- seq(zlim[1],zlim[2],length.out=101)
    rect(0,breaks[-101],1,breaks[-1],col=palette,border=NA)
    axis(4,las=1); mtext("SI (m)",side=3,line=.5,cex=.7)
  }
  invisible(.si_pwv2_restore(p,object$prepared$data))
}

# Climate-component prior checks (not full SI prior predictive simulations).
