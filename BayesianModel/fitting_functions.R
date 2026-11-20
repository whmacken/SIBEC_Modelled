# Standalone interface for site_index_piecewise.stan. Sourcing never fits a model.
.si_pw_file <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
.si_pw_dir <- if(is.null(.si_pw_file)) getwd() else dirname(normalizePath(.si_pw_file))
.si_pw_restore <- function(x,template) {
  if(inherits(template,"data.table")) data.table::as.data.table(x) else x
}
.si_pw_int <- function(x) is.numeric(x) && length(x)==1L && is.finite(x) && x>=1 && x==floor(x)

# Preserve all current climate priors; append interpretable edatopic priors.
si_priors_piecewise <- function(climate_priors=read.csv(file.path(.si_pw_dir,"moisture_optimum_priors.csv"))) {
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
  .si_pw_restore(p,climate_priors)
}
.si_pw_priors <- function(priors) {
  p <- as.data.frame(priors)
  med <- paste0(c("K","V50","w_V","D50","p_CMD"),"_median")
  sds <- c("K_gap_log_sd",paste0(c("V50","w_V","D50","p_CMD"),"_log_sd"))
  cols <- c("P_median","P_share_logit_sd","F_median","F_logit_sd",med,sds,"beta_CMD_mean","beta_CMD_sd")
  if(!all(c("species",cols) %in% names(p)) || !nrow(p)) stop("Invalid prior columns; use si_priors_piecewise().")
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
     any(!is.finite(as.matrix(p[ec])))) stop("Supply finite edatopic priors; use si_priors_piecewise().")
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
.si_pw_data <- function(data,response=FALSE) {
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
prepare_si_data_piecewise <- function(data,priors,n_check=500L,
    save_log_lik=FALSE,prior_only=FALSE,seed=41) {
  d <- .si_pw_data(data,TRUE)
  p <- as.data.frame(priors)
  if(!"species" %in% names(p)) stop("Missing species in prior table.")
  p <- p[p$species %in% unique(d$species),,drop=FALSE]
  pp <- .si_pw_priors(p); p <- pp$table
  sp <- match(d$species,p$species)
  if(anyNA(sp)) stop("Species missing from prior table: ",paste(unique(d$species[is.na(sp)]),collapse=", "))
  if(!.si_pw_int(n_check)) stop("n_check must be a positive integer.")
  for(x in list(save_log_lik,prior_only))
    if(!is.logical(x) || length(x)!=1L || is.na(x)) stop("Flags must be TRUE or FALSE.")
  set.seed(seed); N <- nrow(d); idx <- sort(sample.int(N,min(N,n_check)))
  sd <- list(N=N,S=nrow(p),y=d$site_index,DD5=d$DD5,VPD=d$VPD,CMD=d$CMD,
    snr=d$snr_numeric,smr=d$moisture,species=sp,P_share_centres=qlogis((p$P_median-p$F_median)/(50-p$F_median)),P_share_sd=p$P_share_logit_sd,
    cold_anchor_DD5=200L, prior_centre=unname(pp$centre),prior_sd=unname(pp$sd),
    ed_prior_centre=unname(pp$ed_centre),ed_prior_sd=unname(pp$ed_sd),
    prior_only=as.integer(prior_only),N_check=length(idx),check_idx=idx,
    N_log_lik=if(save_log_lik) N else 0L)
  list(model="dd5_piecewise_smr",stan_data=sd,data=.si_pw_restore(d,data),
    priors=.si_pw_restore(p,priors),species=p$species,check_idx=idx,
    climate_units=c(DD5="degree-days",VPD="kPa",CMD="mm"))
}
fit_si_piecewise <- function(prepared,
    stan_file=file.path(.si_pw_dir,"site_index_piecewise.stan"),
    output_dir="si_fits_piecewise",seed=2026,chains=4,parallel_chains=chains,
    iter_warmup=1000,iter_sampling=1000,adapt_delta=.95,max_treedepth=12,step_size=.001) {
  if(!identical(prepared$model,"dd5_piecewise_smr")) stop("Use prepare_si_data_piecewise().")
  if(!identical(prepared$stan_data$cold_anchor_DD5,200L)) stop("Rebuild prepared data for the DD5=200 anchor.")
  if(is.null(prepared$stan_data$P_share_centres)) stop("Rebuild prepared data for the 50 m bounded-plateau model.")
  if(!requireNamespace("cmdstanr",quietly=TRUE)) stop("cmdstanr is required.")
  S <- prepared$stan_data$S
  init <- list(delta_P=0,rho_P=.5,z_P=rep(0,S),z_climate=matrix(0,S,7),
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
.si_pw_log_site <- function(getv,s,n,smr,VPD,CMD) {
  m <- (smr-3)/3; opt <- getv("O_M",s)
  ed <- ifelse(smr<opt,getv("log_r0",s)*((smr-opt)/opt)^2,
    getv("log_r8",s)*((smr-opt)/(8-opt))^2) +
    (getv("b_NN",s)-getv("log_rA",s))*n+getv("b_NN",s)*n^2
  loss <- if(CMD==0) 0 else log(2)*exp(getv("p_CMD",s)*
    (log(CMD)-log(getv("D50",s))-getv("beta_CMD",s)*m))
  ed-getv("nutrient_qmax",s)+plogis((getv("V50",s)-VPD)/getv("w_V",s),log.p=TRUE)-loss
}

.si_pw_prediction_vars <- c("log_P","nutrient_qmax","cold_floor","K_gap","V50","w_V","D50","p_CMD","beta_CMD","c_M","O_M","log_r0","log_r8","log_rA","b_NN","sigma_species")
# One joint posterior scenario per row, one input observation per column.
draw_si_piecewise <- function(object,newdata,n_draws=100L,
    type=c("mean","predictive"),draw_ids=NULL,seed=NULL) {
  type <- match.arg(type)
  if(!identical(object$prepared$model,"dd5_piecewise_smr")) stop("Requires a piecewise DD5 fit.")
  d <- .si_pw_data(newdata); sp <- match(d$species,object$prepared$species)
  if(anyNA(sp)) stop("Unknown species.")
  if(!.si_pw_int(n_draws)) stop("n_draws must be positive integer.")
  if(!is.null(seed)) set.seed(seed)
  vars <- .si_pw_prediction_vars
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
    eta[,i] <- getv("log_P",s)+.si_pw_log_site(getv,s,n,d$moisture[i],d$VPD[i],d$CMD[i])+
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
predict_si_piecewise <- function(object,newdata,interval=.90,n_draws=1000L,
    type=c("mean","predictive"),seed=73,return_draws=FALSE) {
  type <- match.arg(type)
  if(!is.numeric(interval) || length(interval)!=1L || !is.finite(interval) || interval<=0 || interval>=1)
    stop("interval must be between 0 and 1.")
  y <- draw_si_piecewise(object,newdata,n_draws=n_draws,type=type,seed=seed)
  q <- t(vapply(seq_len(ncol(y)),function(i) quantile(y[,i],
    c((1-interval)/2,.5,(1+interval)/2),names=FALSE),numeric(3)))
  d <- as.data.frame(newdata)
  d$lower <- q[,1]; d$median <- q[,2]; d$upper <- q[,3]
  result <- .si_pw_restore(d,newdata)
  if(return_draws) return(list(summary=result,draws=y,
    posterior_draw_ids=attr(y,"posterior_draw_ids"),interval=interval,type=type))
  result
}
# One output row per input row x species. Uses the same posterior scenarios
# throughout; batch_size limits how many prediction columns are held in memory.
predict_si_validation_piecewise <- function(object,newdata,species=NULL,
    interval=.90,n_draws=1000L,type=c("mean","predictive"),seed=73,batch_size=256L) {
  type <- match.arg(type)
  if(!identical(object$prepared$model,"dd5_piecewise_smr")) stop("Requires a piecewise DD5 fit.")
  fitted <- object$prepared$species
  if(is.null(species)) species <- fitted
  if(!is.character(species) || !length(species) || anyNA(species) ||
     any(!species %in% fitted) || anyDuplicated(species)) stop("Select distinct fitted species codes.")
  if(!is.numeric(interval) || length(interval)!=1L || !is.finite(interval) || interval<=0 || interval>=1)
    stop("interval must be between 0 and 1.")
  if(!.si_pw_int(n_draws) || !.si_pw_int(batch_size)) stop("Draw and batch sizes must be positive integers.")
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
  probe <- d[cols]; probe$species <- species[1]; .si_pw_data(probe)
  m <- if(!is.null(object$draws)) as.matrix(object$draws) else
    as.matrix(object$fit$draws(variables=.si_pw_prediction_vars,format="matrix"))
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
    y <- draw_si_piecewise(local,out[ii,c("species",cols),drop=FALSE],
      draw_ids=seq_along(ids),type=type)
    q[ii,] <- t(vapply(seq_along(ii),function(j) quantile(y[,j],probs,names=FALSE),numeric(3)))
  }
  out$SI_median <- q[,2]; out$SI_lower <- q[,1]; out$SI_upper <- q[,3]
  .si_pw_restore(out,newdata)
}

.si_pw_limits <- function(x,name) {
  if(!is.null(x) && (!is.numeric(x) || length(x)!=2L || any(!is.finite(x)) ||
    x[1]<0 || x[2]<=x[1])) stop(name," must contain two increasing nonnegative values.")
}
.si_pw_climate <- function(DD5,VPD,CMD) {
  for(x in list(DD5,VPD,CMD)) if(!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<0)
    stop("Supply one finite nonnegative value for each fixed climate variable.")
}

# Align every observed point to the chosen reference site/climate. No residual
# noise is generated. DD5 <= 200 stays fixed; higher DD5 is rescaled around the 200 anchor.
adjust_si_data_piecewise <- function(object,data=NULL,species=NULL,
    VPD=1,CMD=200,nutrient="C",moisture=3,interval=.90,n_draws=1000L,
    seed=73,draw_ids=NULL,return_draws=FALSE) {
  if(!identical(object$prepared$model,"dd5_piecewise_smr")) stop("Requires a piecewise DD5 fit.")
  .si_pw_climate(0,VPD,CMD)
  if(length(nutrient)!=1L || length(moisture)!=1L) stop("Choose one SNR and SMR.")
  if(!is.numeric(interval) || length(interval)!=1L || !is.finite(interval) || interval<=0 || interval>=1) stop("interval must be between 0 and 1.")
  if(!.si_pw_int(n_draws)) stop("n_draws must be a positive integer.")
  template <- if(is.null(data)) object$prepared$data else data
  d <- .si_pw_data(template,TRUE)
  if(is.null(species)) species <- unique(d$species)
  if(!is.character(species) || !length(species) || anyNA(species) || any(!species %in% object$prepared$species)) stop("Unknown species.")
  rows <- which(d$species %in% species)
  if(!length(rows)) stop("No observed rows for the selected species.")
  d <- d[rows,,drop=FALSE]
  ref <- .si_pw_data(data.frame(species=species[1],DD5=0,VPD=VPD,CMD=CMD,nutrient=nutrient,moisture=moisture))
  vars <- c("nutrient_qmax","c_M","V50","w_V","D50","p_CMD","beta_CMD","O_M","log_r0","log_r8","log_rA","b_NN")
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
    ref_site <- .si_pw_log_site(getv,s,ref_n,moisture,VPD,CMD)
    cm <- getv("c_M",s)
    for(i in which(d$species==sp)) {
      orig_m <- (d$moisture[i]-3)/3
      x <- if(d$DD5[i]<=200) rep(d$DD5[i],length(cm)) else 200+(d$DD5[i]-200)*exp(cm*(ref_m-orig_m))
      orig_site <- .si_pw_log_site(getv,s,(d$snr_numeric[i]-3)/2,d$moisture[i],d$VPD[i],d$CMD[i])
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
  result <- .si_pw_restore(d,template)
  if(return_draws) list(summary=result,DD5_aligned=xd,SI_adjusted=yd,posterior_draw_ids=draw_ids) else result
}

# Conditional expected-SI curves, not marginal curves or prediction intervals.
plot_si_curves_piecewise <- function(object,variable=c("DD5","VPD","CMD"),
    species=NULL,DD5=2000,VPD=1,CMD=200,nutrient="C",moisture=3,
    grid=NULL,xlim=NULL,ylim=NULL,interval=.90,n_draws=1000L,
    individual_draws=0L,seed=73,plot=TRUE,observations=FALSE,data=NULL,
    adjustment_intervals=FALSE,point_cex=.5) {
  variable <- match.arg(variable)
  if(observations && variable!="DD5") stop("Aligned observation overlays are defined for DD5 only.")
  .si_pw_climate(DD5,VPD,CMD); .si_pw_limits(xlim,"xlim"); .si_pw_limits(ylim,"ylim")
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
  pr <- predict_si_piecewise(object,d,interval=interval,n_draws=n_draws,seed=seed,return_draws=TRUE)
  p <- as.data.frame(pr$summary)
  obs_y <- NULL
  if(observations) {
    pr$observations <- adjust_si_data_piecewise(object,data=data,species=species,
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
        xlab=if(observations) "DD5 aligned to reference SMR" else switch(variable,DD5="DD5",VPD="Summer VPD (kPa)",CMD="CMD (mm)"),
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
plot_si_edatopic_piecewise <- function(object,species,DD5,VPD,CMD,
    SMR=0:8,interval=.90,n_draws=1000L,digits=1,zlim=NULL,
    text_cex=.65,seed=73,plot=TRUE) {
  .si_pw_climate(DD5,VPD,CMD); .si_pw_limits(zlim,"zlim")
  if(!is.character(species) || length(species)!=1L || is.na(species) || !species %in% object$prepared$species)
    stop("Choose one fitted species.")
  if(!is.numeric(SMR) || !length(SMR) || any(!is.finite(SMR)) || any(SMR<0 | SMR>8) ||
     any(SMR!=floor(SMR))) stop("For grid cells, SMR must contain integers in [0,8].")
  SMR <- sort(unique(SMR))
  if(!is.numeric(digits) || length(digits)!=1L || !digits %in% 0:3) stop("digits must be 0-3.")
  if(!is.numeric(text_cex) || length(text_cex)!=1L || !is.finite(text_cex) || text_cex<=0) stop("Invalid text_cex.")
  d <- expand.grid(nutrient=LETTERS[1:5],moisture=SMR,stringsAsFactors=FALSE)
  d$species <- species; d$DD5 <- DD5; d$VPD <- VPD; d$CMD <- CMD
  p <- as.data.frame(predict_si_piecewise(object,d,interval=interval,n_draws=n_draws,seed=seed))
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
  invisible(.si_pw_restore(p,object$prepared$data))
}

# Climate-component prior checks (not full SI prior predictive simulations).
check_si_climate_piecewise_priors <- function(priors,species=NULL,n_draws=1000L,n_curves=40L,
    DD5_grid=seq(0,6000,length.out=201),VPD_grid=seq(0,4,length.out=201),
    CMD_grid=seq(0,2000,length.out=201),SMR=c(2,3,5),seed=87,plot=TRUE) {
  p <- as.data.frame(priors)
  if(!is.null(species)) {
    if(!is.character(species) || !length(species) || anyNA(species) || any(!species %in% p$species)) stop("Unknown species.")
    p <- p[p$species %in% species,,drop=FALSE]
  }
  pp <- .si_pw_priors(p); p <- pp$table
  if(is.null(species)) species <- p$species
  if(!is.character(species) || !length(species) || anyNA(species) || any(!species %in% p$species)) stop("Unknown species.")
  if(!.si_pw_int(n_draws) || !.si_pw_int(n_curves)) stop("Draw counts must be positive integers.")
  for(g in list(DD5_grid,VPD_grid,CMD_grid))
    if(!is.numeric(g) || length(g)<2 || any(!is.finite(g)) || any(g<0) || is.unsorted(g,strictly=TRUE)) stop("Grids must be nonnegative, finite, strictly increasing.")
  if(!is.numeric(SMR) || !length(SMR) || any(!is.finite(SMR)) || any(SMR<0 | SMR>8)) stop("SMR must be in [0,8].")
  set.seed(seed); S <- nrow(p)
  delta <- matrix(rnorm(n_draws*6),n_draws,6)
  rho <- matrix(rbeta(n_draws*6,2,2),n_draws,6)
  cm_bar <- rnorm(n_draws,0,.15); cm_tau <- abs(rnorm(n_draws,0,.10))
  P_shared <- rnorm(n_draws); P_rho <- rbeta(n_draws,2,2)
  samples <- curves <- list()
  if(plot) { old <- par(no.readonly=TRUE); on.exit(par(old)); par(mfrow=c(length(unique(species)),3),mar=c(4,4,3,1)) }
  for(s in seq_len(S)) {
    z <- matrix(rnorm(n_draws*7),n_draws,7)
    z[,2:7] <- sqrt(rho)*delta+sqrt(1-rho)*z[,2:7]
    th <- sweep(sweep(z,2,pp$sd[s,],"*"),2,pp$centre[s,],"+")
    if(!p$species[s] %in% species) next
    cm <- cm_bar+cm_tau*rnorm(n_draws)
    F_m <- 50*plogis(th[,1])
    share0 <- (p$P_median[s]-p$F_median[s])/(50-p$F_median[s])
    P <- F_m+(50-F_m)*plogis(qlogis(share0)+p$P_share_logit_sd[s]*
      (sqrt(P_rho)*P_shared+sqrt(1-P_rho)*rnorm(n_draws)))
    a <- F_m/P
    vals <- cbind(P,200+exp(th[,2]),exp(th[,2:6,drop=FALSE]),th[,7],a,cm,F_m)
    colnames(vals) <- c("P","K_ref","K_gap","V50","w_V","D50","p_CMD","beta_CMD","cold_floor","c_M","F_m")
    samples[[p$species[s]]] <- vals
    ramp <- outer(1/vals[,"K_gap"],pmax(0,DD5_grid-200))
    cold <- a+(1-a)*pmin(ramp,1)
    vpd <- plogis((vals[,"V50"]-outer(rep(1,n_draws),VPD_grid))/vals[,"w_V"])
    cmd <- lapply(SMR,function(smr) {
      ratio <- outer(1/(vals[,"D50"]*exp(vals[,"beta_CMD"]*(smr-3)/3)),CMD_grid)
      exp(-log(2)*ratio^vals[,"p_CMD"])
    }); names(cmd) <- as.character(SMR)
    curves[[p$species[s]]] <- list(DD5=cold,DD5_upper_SI=P*cold,VPD=vpd,CMD=cmd)
    if(plot) {
      ids <- seq_len(min(n_curves,n_draws))
      panel <- function(x,y,xlab) {
        q <- apply(y,2,quantile,c(.05,.5,.95))
        plot(x,q[2,],type="n",ylim=c(0,1),xlab=xlab,ylab="Multiplier",main=p$species[s])
        polygon(c(x,rev(x)),c(q[1,],rev(q[3,])),col="#0072B230",border=NA)
        matlines(x,t(y[ids,,drop=FALSE]),col="#0072B230",lty=1)
        lines(x,q[2,],lwd=2); abline(h=.5,lty=3,col="grey")
      }
      panel(DD5_grid,cold,"DD5 (SMR 3)"); panel(VPD_grid,vpd,"Summer VPD (kPa)")
      plot(CMD_grid,rep(1,length(CMD_grid)),type="n",ylim=c(0,1),xlab="CMD (mm)",ylab="Multiplier",main=p$species[s])
      colors <- grDevices::hcl.colors(length(SMR),"Dark 3")
      for(k in seq_along(SMR)) {
        q <- apply(cmd[[k]],2,quantile,c(.05,.5,.95))
        polygon(c(CMD_grid,rev(CMD_grid)),c(q[1,],rev(q[3,])),col=grDevices::adjustcolor(colors[k],.12),border=NA)
        lines(CMD_grid,q[2,],col=colors[k],lwd=2)
      }
      legend("topright",paste("SMR",SMR),col=colors,lty=1,bty="n",cex=.7)
    }
  }
  invisible(list(parameter_draws=samples,curves=curves,DD5_grid=DD5_grid,VPD_grid=VPD_grid,CMD_grid=CMD_grid,SMR=SMR))
}

# Prior expected SI (m) at DD5=200 over the 45 edatopes. No residual draws:
# VPD=CMD=0 gives the least penalized fixed nonnegative drought climate.
check_si_piecewise_cold200_priors <- function(priors,species,VPD=0,CMD=0,
    n_draws=1000L,seed=91,plot=TRUE) {
  .si_pw_climate(200,VPD,CMD)
  if(!.si_pw_int(n_draws) || n_draws<2L) stop("Use at least two prior draws for the SI grid.")
  if(!is.character(species) || !length(species) || anyNA(species)) stop("Specify species codes.")
  species <- unique(species)
  p <- as.data.frame(priors)
  # Draw independent parameter blocks with separate streams; keep the shared
  # population components aligned across species in each simulated scenario.
  climate <- check_si_climate_piecewise_priors(p,species,n_draws=n_draws,
    DD5_grid=c(0,200),VPD_grid=c(0,1),CMD_grid=c(0,1),seed=seed,plot=FALSE)
  ed <- check_si_piecewise_priors(p,species,n_draws=n_draws,
    seed=(as.double(seed)+1)%%.Machine$integer.max,plot=FALSE)
  set.seed((as.double(seed)+2)%%.Machine$integer.max)

  out <- list(); draws <- list(); max_cold <- list()
  if(plot) { old <- par(no.readonly=TRUE); on.exit(par(old))
    nc <- min(2L,length(species)); par(mfrow=c(ceiling(length(species)/nc),nc),mar=c(4,4,4,1)) }
  for(sp in species) {
    v <- climate$parameter_draws[[sp]]; grid <- ed$surfaces[[sp]]$grid
    A <- v[,"P"]; m <- (grid$moisture-3)/3
    cold <- matrix(v[,"cold_floor"],n_draws,nrow(grid))
    d50 <- v[,"D50"]*exp(outer(v[,"beta_CMD"],m))
    water <- exp(-log(2)*(CMD/d50)^v[,"p_CMD"])
    si <- A*ed$surfaces[[sp]]$draws*cold*water*plogis((v[,"V50"]-VPD)/v[,"w_V"])
    q <- t(apply(si,2,quantile,c(.05,.5,.95),names=FALSE))
    grid$species <- sp; grid$DD5 <- 200; grid$VPD <- VPD; grid$CMD <- CMD
    grid$lower <- q[,1]; grid$median <- q[,2]; grid$upper <- q[,3]
    out[[sp]] <- grid; draws[[sp]] <- si; max_cold[[sp]] <- apply(cold,1,max)
    if(plot) {
      image(1:5,0:8,matrix(grid$median,5,9),axes=FALSE,asp=1,ylim=c(8.5,-.5),
        xlab="SNR",ylab="SMR",col=hcl.colors(100,"YlOrRd",rev=TRUE),
        main=paste(sp,"prior median SI (m)"))
      axis(1,at=1:5,labels=LETTERS[1:5]); axis(2,at=0:8,las=1,pos=.5)
      mtext(paste("DD5 200 | VPD",VPD,"kPa | CMD",CMD,"mm"),cex=.65)
      for(i in seq_len(nrow(grid))) text(match(grid$nutrient[i],LETTERS[1:5]),
        grid$moisture[i],sprintf("%.2f",grid$median[i]),cex=.65,
        col=if(grid$median[i] > mean(range(grid$median))) "white" else "black")
    }
  }
  invisible(list(summary=.si_pw_restore(do.call(rbind,out),priors),
    SI_draws=draws,max_cold_multiplier=max_cold))
}

# Prior edatopic response, independent of climate. Ratios apply at all SNRs.
check_si_piecewise_priors <- function(priors,species=NULL,n_draws=1000L,
    n_curves=40L,seed=91,plot=TRUE) {
  p <- as.data.frame(priors)
  if(!is.null(species)) {
    if(!is.character(species) || !length(species) || anyNA(species) || any(!species %in% p$species)) stop("Unknown species.")
    p <- p[p$species %in% species,,drop=FALSE]
  }
  pp <- .si_pw_priors(p); p <- pp$table
  if(is.null(species)) species <- p$species
  if(!is.character(species) || !length(species) || anyNA(species) || any(!species %in% p$species)) stop("Unknown species.")
  if(!.si_pw_int(n_draws) || !.si_pw_int(n_curves)) stop("Draw counts must be positive integers.")
  set.seed(seed)
  delta <- matrix(rnorm(n_draws*5),n_draws,5)
  rho <- matrix(rbeta(n_draws*5,2,2),n_draws,5)
  parameters <- surfaces <- list(); rows <- list()
  moisture <- seq(0,8,length.out=161)
  if(plot) { old <- par(no.readonly=TRUE); on.exit(par(old)); par(mfrow=c(length(unique(species)),2),mar=c(4,4,3,1)) }
  for(s in seq_len(nrow(p))) {
    z <- sqrt(rho)*delta+sqrt(1-rho)*matrix(rnorm(n_draws*5),n_draws,5)
    th <- sweep(sweep(z,2,pp$ed_sd[s,],"*"),2,pp$ed_centre[s,],"+")
    if(!p$species[s] %in% species) next
    opt <- 8*plogis(th[,1]); lr0 <- plogis(th[,2],log.p=TRUE); lr8 <- plogis(th[,3],log.p=TRUE)
    pars <- data.frame(O_M=opt,r0=exp(lr0),r8=exp(lr8),rA=exp(th[,4]),b_NN=th[,5])
    parameters[[p$species[s]]] <- pars
    for(v in names(pars)) {
      q <- quantile(pars[[v]],c(.05,.5,.95),names=FALSE)
      rows[[length(rows)+1L]] <- data.frame(species=p$species[s],parameter=v,lower=q[1],median=q[2],upper=q[3])
    }
    f <- vapply(moisture,function(x) exp(ifelse(x<opt,lr0*((x-opt)/opt)^2,
      lr8*((x-opt)/(8-opt))^2)),numeric(n_draws))
    grid <- expand.grid(nutrient=LETTERS[1:5],moisture=0:8,stringsAsFactors=FALSE)
    lin <- th[,5]-th[,4]; b <- th[,5]
    qmax <- pmax(b-lin,b+lin)
    inside <- b<0 & abs(lin)< -2*b
    qmax[inside] <- -lin[inside]^2/(4*b[inside])
    surface <- vapply(seq_len(nrow(grid)),function(i) {
      n <- (match(grid$nutrient[i],LETTERS[1:5])-3)/2; x <- grid$moisture[i]
      exp(ifelse(x<opt,lr0*((x-opt)/opt)^2,lr8*((x-opt)/(8-opt))^2)+
        (th[,5]-th[,4])*n+th[,5]*n^2-qmax)
    },numeric(n_draws))
    surfaces[[p$species[s]]] <- list(grid=grid,draws=surface)
    if(plot) {
      q <- apply(f,2,quantile,c(.05,.5,.95))
      plot(moisture,q[2,],type="n",ylim=c(0,1),main=paste(p$species[s],"moisture"),xlab="SMR",ylab="Fraction of species optimum")
      polygon(c(moisture,rev(moisture)),c(q[1,],rev(q[3,])),col="#0072B230",border=NA)
      matlines(moisture,t(f[seq_len(min(n_curves,n_draws)),,drop=FALSE]),col="#0072B230",lty=1)
      lines(moisture,q[2,],lwd=2); abline(h=.02,col="#D55E00",lty=3)
      med <- apply(surface,2,median)
      palette <- hcl.colors(100,"YlGnBu",rev=TRUE)
      image(1:5,0:8,matrix(med,5,9),axes=FALSE,ylim=c(8.5,-.5),asp=1,
        col=palette,xlab="SNR",ylab="SMR",
        main=paste(p$species[s],"median relative productivity"))
      axis(1,at=1:5,labels=LETTERS[1:5]); axis(2,at=0:8,las=1,pos=.5)
      colors <- palette[pmin(100,pmax(1,1+floor(99*(med-min(med))/max(diff(range(med)),1e-10))))]
      for(i in seq_len(nrow(grid))) {
        lum <- sum(col2rgb(colors[i])/255*c(.2126,.7152,.0722))
        text(match(grid$nutrient[i],LETTERS[1:5]),grid$moisture[i],sprintf("%.2f",med[i]),
          cex=.65,col=if(lum<.5) "white" else "black")
      }
    }
  }
  invisible(list(summary=.si_pw_restore(do.call(rbind,rows),priors),
    parameter_draws=parameters,surfaces=surfaces))
}



# Prior expected SI at the baseline-best edatope for EACH draw. With nonzero
# drought effects the best realized SI can occur at a different SMR; P still
# bounds expected SI everywhere. Summaries use all draws, thin lines a subset.
check_si_piecewise_plateau_priors <- function(priors,species=NULL,VPD=0,CMD=0,
    DD5_grid=seq(0,6000,by=20),n_draws=1000L,n_curves=40L,seed=91,plot=TRUE) {
  .si_pw_climate(0,VPD,CMD)
  if(!.si_pw_int(n_draws) || n_draws<2L || !.si_pw_int(n_curves)) stop("Need at least two draws and one plotted curve.")
  cl <- check_si_climate_piecewise_priors(priors,species=species,n_draws=n_draws,n_curves=n_curves,
    DD5_grid=DD5_grid,VPD_grid=c(0,1),CMD_grid=c(0,1),seed=seed,plot=FALSE)
  spp <- names(cl$parameter_draws)
  p <- as.data.frame(priors); p <- p[match(spp,p$species),,drop=FALSE]
  set.seed((as.double(seed)+1)%%.Machine$integer.max)
  shared <- rnorm(n_draws); rho <- rbeta(n_draws,2,2)
  curves <- rows <- list()
  if(plot) { old <- par(no.readonly=TRUE); on.exit(par(old)); nc <- min(2L,length(spp))
    par(mfrow=c(ceiling(length(spp)/nc),nc),mar=c(4,4,3,1)) }
  for(i in seq_along(spp)) {
    sp <- spp[i]; v <- cl$parameter_draws[[sp]]
    opt <- 8*plogis(qlogis(p$O_M_median[i]/8)+p$O_M_logit_sd[i]*
      (sqrt(rho)*shared+sqrt(1-rho)*rnorm(n_draws)))
    m <- (opt-3)/3
    gap <- v[,"K_gap"]*exp(v[,"c_M"]*m)
    u <- pmin(outer(1/gap,pmax(0,DD5_grid-200)),1)
    f <- v[,"cold_floor"]+(1-v[,"cold_floor"])*u
    wet <- exp(-log(2)*(CMD/(v[,"D50"]*exp(v[,"beta_CMD"]*m)))^v[,"p_CMD"])
    si <- v[,"P"]*wet*plogis((v[,"V50"]-VPD)/v[,"w_V"])*f
    curves[[sp]] <- si
    rows[[i]] <- data.frame(species=sp,P_median=median(v[,"P"]),
      P_99=unname(quantile(v[,"P"],.99)),Pr_P_above_50=mean(v[,"P"]>50),
      Pr_P_above_49=mean(v[,"P"]>49),F_median=median(v[,"F_m"]),
      K_ref_median=median(v[,"K_ref"]))
    if(plot) {
      q <- apply(si,2,quantile,c(.05,.5,.95)); ids <- seq_len(min(n_curves,n_draws))
      plot(DD5_grid,q[2,],type="n",ylim=c(0,max(52,q,si[ids,])),main=sp,
        xlab="DD5",ylab="Prior expected SI (m)")
      polygon(c(DD5_grid,rev(DD5_grid)),c(q[1,],rev(q[3,])),col="#0072B230",border=NA)
      matlines(DD5_grid,t(si[ids,,drop=FALSE]),col="#0072B240",lty=1)
      lines(DD5_grid,q[2,],lwd=2); abline(h=50,lty=2,col="#B54A00")
      mtext(paste("Best baseline edatope | VPD",VPD,"kPa | CMD",CMD,"mm"),cex=.65)
    }
  }
  invisible(list(summary=.si_pw_restore(do.call(rbind,rows),priors),SI_draws=curves,DD5_grid=DD5_grid))
}

# Joint posterior diagnostic: where does each thermal response level off?
check_si_piecewise_plateaus <- function(object,moisture=4) {
  if(is.null(object$prepared$stan_data$P_share_centres))
    stop("Requires a refit with the bounded plateau.")
  if(!is.numeric(moisture) || length(moisture)!=1L || !is.finite(moisture) ||
     moisture<0 || moisture>8) stop("Choose one SMR between 0 and 8.")
  m <- if(!is.null(object$draws)) as.matrix(object$draws) else
    as.matrix(object$fit$draws(c("P","F_m","K_gap","c_M"),format="matrix"))
  d <- object$prepared$data
  rows <- lapply(seq_along(object$prepared$species),function(i) {
    s <- object$prepared$species[i]
    v <- function(n) m[,paste0(n,"[",i,"]")]
    K <- 200+v("K_gap")*exp(v("c_M")*(moisture-3)/3)
    max_dd5 <- max(d[["DD5"]][d[["species"]]==s])
    q <- function(x,p) unname(quantile(x,p))
    data.frame(species=s,SMR=moisture,F_median=q(v("F_m"),.5),
      P_lower=q(v("P"),.05),P_median=q(v("P"),.5),P_upper=q(v("P"),.95),
      Pr_P_above_49=mean(v("P")>49),K_lower=q(K,.05),K_median=q(K,.5),
      K_upper=q(K,.95),max_observed_DD5=max_dd5,Pr_K_beyond_max_DD5=mean(K>max_dd5))
  })
  .si_pw_restore(do.call(rbind,rows),d)
}

plot_si_residuals <- function(
    object,
    data,
    variable = "PAS",
    species = NULL,
    response_col = "site_index",
    n_draws = 500L,
    seed = 73,
    log_x = FALSE,
    standardised = FALSE,
    smooth = TRUE,
    span = 0.9
) {
  dt <- data.table::copy(data.table::as.data.table(data))
  
  required <- unique(c(
    "species", "DD5", "VPD", "CMD", "moisture", "nutrient",
    variable, response_col
  ))
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols)) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dt[, observation_row := .I]
  
  if (!is.null(species)) {
    selected_species <- species
    dt <- dt[species %in% selected_species]
  }
  
  if (!is.numeric(dt[[variable]]) ||
      !is.numeric(dt[[response_col]])) {
    stop("The diagnostic variable and observed SI must be numeric.")
  }
  
  keep <- is.finite(dt[[variable]]) &
    is.finite(dt[[response_col]]) & dt[[response_col]] > 0
  
  if (any(!keep)) {
    warning(sum(!keep), " rows omitted because the variable or SI is invalid.")
  }
  dt <- dt[keep]
  if (!nrow(dt)) stop("No observations remain.")
  
  if (log_x && any(dt[[variable]] < 0)) {
    stop("log_x = TRUE requires non-negative values.")
  }
  
  predictions <- draw_si_piecewise(
    object,
    newdata = dt,
    n_draws = n_draws,
    type = "mean",
    seed = seed
  )
  
  draw_ids <- attr(predictions, "posterior_draw_ids")
  log_mu <- attr(predictions, "log_mean")
  
  posterior <- if (!is.null(object$draws)) {
    as.matrix(object$draws)
  } else {
    as.matrix(object$fit$draws(
      variables = "sigma_species",
      format = "matrix"
    ))
  }
  
  species_index <- match(dt$species, object$prepared$species)
  sigma_names <- paste0("sigma_species[", species_index, "]")
  sigma <- posterior[draw_ids, sigma_names, drop = FALSE]
  
  # Correct for the mean parameterisation of the lognormal model.
  log_location <- log_mu - 0.5 * sigma^2
  residuals <- sweep(
    -log_location, 2, log(dt[[response_col]]), FUN = "+"
  )
  
  dt[, log_residual := apply(residuals, 2, median)]
  dt[, standardised_residual := apply(residuals / sigma, 2, median)]
  
  plotted_residuals <- if (standardised) residuals / sigma else residuals
  dt[, residual_lower := apply(
    plotted_residuals, 2, quantile, probs = 0.05, names = FALSE
  )]
  dt[, residual_upper := apply(
    plotted_residuals, 2, quantile, probs = 0.95, names = FALSE
  )]
  
  dt[, predictor_value := get(variable)]
  dt[, predictor_plot := if (log_x) {
    log1p(predictor_value)
  } else {
    predictor_value
  }]
  dt[, residual_plot := if (standardised) {
    standardised_residual
  } else {
    log_residual
  }]
  
  p <- ggplot2::ggplot(
    dt, ggplot2::aes(predictor_plot, residual_plot)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = 2, colour = "grey45"
    ) +
    ggplot2::geom_point(alpha = 0.20, size = 0.8) +
    ggplot2::geom_rug(sides = "b", alpha = 0.08) +
    ggplot2::facet_wrap(~ species) +
    ggplot2::labs(
      x = if (log_x) paste0("log(1 + ", variable, ")") else variable,
      y = if (standardised) "Standardised log residual" else "Log residual",
      subtitle = "Posterior median residuals"
    ) +
    ggplot2::theme_bw()
  
  if (smooth) {
    # Only smooth species with enough observations and distinct x values.
    eligible <- dt[
      , .(n = .N, n_x = data.table::uniqueN(predictor_plot)),
      by = species
    ][n >= 10 & n_x >= 5, species]
    
    p <- p + ggplot2::geom_smooth(
      data = dt[species %in% eligible],
      method = "loess",
      formula = y ~ x,
      span = span,
      se = FALSE,
      colour = "steelblue",
      linewidth = 0.9
    )
  }
  
  print(p)
  
  invisible(list(
    data = dt,
    plot = p,
    variable = variable,
    posterior_draw_ids = draw_ids
  ))
}
