// DD5 rise, VPD decline, generalized CMD decline with fixed-site SMR interaction.
// Asymmetric moisture response, nutrient quadratic; no nutrient-moisture interaction.
// No BGC effects. DD5 in degree-days, VPD in kPa, CMD and D50 in mm.
functions {
  real si_eta(real T, real V, real D, real n, real m,
              real log_P, real floor_frac, real gap, real v50, real wv,
              real log_d50, real p, real beta, real c_M,
              real opt, real log_r0, real log_r8, real log_rA, real b_NN, real qmax) {
    real smr = 3+3*m;
    real moisture_effect;
    real ed;
    if(smr < opt) moisture_effect = log_r0*square((smr-opt)/opt);
    else moisture_effect = log_r8*square((smr-opt)/(8-opt));
    // q(-1)=log_rA, q(0)=0; same nutrient ratio at every SMR.
    ed = moisture_effect - qmax + (b_NN-log_rA)*n + b_NN*square(n);
    real cmd_loss = 0;
    // Branch avoids log(0) and undefined gradients of 0^p at CMD=0.
    if (D>0) cmd_loss = log(2.0)*exp(p*(log(D)-log_d50-beta*m));
    real ramp=fmin(1, fmax(0,T-200)/(gap*exp(c_M*m)));
    return log_P+ed+log(floor_frac+(1-floor_frac)*ramp)
      +log_inv_logit((v50-V)/wv)-cmd_loss;
  }
}
data {
  int<lower=1> N;
  int<lower=1> S;
  vector<lower=0>[N] y;
  vector<lower=0>[N] DD5;
  vector<lower=0>[N] VPD;
  vector<lower=0>[N] CMD;
  vector<lower=1,upper=5>[N] snr;
  vector<lower=0,upper=8>[N] smr;
  array[N] int<lower=1,upper=S> species;
  vector[S] P_share_centres;
  vector<lower=0>[S] P_share_sd;
  // Column 1: logit(F_m / 50) for the best-site cold floor in metres.
  // Others: log(K_ref-200), log V50, log w_V, log D50, log p, beta.
  int<lower=200,upper=200> cold_anchor_DD5;
  matrix[S,7] prior_centre;
  matrix<lower=0>[S,7] prior_sd;
  // logit(optimum/8), logit(r0), logit(r8), log(rA), nutrient curvature.
  matrix[S,5] ed_prior_centre;
  matrix<lower=0>[S,5] ed_prior_sd;
  int<lower=0,upper=1> prior_only;
  int<lower=1,upper=N> N_check;
  array[N_check] int<lower=1,upper=N> check_idx;
  int<lower=0,upper=N> N_log_lik;
}
transformed data {
  vector[N] n=(snr-rep_vector(3,N))/2;
  vector[N] m=(smr-rep_vector(3,N))/3;
  for(s in 1:S) for(j in 1:5)
    if(ed_prior_sd[s,j]<=0) reject("Edatopic prior SDs must be positive.");
  if(min(y)<=0) reject("SI must be positive.");
  for(s in 1:S) for(j in 1:7)
    if(prior_sd[s,j]<=0) reject("Prior SDs must be positive.");
}
parameters {
  real delta_P;
  real<lower=0,upper=1> rho_P;
  vector[S] z_P;
  matrix[S,7] z_climate;
  // Cold anchor prior is independent by species; six other parameters share
  // standardized deviations, preserving the supplied marginal prior SDs.
  vector[6] delta_climate;
  vector<lower=0,upper=1>[6] rho_climate;
  vector[5] delta_ed;
  vector<lower=0,upper=1>[5] rho_ed;
  matrix[S,5] z_ed;
  real c_M_bar;
  real<lower=0> tau_c_M;
  vector[S] z_c_M;
  real log_sigma_bar;
  real<lower=0> tau_log_sigma;
  vector[S] z_log_sigma;
}
transformed parameters {
  vector[S] plateau_logit=P_share_centres+P_share_sd .* (rep_vector(sqrt(rho_P)*delta_P,S)+sqrt(1-rho_P)*z_P);
  vector<lower=0,upper=50>[S] F_m;
  vector<lower=0,upper=50>[S] P;
  vector[S] log_P;
  vector[S] nutrient_qmax;
  matrix[S,7] theta;
  vector<lower=200>[S] K_ref;
  vector<lower=0,upper=1>[S] cold_floor;
  vector<lower=0>[S] K_gap;
  vector<lower=0>[S] V50;
  vector<lower=0>[S] w_V;
  vector<lower=0>[S] D50;
  vector<lower=0>[S] p_CMD;
  vector[S] beta_CMD;
  vector[S] c_M=c_M_bar+tau_c_M*z_c_M;
  matrix[S,5] theta_ed;
  vector<lower=0,upper=8>[S] O_M;
  vector[S] log_r0;
  vector[S] log_r8;
  vector[S] log_rA;
  vector<lower=0,upper=1>[S] r0;
  vector<lower=0,upper=1>[S] r8;
  vector<lower=0>[S] rA;
  vector[S] b_NN;
  vector<lower=0>[S] sigma_species=exp(log_sigma_bar+tau_log_sigma*z_log_sigma);
  for(s in 1:S) {
    theta[s,1]=prior_centre[s,1]+prior_sd[s,1]*z_climate[s,1];
    for(j in 2:7)
      theta[s,j]=prior_centre[s,j]+prior_sd[s,j]*(
        sqrt(rho_climate[j-1])*delta_climate[j-1]
        +sqrt(1-rho_climate[j-1])*z_climate[s,j]);
    K_gap[s]=exp(theta[s,2]);
    F_m[s]=50*inv_logit(theta[s,1]);
    P[s]=F_m[s]+(50-F_m[s])*inv_logit(plateau_logit[s]);
    log_P[s]=log(P[s]);
    cold_floor[s]=F_m[s]/P[s]; // derived only; lifting P cannot lift the floor in metres
    K_ref[s]=cold_anchor_DD5+K_gap[s];
    V50[s]=exp(theta[s,3]); w_V[s]=exp(theta[s,4]);
    D50[s]=exp(theta[s,5]); p_CMD[s]=exp(theta[s,6]);
    beta_CMD[s]=theta[s,7];
    for(j in 1:5)
      theta_ed[s,j]=ed_prior_centre[s,j]+ed_prior_sd[s,j]*(
        sqrt(rho_ed[j])*delta_ed[j]+sqrt(1-rho_ed[j])*z_ed[s,j]);
    O_M[s]=8*inv_logit(theta_ed[s,1]);
    log_r0[s]=log_inv_logit(theta_ed[s,2]);
    log_r8[s]=log_inv_logit(theta_ed[s,3]);
    log_rA[s]=theta_ed[s,4];
    r0[s]=exp(log_r0[s]); r8[s]=exp(log_r8[s]); rA[s]=exp(log_rA[s]);
    b_NN[s]=theta_ed[s,5];
    {
      real b=b_NN[s];
      real a=b-log_rA[s];
      // Maximum on continuous n in [-1,1] (including SNR between classes).
      real peak=fmax(b-a,b+a);
      if(b<0 && abs(a)< -2*b) peak= -square(a)/(4*b);
      nutrient_qmax[s]=peak;
    }
  }
}
model {
  delta_P ~ std_normal();
  rho_P ~ beta(2,2);
  z_P ~ std_normal();
  to_vector(z_climate) ~ std_normal();
  delta_climate ~ std_normal();
  rho_climate ~ beta(2,2);
  delta_ed ~ std_normal();
  rho_ed ~ beta(2,2);
  to_vector(z_ed) ~ std_normal();
  c_M_bar ~ normal(0,.15);
  tau_c_M ~ normal(0,.10);
  z_c_M ~ std_normal();
  log_sigma_bar ~ normal(log(.25),.5);
  tau_log_sigma ~ normal(0,.35);
  z_log_sigma ~ std_normal();
  if(!prior_only) {
    vector[N] location;
    vector[N] sd_log;
    for(i in 1:N) {
      int s=species[i];
      location[i]=si_eta(DD5[i],VPD[i],CMD[i],n[i],m[i],log_P[s],
        cold_floor[s],K_gap[s],V50[s],w_V[s],theta[s,5],p_CMD[s],beta_CMD[s],c_M[s],O_M[s],log_r0[s],log_r8[s],log_rA[s],b_NN[s],nutrient_qmax[s])
        -.5*square(sigma_species[s]);
      sd_log[i]=sigma_species[s];
    }
    y ~ lognormal(location,sd_log);
  }
}
generated quantities {
  vector[N_check] log_mu_check;
  vector[N_check] log_y_rep;
  vector[N_check] log_location_check;
  vector[N_check] log_sd_check;
  vector[N_log_lik] log_lik;
  for(j in 1:N_check) {
    int i=check_idx[j];
    int s=species[i];
    log_mu_check[j]=si_eta(DD5[i],VPD[i],CMD[i],n[i],m[i],log_P[s],
      cold_floor[s],K_gap[s],V50[s],w_V[s],theta[s,5],p_CMD[s],beta_CMD[s],c_M[s],O_M[s],log_r0[s],log_r8[s],log_rA[s],b_NN[s],nutrient_qmax[s]);
    log_location_check[j]=log_mu_check[j]-.5*square(sigma_species[s]);
    log_sd_check[j]=sigma_species[s];
    log_y_rep[j]=normal_rng(log_location_check[j],sigma_species[s]);
  }
  if(N_log_lik>0) for(i in 1:N_log_lik) {
    int s=species[i];
    real eta=si_eta(DD5[i],VPD[i],CMD[i],n[i],m[i],log_P[s],
      cold_floor[s],K_gap[s],V50[s],w_V[s],theta[s,5],p_CMD[s],beta_CMD[s],c_M[s],O_M[s],log_r0[s],log_r8[s],log_rA[s],b_NN[s],nutrient_qmax[s]);
    log_lik[i]=lognormal_lpdf(y[i]|eta-.5*square(sigma_species[s]),sigma_species[s]);
  }
}

