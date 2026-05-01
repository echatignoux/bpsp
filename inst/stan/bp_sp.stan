functions {
#include sp_fun.stan
}

data{
  int<lower=0> prior_only;  
  int<lower=1> n_obs;       // total number of observations
  array[n_obs,2] int inc;         // response
  vector[n_obs] off;        // offset (log scale)
  // geo effects
  int<lower=1> n_geo;          // total number of geo
  array[n_obs] int<lower=1> id_geo;  // grouping indicator per observation
  // Age
  int<lower=1> n_age;  // total number of observations
  array[n_obs] int<lower=1> id_age;  // grouping indicator per observation
  int d_z; // Dimension of spline basis
  matrix[n_age, d_z] Z;  // Spline basis
  // Ecological covariates (W) 
  int n_w; // Number of covar
  int d_W; // Dim of covar basis
  array[n_w] int d_w; // Dimension of each covar basis
  array[n_w,3] int has_w; // Idensitication of which covars for I, Se and PPV
  array[n_w] int is_pen_w; // Are covars pernalized 
  matrix[n_geo, d_W] W;  // Covar basis
  // Geography over the whole territory
  int<lower=1> n_comp;
  vector[n_comp] scales;
  array[n_comp] int comp_size; // observational units per group
  array[n_geo] int comp_idx; // index of geo, ordered by group
  int<lower=1> n_edges; 
  array[n_edges] int<lower=1, upper=n_obs> node1;
  array[n_edges] int<lower=1, upper=n_obs> node2;
  // Registry area indexes
  int n_geo_ra; 
  int n_geo_ora; 
  array[n_geo_ra] int<lower=1> geo_ra;  // grouping indicator per observation
  array[n_geo_ora] int<lower=1> geo_ora;  // grouping indicator per observation
  int n_obs_ra; 
  int n_obs_ora;
  array[n_obs_ra] int<lower=1> id_geo_ra;  // grouping indicator per observation
  array[n_obs_ora] int<lower=1> id_geo_ora;  // grouping indicator per observation
  array[n_obs_ra] int id_obs_ra; 
  array[n_obs_ora] int id_obs_ora;
  // Form of random effect
  array[3] int re_typ;
  // Restrict to RA data ?
  int ra_only;
  // Penalty for age spline?
  int is_pen;
  // Correlated ui?
  int is_cor;
  //Prior ru
  array[2] real p_ru; // index of geo, ordered by group

}

transformed data {
  array[1] int d_x;
  int has_cov;
  array[3] int n_w_isp; // NUmber of cov for I, se and PPV
  int n_pen_w; // NUmber of penalized covars
  d_x[1]=d_z+1;
  matrix[n_age, d_x[1]] X;  // Spline basis + intercept
  vector[n_age] w_pop; 
  
  array[1] int is_pen_z;
  array[1] int has_z;
  has_z[1] =1;
  is_pen_z[1] =is_pen;

  X = append_col(rep_vector(1,n_age),Z);
  for (i in 1:n_age)
    w_pop[i] = sum(exp(off[(n_geo*(i-1)+1):(n_geo*i)]))/sum(exp(off));
  
  has_cov = d_W>0;
  n_pen_w = sum(is_pen_w); 
  for (i in 1:3)
    n_w_isp[i] = sum(has_w[,i]);
}

parameters{
  array[3] real int_isp ;
  /*Params re*/
  vector[re_typ[1] > 1 ? n_geo : 0] c_icar; /* Structured cj*/
  vector[re_typ[1] != 2 ? n_geo : 0] c_iid; /* Unstructured b1j*/    
  vector[re_typ[2] > 1 ? n_geo : 0] u1_icar; /* Structured b1j*/
  vector[re_typ[2] != 2 ? n_geo_ra : 0] u1_iid; /* Unstructured b1j*/   
  vector[re_typ[3] > 1 ? n_geo : 0] u2_icar; /* Structured b1j*/
  vector[re_typ[3] != 2 ? n_geo_ra : 0] u2_iid; /* Unstructured b1j*/   
  real<lower=0, upper = 2> s_c;
  array[2] real<lower=0, upper = 2> s_u;
  real<lower=0, upper = 1> phi_c; /*propotion of structured variance*/
  array[2] real<lower=0, upper = 1> phis_u; /*propotion of structured variance*/
  vector[n_geo] d_star ;
  array[(!ra_only) ? 1 : 0] real<lower=0> ga_ora ;
  array[(is_cor) ? 1 : 0] real<lower=0,upper=1> ru;

  /*Params splines*/
  array[(is_pen) ? 1 : 0] real<lower=0> s_f;
  matrix[d_z,3] b_zf;  

  /*Params cov*/
  array[n_pen_w] real<lower=0> s_w;
  matrix[d_W,3] b_wf;  

}

transformed parameters{
  vector[n_geo] c ;
  matrix[n_geo_ra,2] u ;
  array[2] real<lower=0, upper=1> w_sp ;
  array[(!ra_only) ? 1 : 0] real<lower=0> s_ora ;
  array[1] real<lower=0,upper=1> r_u;
  matrix[d_x[1],3] b_f;
  matrix[n_age,3] r_a ; // i, se, vpp
  matrix[d_W,3] b_w;
  matrix[n_geo,3] r_w ; // I, se, vpp

  if (is_cor) {
    r_u[1] = ru[1];
  } else {
    r_u[1] = 0;
  }

  // splines for age effects for i, se, vpp
  for (j in 1:3){
    b_f[1,j] = int_isp[j];
    b_f[2:d_x[1],j] = set_beta(b_zf[,j],rep_array(d_z,1),rep_array(1,1),is_pen_z, s_f);
    r_a[,j] = X * b_f[,j];
  }

  // covars for i, se, vpp
  for (j in 1:3){
    if (has_cov) {
      b_w[,j] = set_beta(b_wf[,j],d_w,has_w[,j],is_pen_w,s_w);
      r_w[,j] = W * b_w[,j];
    } else {
      r_w[,j] = rep_vector(0,n_geo);
    }
  }
    

  // mean se and vpp
  w_sp[1] = sum( w_pop .* exp(r_a[,1] ) .* inv_logit(r_a[,2]) ) / sum( w_pop .* exp(r_a[,1]) ) ; // Mean se
  w_sp[2] = sum( w_pop .* exp(r_a[,1] ) .* inv_logit(r_a[,2]) ) / sum( w_pop .* exp(r_a[,1]) .* inv_logit(r_a[,2]) ./ inv_logit(r_a[,3] ) ) ; // Mean vpp

  if (!ra_only){
    s_ora[1] = ga_ora[1] ;
  }
  
  // Set re
  // c
  c = set_re(re_typ[1], c_icar, c_iid, phi_c, comp_size, comp_idx, scales);
  // u1 and u2
  u[,1] = set_re_ra(re_typ[2], u1_icar, u1_iid, phis_u[1], comp_size, comp_idx, scales, geo_ra);
  u[,2] = set_re_ra(re_typ[3], u2_icar, u2_iid, phis_u[2], comp_size, comp_idx, scales, geo_ra);
  u[,2] = r_u[1]*u[,1] + sqrt(1-r_u[1]^2)*u[,2]; // Correlated u's

}

model{
  matrix[n_obs_ra,3] isp_ij ; // log-Incidence, logit-Se, PPV
  matrix[n_obs_ra,3] llbd ; // FN, FP and TP log mean
  vector[(!ra_only) ? n_obs_ora:0] lmu;//mu2

  // BP likelihood in RA
  isp_ij[,1] = off[id_obs_ra] + r_a[id_age[id_obs_ra],1] + s_c * c[geo_ra][id_geo_ra] + r_w[geo_ra,1][id_geo_ra];
  for (j in 2:3){
    isp_ij[,j]= inv_logit(r_a[id_age[id_obs_ra],j] + s_u[(j-1)]*u[,(j-1)][id_geo_ra] + r_w[geo_ra,j][id_geo_ra] );
  }
  /* log-mean of FN, FP and TP*/
  llbd[,1] = isp_ij[,1] + log(1 - isp_ij[,2]);
  llbd[,2] = isp_ij[,1] + log(1 - isp_ij[,3]) + log(isp_ij[,2]) - log(isp_ij[,3]);
  llbd[,3] = isp_ij[,1] + log(isp_ij[,2]);

  // Poisson likelihood in ORA
  // log(mu^P)=log(r^I)+log(SE_i/PPV_i)+dj 
  // + spatial covariates
  if (!ra_only){
    lmu = off[id_obs_ora] + 
      r_a[id_age[id_obs_ora],1]+
      log(inv_logit(r_a[id_age[id_obs_ora],2]+r_w[geo_ora,2][id_geo_ora])) -
      log(inv_logit(r_a[id_age[id_obs_ora],3]+r_w[geo_ora,3][id_geo_ora])) +
      s_ora[1] * c[geo_ora][id_geo_ora] +
      r_w[geo_ora,1][id_geo_ora] ;
  }

  if ( re_typ[1] > 1 )
    c_icar ~ icar_normal(node1, node2, comp_size, comp_idx);
  if ( re_typ[1] != 2 )
    target += std_normal_lpdf(c_iid);
  if ( re_typ[2] > 1 ) {
    u1_icar ~ icar_normal(node1, node2, comp_size, comp_idx);
    target += normal_lpdf(sum(u1_icar[geo_ra]) | 0, 0.005 * n_geo_ra) ;//soft sum to 0 cstr on RA
  }
  if ( re_typ[2] != 2 ) {
    target += std_normal_lpdf(u1_iid);
  }
  if ( re_typ[3] > 1 ) {
    u2_icar ~ icar_normal(node1, node2, comp_size, comp_idx);
    target += normal_lpdf(sum(u2_icar[geo_ra]) | 0, 0.005 * n_geo_ra) ;//soft sum to 0 cstr on RA
  }
  if ( re_typ[3] != 2 ) {
    target += std_normal_lpdf(u2_iid);
  }  // target += beta_lpdf( r_u | 1, 1 ) ;
  if (is_cor)
    // target += normal_lpdf( ru | .5, .5 ) ;
    target += beta_lpdf( ru | p_ru[1], p_ru[2] ) ;
  target += std_normal_lpdf(d_star);
  for (j in 1:3){
    target += std_normal_lpdf(b_zf[,j]);
    target += std_normal_lpdf(b_wf[,j]);
  }
  if (is_pen)
    target += student_t_lpdf(s_f | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(s_w | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += normal_lpdf(int_isp|0,5);
  // target += beta_lpdf( phi | 0.2, 0.2 ) ;
  target += normal_lpdf( phi_c | .5, .5 ) ;
  target += normal_lpdf( phis_u | .5, .5 ) ;
  target += student_t_lpdf(s_c | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(s_u | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  if (!ra_only){
    // target += normal_lpdf(ga_ora|1,.5);
    target += student_t_lpdf(ga_ora | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  }
  // likelihood
  if (!prior_only){
    for (i in 1:n_obs_ra){
      target += bp_lpmf({inc[id_obs_ra[i],1],inc[id_obs_ra[i],2]}|llbd[i,1],llbd[i,2],llbd[i,3]);
      // target += bipois_lpmf({inc[id_obs_ra[i],1],inc[id_obs_ra[i],2]}|llbd[i,1],llbd[i,2],llbd[i,3]);
    }
    if (!ra_only)
      target += poisson_log_lpmf(inc[id_obs_ora,2]|lmu);
  }
}

generated quantities{
  real s_e;
  real s_d;
  real s_ds;
  array[2] real<lower=0> s_nu; 
  real kappa;
  real psi;
  matrix[n_age,4] rho;
  matrix[n_geo,2] theta;
  vector[n_geo_ra] val_theta;
  matrix[n_geo_ra,2] nu;
  array[2] real<lower=0,upper=1> phi_u; //, upper=1
  phi_u[1] = phis_u[1];
  phi_u[2] = r_u[1]^2*phis_u[1]+(1-r_u[1]^2)*phis_u[2];
 
  // Variance of nuj (order 3 DL)
  for (i in 1:2){
    s_nu[i] = sqrt(
		   (1-w_sp[i])^2*s_u[i]^2 + 
		   ((1-w_sp[i])^2*w_sp[i]*(2*w_sp[i]-1) + (w_sp[i]*(1-w_sp[i]))^2/2)*s_u[i]^4+
		   15*((1-w_sp[i])*w_sp[i]*(2*w_sp[i]-1)/6)^2*s_u[i]^6
		   );
  }
 
  // Variance of P and P errors
  s_e = sqrt(s_nu[1]^2+s_nu[2]^2-2*r_u[1]*s_nu[1]*s_nu[2]);  
  s_d = sqrt(s_c^2+s_e^2);  
  if (!ra_only) 
    s_ds = ga_ora[1]; 
  else 
    s_ds = 0;
  //Shrinkage and inflation factors
  kappa = (s_c/s_d)^2;   //Shrinkage factor for c_j ora
  psi = s_ds/s_d; //Variance inflation outside RA

  //Age rates (I, P, se, ppv)
  rho[,1] = exp(r_a[,1]);
  rho[,2] = exp(r_a[,1]) .* (inv_logit(r_a[,2])./inv_logit(r_a[,3]));
  for (j in 1:2) {
    rho[,j+2] = inv_logit(r_a[,j+1]);
  }

  // spatial se and vpp
  for (j in 1:2) {
    nu[,j] = nuj(w_sp[j],s_u[j]*u[,j]+r_w[geo_ra,j+1]);
  }

  // spatial risks for I and P
  theta[geo_ra,1]  =  exp(s_c * c[geo_ra] + r_w[geo_ra,1]);
  theta[geo_ra,2]  =  theta[geo_ra,1] .* nu[,1] ./ nu[,2] ;
  if (!ra_only){
    theta[geo_ora,1]  =  exp( kappa * s_ds * c[geo_ora] + sqrt(kappa * (1-kappa) ) * s_ds * d_star[geo_ora] + r_w[geo_ora,1] ) ;
    theta[geo_ora,2]  =  exp( s_ds * c[geo_ora] + r_w[geo_ora,1] + 
                              log(nuj(w_sp[1],r_w[geo_ora,2])) - 
                              log(nuj(w_sp[2],r_w[geo_ora,3])) ) ;
  }
  // spatial risks for I predicted from Theta2 in RA
  val_theta  =  s_c * c[geo_ra] + 
    log(nuj(w_sp[1],s_u[1]*u[,1])) - 
    log(nuj(w_sp[2],s_u[2]*u[,2])); //dj
  val_theta  =  exp( kappa * val_theta  + // shrunken d_j
		     sqrt(kappa*(1-kappa)) * s_d * d_star[geo_ra] + // + Error for accurate variance
		     r_w[geo_ra,1]) ; // + covariates
}

