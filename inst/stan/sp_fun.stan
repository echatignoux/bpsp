/*------------------------------------------------------------------------------------------------*/
/*  Stan functions for spatial random effects priors and Poisson-bivariate density                */
/*                                                                                                */
/*  Spatial random effects code was adapter from Mitzi Morris                                     */
/*   (see https://mc-stan.org/users/documentation/case-studies/icar_stan.html)                    */
/*                                                                                                */
/*  Coding of Bivariate-Poisson distribution was inspired by code used in the footBayes R package */
/*   (see https://github.com/LeoEgidi/footBayes/tree/master)                                      */
/*------------------------------------------------------------------------------------------------*/

/*================================*/
/* Functions for spatial  analysis*/
/*================================*/

/**
 * Log probability of the intrinsic conditional autoregressive (ICAR) prior
 * Soft constrains has an exponetial decay with size of the area (at around 0.01 for small 
 * areas; 0.001 for large). This (empirically) allows better performances.
 * 
 * @param icar Vector of parameters for spatial smoothing (on unit scale)
 * @param node1 
 * @param node2
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 *
 * @return Log probability density of ICAR prior up to additive constant
 **/
real icar_normal_lpdf(vector phi, 
		      array[] int node1, 
		      array[] int node2, 
		      array[] int comp_size, 
		      array[] int comp_idx) {
  int n_comp = size(comp_size) ;
  real lp;
  int pos=1;
  lp = -0.5 * dot_self(phi[node1] - phi[node2]);
  for (j in 1:n_comp) {
    if (comp_size[j] > 1) {
      lp += normal_lpdf(sum(phi[segment(comp_idx, pos, comp_size[j])]) | 0, 0.01 * (0.1+2*exp(-0.1*log(2)*comp_size[j])) * comp_size[j]) ;
    } else {
      lp += std_normal_lpdf(phi[segment(comp_idx, pos, comp_size[j])]);
    }
    pos += comp_size[j];
  }
  return lp;
}

/**
 * Combine local and global partial-pooling components into the convolved BYM2 term.
 *
 * @param icar local (spatially autocorrelated) component
 * @param iid global component
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 * @param phi proportion of convolution that is spatially autocorrelated
 * @param scale_factor The scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 *
 * @return BYM2 convolution vector
 */
vector convolve_bym2(vector icar, 
		     vector iid,
		     array[] int comp_size, 
		     array[] int comp_idx,
		     real phi, 
		     vector scale_factor
		     ) {
  int n = size(icar) ;
  int n_comp = size(comp_size) ;
  vector[n] convolution;
  int pos=1;
  for (j in 1:n_comp) {
    if (comp_size[j] == 1) {
      convolution[ segment(comp_idx, pos, comp_size[j]) ] = iid[ segment(comp_idx, pos, comp_size[j]) ];
    } else {
      convolution[ segment(comp_idx, pos, comp_size[j]) ] =
	sqrt(phi) * inv(sqrt(scale_factor[j])) * icar[ segment(comp_idx, pos, comp_size[j]) ] +
	sqrt(1 - phi) * iid[ segment(comp_idx, pos, comp_size[j]) ]
	;
    }
    pos += comp_size[j];
  }
  return convolution;
}


/**
 * Scale underliing icar graphe
 *
 * @param icar unscaled icar
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 * @param scale_factor Scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 * @return Vector of scaled re
 **/
vector scale_icar(vector icar,
		  array[] int comp_size, 
		  array[] int comp_idx,
		  vector scale_factor
		  ) {
  int n = size(icar) ;
  int n_comp = size(comp_size) ;
  vector[n] sc_icar;
  int pos=1;
  for (j in 1:n_comp) {
    if (comp_size[j] == 1) {
      sc_icar[ segment(comp_idx, pos, comp_size[j]) ] = icar[ segment(comp_idx, pos, comp_size[j]) ];
    } else {
      sc_icar[ segment(comp_idx, pos, comp_size[j]) ] =
	inv(sqrt(scale_factor[j])) * icar[ segment(comp_idx, pos, comp_size[j]) ] ;
    }
    pos += comp_size[j];
  }
  return sc_icar;
}
 
/**
 * Setting up random effects
 *
 * @param re_typ type of re (1 = iid, 2 = icar, 3 = bym)
 * @param icar local (spatially autocorrelated) component
 * @param iid  Unstructured component
 * @param phi part of spatially autocorrelated coef in the BYM formulation
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 * @param scale_factor Scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 * @return Vector of random effects
 **/
vector set_re(
	      int re_typ,
	      vector icar,
	      vector iid,
	      real phi,
	      array[] int comp_size,
	      array[] int comp_idx,
	      vector scale_factor
	      ) {
  int n_geo = sum(comp_size) ;
  int n_comp = size(comp_size) ;
  vector[n_geo] re ;
  /* re_k */
  if (re_typ == 0 ){
    re =  rep_vector(0,n_geo);
  }
  if (re_typ == 1 ){
    re =  iid ;
  }
  if (re_typ == 2 ){
    re = scale_icar( icar, comp_size, comp_idx, scale_factor);
  }
  if (re_typ == 3 ){
    re = convolve_bym2(icar, iid, comp_size, comp_idx,
                       phi, scale_factor);
  }
  return(re);
}

/**
 * Same as set_re, but dealing with filtering the ICAR random
 * effects, defined on the full graph structure, to the registry area
 *
 * @param re_typ type of re (1 = iid, 2 = icar, 3 = bym)
 * @param icar local (spatially autocorrelated) component
 * @param iid  Unstructured component
 * @param phi part of spatially autocorrelated coef in the BYM formulation
 * @param comp_size number of observational units in each group
 * @param comp_idx index of observations in order of their group membership
 * @param scale_factor Scaling factor for the ICAR variance (see scale_c R function, using R-INLA); 
 * @param geo_ra Index of RA locations amoung all geographical units
 * @return Vector of random effects
 **/

vector set_re_ra(
		 int re_typ,
		 vector icar,
		 vector iid,
		 real phi,
		 array[] int comp_size,
		 array[] int comp_idx,
		 vector scale_factor,
		 array[] int geo_ra
		 ) {
  int n_geo = sum(comp_size) ;
  int n_comp = size(comp_size) ;
  int n_geo_ra = size(geo_ra) ;
  vector[n_geo_ra] re ;
  vector[n_geo] iid_filled ;
  
  /* re_k */
  if (re_typ == 0 ){
    re =  rep_vector(0,n_geo_ra);
  }
  if (re_typ == 1 ){
    re =  iid ;
  }
  if (re_typ == 2 ){
    re = scale_icar( icar,  comp_size, comp_idx, scale_factor)[geo_ra];
  }
  if (re_typ == 3 ){
    iid_filled =  rep_vector(0,n_geo);
    iid_filled[geo_ra] = iid;
    re = convolve_bym2(icar, 
                       iid_filled, 
                       comp_size, 
                       comp_idx,
                       phi, 
                       scale_factor)[geo_ra];
  }
  return(re);
}

/*================================*/
/*     Bivariate-Poisson pdf      */
/*================================*/


/**
 * Log probability of bivariate poisson distribution
 *
 * @param r vector of observations (Y1,Y2)
 * @param lmu_x1 log intensity for x1
 * @param lmu_x2 log intensity for x2
 * @param lmu_x3 log intensity for x3 (common cases)
 *
 * @return Log probability density of bivariate poisson distribution
 **/
real bipois_lpmf(array[] int r , real lmu_x1, real lmu_x2, real lmu_x3) {
  real ss;
  real log_s;
  real lmus;
  int  miny;
  
  miny = min(r[1], r[2]);
  
  ss = poisson_log_lpmf(r[1] | lmu_x1) + poisson_log_lpmf(r[2] | lmu_x2) - exp(lmu_x3);
  if(miny > 0) {
    lmus = -lmu_x1-lmu_x2+lmu_x3;
    log_s = ss; 
    /* log_s = 0; */
    
    for(k in 1:miny) {
      log_s +=
	lmus +
	log(r[1] - k + 1) +
	log(r[2] - k + 1) +
	- log(k);
      ss = log_sum_exp(ss, log_s);
    }
  }
  return(ss);
}

/**
 * Laplace approximation to the log probability of bivariate poisson distribution
 *
 * @param r vector of observations (Y1,Y2)
 * @param lmu_x1 log intensity for x1
 * @param lmu_x2 log intensity for x2
 * @param lmu_x3 log intensity for x3 (common cases)
 *
 * @return Log probability density of bivariate poisson distribution
 **/
real bipoislap_lpmf(array[] int r , real lmu_x1, real lmu_x2, real lmu_x3) {
  real ss;
  real z;
  real mc;
  real ms;
  real logbp;
  real b;
  real delta;
  real xx;
  real yy;
  int  miny;
  miny = min(r[1], r[2]);

  logbp=0;
  z = exp(-lmu_x1-lmu_x2+lmu_x3) ;
  ss = poisson_log_lpmf(r[1] | lmu_x1) + poisson_log_lpmf(r[2] | lmu_x2) - exp(lmu_x3);
  xx = r[1]+.5;
  yy = r[2]+.5;
  b = -(xx+yy)*z-1 ;
  delta = b^2-4*z*(xx*yy*z-0.5);
  mc = (-b - sqrt(delta))/(2*z);
  ms = abs(trigamma(r[1]-mc+1)+trigamma(r[2]-mc+1)+trigamma(mc+1));
  
  logbp = 0.5*log(2*pi()/ms) + 
    - lgamma(r[1] - mc + 1) +
    - lgamma(r[2] - mc + 1) - lgamma(mc + 1) +
    mc*log(z)  + 
    lgamma(r[1]+1) + lgamma(r[2]+1); 
  
  logbp +=  ss + log(normal_cdf( miny | mc, sqrt(ms)) - normal_cdf( 0 | mc, sqrt(ms))) ;
  
  return(logbp);
}


/**
 * Log probability of bivariate poisson distribution
 * combinaison of exact (min(Y1,Y2)<50) and laplace
 * approximation (min(Y1,Y2)>=50)
 *
 * @param r vector of observations (Y1,Y2)
 * @param lmu_x1 log intensity for x1
 * @param lmu_x2 log intensity for x2
 * @param lmu_x3 log intensity for x3 (common cases)
 *
 * @return Log probability density of bivariate poisson distribution
 **/

real bp_lpmf(array[] int r , real lmu_x1, real lmu_x2, real lmu_x3) {
  real ll;
  int  miny;
  miny = min(r[1], r[2]);
  if ((miny > 50)) {
    ll = bipoislap_lpmf(r|lmu_x1,lmu_x2,lmu_x3);
  } else {
    ll = bipois_lpmf(r|lmu_x1,lmu_x2,lmu_x3);
  }
  return(ll);
}


/*================================*/
/*          MISC                  */
/*================================*/


/**
 * Convert logit-scale RE to log-scale
 *
 * @param w : mean Se or PPV
 * @param u : spatial Se or PPV
 */
 vector nuj(real w, vector u) {// Convert logit risk to risk (1st order Taylor)
    vector[size(u)] nu;
    nu = inv(w+(1-w)*exp(-u));
    return nu;
  }

/**
 * Create vector of coef including eventual penalisation
 * to use for multiplication with the W model matrix of covariates
 * W regroup all the covariates (thus might not concern each component i, se and ppv)
 *
 * @param n_w : number of covariates in the matrix 
 * @param d_W : n col of W
 * @param d_w : d_w[k] = index for W columns related to covariate k
 * @param b_wf : unpenalized coefs (size d_W)  
 * @param s_w : penalisaiton variance
 * @param has_w : is covariate k in the 
 * @param Z spline basis
 * @param sig_nl spline sd (i.e. penalty) for nl part
 *
 * @return coefficients to apply to W
 */
vector set_beta(
		vector b_wf, // unpenalized coefs
		array[] int d_w, //n col covar 
		array[] int has_w, // has the covar?
		array[] int is_pen_w, // is the covar penalized?
		array[] real s_w // penalisation if penalized
		) {// Set covar basis
                 
  int d_W = size(b_wf);
  int n_w = size(d_w);
  vector[d_W] b_w;
  int pos_w=1;
  b_w = rep_vector(0,d_W);
  for (k in 1:n_w){// loop over covars
    if (has_w[k]) {
      b_w[pos_w:d_w[k]]=b_wf[pos_w:d_w[k]];// unpenalized fixed and linear part of spline
      if (d_w[k]>1 && is_pen_w[k]) //penalize if needed 
	b_w[(pos_w+1):d_w[k]] =  s_w[k] * b_wf[(pos_w+1):d_w[k]] ;
      pos_w = pos_w + d_w[k] ;
    }
  }
  return b_w;
}

