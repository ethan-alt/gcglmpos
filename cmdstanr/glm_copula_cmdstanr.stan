functions{ 
  real qnorm_stan(real p) {
    real r;
    real val;
    real q = p - 0.5;
    
    if (fabs(q) <= .425) {
      r = .180625 - q * q;
      val = q * (((((((r * 2509.0809287301226727 +
      33430.575583588128105) * r + 67265.770927008700853) * r +
      45921.953931549871457) * r + 13731.693765509461125) * r +
      1971.5909503065514427) * r + 133.14166789178437745) * r +
      3.387132872796366608)
      / (((((((r * 5226.495278852854561 +
      28729.085735721942674) * r + 39307.89580009271061) * r +
      21213.794301586595867) * r + 5394.1960214247511077) * r +
      687.1870074920579083) * r + 42.313330701600911252) * r + 1.0
      );
    }
    else { /* closer than 0.075 from {0,1} boundary */
      if (q > 0) 
        r = 1.0 - p;
      else r = p;
      
      r = sqrt(-log(r));
      if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
        r += -1.6;
        val = (((((((r * 0.00077454501427834140764 +
        .0227238449892691845833) * r + .24178072517745061177) *
        r + 1.27045825245236838258) * r +
        3.64784832476320460504) * r + 5.7694972214606914055) * r + 4.6303378461565452959) * r +
        1.42343711074968357734) / (((((((r *
        0.00000000105075007164441684324 + 0.0005475938084995344946) *
        r + .0151986665636164571966) * r +
        .14810397642748007459) * r + .68976733498510000455) *
        r + 1.6763848301838038494) * r +
        2.05319162663775882187) * r + 1.);
      }
      else { /* very close to  0 or 1 */
        r += -5.;
        val = (((((((r * 0.000000201033439929228813265 +
        0.0000271155556874348757815) * r +
        .0012426609473880784386) * r + .026532189526576123093) *
        r + .29656057182850489123) * r +
        1.7848265399172913358) * r + 5.4637849111641143699) *
        r + 6.6579046435011037772) / (((((((r *
        0.00000000000000204426310338993978564 + 0.00000014215117583164458887)*
        r + 0.000018463183175100546818) * r +
        0.0007868691311456132591) * r + .0148753612908506148525)
        * r + .13692988092273580531) * r +
        .59983220655588793769) * r + 1.);
      }
      if(q < 0.0) 
        val = -val;
    }
    return val;
  }
  vector qnorm_stan_vec (vector p) {
    int N = num_elements(p);
    vector[N] out;
    for (n in 1:N)
    out[n] = qnorm_stan(p[n]);
    return out;
  }
  real gaussian_copula_cholesky_lpdf(
    matrix U, matrix L
    ) {
      int N                  = rows(U);
      int J                  = cols(U);
      matrix[N,J] Q          = to_matrix( qnorm_stan_vec(to_vector(U)), N, J );
      matrix[J,J] Gammainv   = chol2inv(L);
      return -N * sum(log(diagonal(L))) - 0.5 * sum( add_diag(Gammainv, -1.0) .* crossprod(Q) );
   }
    matrix get_bounds_bernoulli(int[] y, vector p) {
      int N = size(y);
      matrix[size(y), 2] bounds;
      vector[N] onemp = 1.0 - p;
      bounds[, 1] = rep_vector(0.0, N);
      bounds[, 2] = rep_vector(1.0, N);
      for ( i in 1:N ) {
        if ( y[i] == 0 )
          bounds[i,2] = onemp[i];
        else
          bounds[i,1] = onemp[i];
      }
      return bounds;
    }
    matrix get_bounds_poisson(int[] y, vector yp1, vector lambda) {
      int N = size(y);
      matrix[size(y), 2] bounds;
      bounds[,2] = gamma_q(yp1, lambda);
      bounds[,1] = rep_vector(0.0, N);
      for ( i in 1:N ) {
        if(y[i] > 0)
          bounds[i,1] = gamma_q(y[i], lambda[i]);
      }
      return bounds;
    }
}
data {
  int<lower=0> N;                  // number of observations
  int<lower=0> J;                  // total number of outcomes
  int<lower=0> Jn;                 // number of normal outcomes
  int<lower=0> Jb;                 // number of binomial outcomes
  int<lower=0> Jp;                 // number of poisson outcomes
  int<lower=0> K;                  // number of covariates for concatenated design matrices (X1, .., XK)
  matrix[N, Jn]        Yn;         // normal outcomes
  int<lower=0,upper=1> Yb[N, Jb];  // bernoulli outcomes
  int<lower=0>         Yp[N, Jp];  // poisson outcomes
  matrix[N,K]  X;                  // concatenated design matrices (X1, ..., XK)
  int<lower=0> Xindx[J,2];         // Jx2 integer array giving the start and end indexes of X matrix for each outcome
  int<lower=0> Kj[J];              // J-dim integer array giving how many covariates per outcome
  matrix[N,J] mu0;                 // matrix giving prior prediction for each response (for Chen-Ibrahim (2003) conjugate prior of GLM)
  vector[J] lambda;                // precision parameters for conjugate prior of Chen and Ibrhaim (2003)
}
transformed data {
  matrix<lower=1.0>[N,Jp] Ypp1 = to_matrix(Yp);
  vector[K] Xj_mu0j;
  Ypp1 += 1.0;
  for ( j in 1:J ) {
    int start = Xindx[j,1];
    int end   = Xindx[j,2];
    matrix[N,Kj[j]] X_j = X[,start:end];
    Xj_mu0j[start:end] = X_j' * mu0[,j];
  }
}
parameters {
  vector[K] beta;                             // long vector of all regression coefficients
  cholesky_factor_corr[J] L;                  // Cholesky decomposition of JxJ correlation matrix
  real<lower=0> sigmasq[Jn];                  // array of variances
  matrix<lower=0,upper=1>[N,Jb] uraw_b;       // latent variables for bernoulli outcomes
  matrix<lower=0,upper=1>[N,Jp] uraw_p;       // latent variables for Poisson outcomes
}
model {
  // define helper quantities
  int j  = 1;              // counter for outcomes
  int jn = 1;              // counter for normal outcomes
  int jb = 1;              // counter for bernoulli outcomes
  int jp = 1;              // counter for poisson outcomes
  matrix[N,J] U;           // matrix of copula random variables
  vector[N] mu_j;          // vector of means (differs by outcome)
  matrix[N,2] bounds;      // N x 2 matrix giving (lower_bound, upper_bound)
  vector[N] length_bounds; // N-dim vector giving upper_bound - lower_bound
  
  
  // ADD CAUCHY PRIOR FOR NORMAL DISPERSION PARAMETER
  if ( Jn > 0 )
    sigmasq ~ cauchy(0.0, 20.0);
  
  // FOR NORMAL RESPONSES
  //    1. ADD LOG LIKELIHOOD TO TARGET
  //    2. ADD LOG PRIOR TO TARGET
  //    3. CONSTRUCT UNIFORM RANDOM VARIABLE: 
  //       U = PHI(Z); Z = (Y - MU) / SIGMA ~ N(0,1)
  while ( jn <= Jn ) {
    int start = Xindx[j,1];        // start index of design matrix
    int end   = Xindx[j,2];        // end index of design matrix
    matrix[N,Kj[j]] X_j     = X[:,start:end];    // design matrix for jth outcome
    vector[Kj[j]] beta_j    = beta[start:end];   // subset of beta pertaining to j-th outcome
    vector[Kj[j]] y0tilde_j = Xj_mu0j[start:end]; // y_{0j}'X_j
    vector[N]     y_j       = Yn[:,jn];          // jth response variable
    
    // Compute mean of jth outcome
    mu_j = X_j * beta_j;
    
    // normal likelihood
    y_j ~ normal(mu_j, sqrt(sigmasq[j]));
    
    // conjugate linear model prior
    target += lambda[j] * inv(sigmasq[j]) * ( dot_product(y0tilde_j, beta_j) - 0.5 * sum(square(mu_j)) );
    
    // Create uniform random variable
    U[:,j] = Phi_approx( inv_sqrt(sigmasq[j]) * (y_j - mu_j) );
    
    // Increment counters for normal and all outcomes
    jn += 1;
    j  += 1;
  }
  // FOR BERNOULLI RESPONSES
  //    1. SET U[i,j] = LB + (UB - LB) * U_RAW ~ U(LB, UB)
  //    2. JACOBIAN ADJUSTMENT: ADD LOG(UB - LB) TO TARGET
  //    3. ADD LOG CONJUGATE PRIOR FOR REGRESSION COEFFICIENTS
  while ( jb <= Jb ) {
    int start = Xindx[j,1];        // start index of design matrix
    int end   = Xindx[j,2];        // end index of design matrix
    matrix[N,Kj[j]] X_j     = X[:,start:end];               // design matrix for jth outcome
    vector[Kj[j]] beta_j    = beta[start:end];              // subset of beta pertaining to j-th outcome
    vector[Kj[j]] y0tilde_j = Xj_mu0j[start:end];            // y_{0j}'X_j
    
    // compute mean
    mu_j = inv_logit( X_j * beta_j );
    
    // construct uniform random variables and jacboian adjustments
    bounds = get_bounds_bernoulli(Yb[,jb], mu_j);
    length_bounds = bounds[, 2] - bounds[, 1];
    target += log(length_bounds);
    U[,j] = bounds[, 1] + length_bounds .* uraw_b[,jb];
    
    // conjugate logistic regression prior
    target += lambda[j] * ( dot_product(y0tilde_j, beta_j) + sum(log1m(mu_j) ) );
    
    // Increment counters for Bernoulli and all outcomes.
    jb += 1;
    j  += 1;
  }
  // FOR POISSON RESPONSES
  //    1. COMPUTE U = LB + (UB - LB) * U_RAW ~ U(LB, UB)
  //    2. JACOBIAN ADJUSTMENT: ADD LOG(UB - LB) TO TARGET
  //    3. ADD LOG CONJUGATE PRIOR FOR REGRESSION COEFFICIENTS
  while ( jp <= Jp ) {
    int start = Xindx[j,1];        // start index of design matrix
    int end   = Xindx[j,2];        // end index of design matrix
    matrix[N,Kj[j]] X_j     = X[:,start:end];       // design matrix for jth outcome
    vector[Kj[j]] beta_j    = beta[start:end];      // subset of beta pertaining to j-th outcome
    vector[Kj[j]] y0tilde_j = Xj_mu0j[start:end];   // y_{0j}'X_j
    
    // compute mean
    mu_j = exp( X_j * beta_j );                     // mean of jth outcome
    
    // construct uniform random variables and jacboian adjustments
    bounds = get_bounds_poisson(Yp[,jp], Ypp1[,jp], mu_j);
    length_bounds = bounds[, 2] - bounds[, 1];
    U[, j] = bounds[, 1] + length_bounds .* uraw_p[,jp];  // Uniform draw
    target += log(length_bounds);
    
    // conjugate poisson regression prior
    target += lambda[j] * ( dot_product(y0tilde_j, beta_j) - sum(mu_j) );
    
    // Increment counter for Poisson and all outcomes
    jp += 1;
    j  += 1;
  }
  
  // Increase target probability with Gaussian copula density
  L ~ lkj_corr_cholesky(1.0);
  target += gaussian_copula_cholesky_lpdf(U | L);
}
generated quantities {
  corr_matrix[J] Gamma = multiply_lower_tri_self_transpose(L);
}