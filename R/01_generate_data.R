remove(list = ls())
library(mvtnorm)

Ndatasets = 10000
Nmax      = 650

## get copula functions
source('~/Projects/Paper2/cmdstanr/glm_copula_rfuns.R')


set.seed(123)

corrs     = c(-0.40, -0.20, 0.0, 0.20, 0.40, 0.80)
corrs.lbl = c('Moderately negative', 'Lowly negative', 'Independent', 'Lowly positive', 'Moderately positive', 'Highly positive')
names(corrs) = corrs.lbl
corrs.path = c('mod_neg', 'low_neg', 'independent', 'low_pos', 'mod_pos', 'high_pos')

## Obtain vivli results
setwd('~/Projects/Paper2/Vivli/Results')
  ## covariates
  xsmpl   = readRDS('covariate_samples.rds')
  xnames  = xsmpl$xnames
  distvec = xsmpl$dist
  ## responses
  ysmpl        = readRDS('response_samples.rds')   ## sample from vivli data
  ysmpl.h1     = readRDS('sample_h1.rds')          ## sample from vivli data containing only positive treatment effects
  ysmpl.h0     = readRDS('sample_h0.rds')          ## sample from vivli data excluding treatment effects (null validation prior)
  ynames       = ysmpl$ynames
  formula.list = ysmpl$formula.list
  family.list  = ysmpl$family.list
  
  
## Input frequentist
freqfit = readRDS('frequentist_analysis.rds')

gc.sample.x = function(parms, xnames, distvec, N = Nmax, scale = T) {
  P       = length(xnames)
  ## Extract correlation from parameters
  corr    = parms[grepl('corr\\[', names(parms))]
  corr    = corrvec2corrmat(corr)
  ## Draw sample from Gaussian Copula
  X = pnorm(mvtnorm::rmvnorm(N, sigma = corr))
  colnames(X) = xnames
  ## Loop through X's and perform inverse-CDF
  for( p in 1:P ) {
    ## Get name of covariate, extract parameters, and distribution name
    name_p  = xnames[p]
    parms_p = parms[paste0(name_p, '[', c('mu', 'sigma'), ']')]
    dist_p  = distvec[name_p]
    ## Generate covarate from GC based on proper quantile function
    ## depending on distribution
    if ( dist_p == 'lognormal' )
      X[,p] = qlnorm(X[,p], meanlog = parms_p[1], sdlog = parms_p[2])
    if (dist_p == 'normal')
      X[,p] = qnorm(X[,p], mean = parms_p[1], sd = parms_p[2])
  }
  ## None of the variables should be negative nor exceed 100
  X = ifelse(X < 0, 0, X)
  X = ifelse(X > 100, 100, X)
  ## Center/scale if requested (Default)
  if(scale)
    X = as.data.frame(scale(X))
  X
}


gc.sample.y = function(parms, ynames, formula.list, family.list, data, corr = NULL) {
  J = length(formula.list)
  N = nrow(data)
  ## Extract correlation (or generate if applicable)
  if(is.null(corr)) {
    corr = parms[grepl('corr', names(parms))]
    corr = corrvec2corrmat(corr)
  } else {
    corr = matrix(corr, J, J) + diag(1 - corr, J)
  }
  ## Sample from Gaussian copula
  Y = pnorm( mvtnorm::rmvnorm(N, sigma = corr) )
  colnames(Y) = ynames
  ## Loop through outcomes; extract right parameters; generate outcome
  for ( j in 1:J ){
    yname_j       = ynames[j]
    formula_j     = formula.list[[j]]
    family_j      = family.list[[j]]
    rhs.formula_j = formula_j[-2]
    X_j           = model.matrix(rhs.formula_j, data)
    parms_j       = parms[grepl( paste0(yname_j, '\\['), names(parms))]
    if (family_j$family == 'gaussian') {
      indx.sigmasq = which(grepl('sigmasq', names(parms_j)))
      sigmasq_j    = parms_j[indx.sigmasq]   ## extract sigmasq
      parms_j      = parms_j[-indx.sigmasq]  ## now only regression coeffs
    }
    parms_j = parms_j[paste0(yname_j, '[', colnames(X_j), ']')]  ## convert to beta; make sure order is correct based on X_j
    mean_j  = family_j$linkinv(X_j %*% parms_j)
    if (family_j$family == 'gaussian'){
      Y[,j]        = qnorm(Y[,j], mean = mean_j, sd = sqrt(sigmasq_j))
    } else if (family_j$family == 'binomial') {
      Y[,j] = qbinom(Y[,j], size = 1, prob = family_j$linkinv(X_j %*% parms_j))
    } else {
      Y[,j] = qpois(Y[,j], lambda = family_j$linkinv(X_j %*% parms_j))
    }
  }
  data.frame(Y)
}

## Construct data sets
for ( i in seq_len(Ndatasets) ) {
  ## Create covariates
  parms = as.numeric(xsmpl$samples[i, ])
  names(parms) = colnames(xsmpl$samples)
  xdata = gc.sample.x(parms = parms, xnames = xnames, distvec = distvec, N = Nmax, scale = T)
  
  ## Draw treatment indicator
  xdata$TRT01PN = rbinom(Nmax, size = 1, prob = 0.50)
  
  ## Get responses for POS
  parms = as.numeric(ysmpl$samples[i, ])
  names(parms) = colnames(ysmpl$samples)
  ydata = gc.sample.y(parms, ynames, ysmpl$formula.list, ysmpl$family.list, data = xdata)
  saveRDS(
    list(
      'data'         = list(data = cbind(ydata, xdata)),
      'formula.list' = ysmpl$formula.list,
      'family.list'  = ysmpl$family.list,
      'vprior'       = 'posterior',
      'corr'         = 'posterior'
    ),
    file = file.path('~/Projects/Paper2/vprior/pos', paste0('data_', i, '.rds'))
  )
  
  ## Loop through correlations (for type I error / power)
  for (j in seq_along(corrs)) {
    ## Get fixed correlation
    corr = corrs[j]
    ## Sample from type I error
    parms = as.numeric(ysmpl.h0$samples[i, ])
    names(parms) = colnames(ysmpl.h0$samples)
    ydata = gc.sample.y(parms, ysmpl.h0$ynames, ysmpl.h0$formula.list, ysmpl.h0$family.list, data = xdata, corr = corr)
    saveRDS(
      list(
        'data' = cbind(ydata, xdata),
        'formula.list' = ysmpl.h0$formula.list,
        'family.list'  = ysmpl.h0$family.list,
        'vprior'       = 'type1error',
        'corr'         = corr
      ),
      file = file.path('~/Projects/Paper2/vprior/type1error', corrs.path[j], paste0('data_', i, '.rds'))
    )
    
    ## Sample from power
    parms = as.numeric(ysmpl.h1$samples[i, ])
    names(parms) = colnames(ysmpl.h1$samples)
    ydata = gc.sample.y(parms, ysmpl.h1$ynames, ysmpl.h1$formula.list, ysmpl.h1$family.list, data = xdata, corr = corr)
    saveRDS(
      list(
        'data' = cbind(ydata, xdata),
        'formula.list' = ysmpl.h0$formula.list,
        'family.list'  = ysmpl.h0$family.list,
        'vprior'       = 'power',
        'corr'         = corr
      ),
      file = file.path('~/Projects/Paper2/vprior/power', corrs.path[j], paste0('data_', i, '.rds'))
    )
  }
}