

## required packages
library(cmdstanr)
library(coda)
library(posterior)
ncores = max(1, parallel::detectCores() - 1)
options(mc.cores = ncores)

## read stan model
setwd('~/Projects/Paper2/cmdstanr')
mod = cmdstan_model('glm_copula_cmdstanr.stan')
#' Bayesian copula GLM using stan
#' 
#' Samples from the posterior distribution of a Gaussian copula GLM.
#' Uses the conjugate prior of Chen and Ibrahim (2003) for regression
#' coefficients, a cauchy(20) prior for dispersion parameters, and
#' a LKJ(1.0) prior for the correlation matrix. Currently, only normal,
#' logistic, and Poisson (log link) models are supported.
#' 
#' @param formula.list list of formulas
#' @param family.list list of families
#' @param data the data set
#' @param mu0 a matrix giving a prior prediction parameter for conjugate prior. 
#' If NULL, default values are 0, 0.5, and 1 for normal, Bernoulli, and Poisson 
#' responses, respectively.
#' @param lambda prior precision parameter for conjugate prior.
#' @param init list of dimension number of chains each itself a named list 
#' giving initial values. Highly recommended to set as 'mle', which uses the 
#' maximum likelihood estimate for initial values of the regression coefficients.
#' @param ... other parameters to pass onto cmdstanr:::sample
#' 
#' @return a list of lists, first element of list is posterior samples.
copula.glm = function(formula.list, family.list, data, mu0 = NULL, lambda = .01, init = 'mle', ...) {
  N  = nrow(data)
  J  = length(formula.list)
  
  ## reorder variables if necessary
  n.indx    = sapply(family.list, function(x) x$family) == 'gaussian'
  b.indx    = sapply(family.list, function(x) x$family) == 'binomial'
  p.indx    = sapply(family.list, function(x) x$family) == 'poisson'
  new.order = c(which(n.indx), which(b.indx), which(p.indx))
  
  ## total number of endpoints and matrix of each outcome type
  Jn = sum(n.indx)
  Jb = sum(b.indx)
  Jp = sum(p.indx)
  Y  = sapply(formula.list, function(x) data[, all.vars(x)[1]])
  Yn = Y[, n.indx, drop = F]
  Yb = Y[, b.indx, drop = F]
  Yp = Y[, p.indx, drop = F]
  
  args = list(...)
  chains = args$chains
  if (is.null(chains) )
    chains = 4  ## stan default
  
  
  if ( !all( new.order == 1:J ) ) {
    message("Rearranging order of dependent variables to normal, binomial, poisson")
    formula.list = formula.list[new.order]
    family.list  = family.list[new.order]
  }
  ## check if two-sided formula
  if ( all( sapply(formula.list, length) != 3 ) )
    stop("Each element in formula.list must be a two-sided formula")
  
  
  ## get design matrices and number of covariates per regression
  Xlist = lapply(formula.list, model.matrix, data = data)
  Xbig  = do.call(cbind, Xlist)
  K     = ncol(Xbig)
  Kj    = sapply(Xlist, ncol)
  
  ## conduct frequentist analysis to get good starting values (if specified)
  if ( init == 'mle' ) {
    fitlist       = mapply(function(x, y) glm(formula = x, family = y, data = data), x = formula.list, y = family.list, SIMPLIFY = F)
    beta.start    = sapply(fitlist, coef)
    sigmasq.start = sapply( fitlist[sapply(family.list, function(x) x$family) == 'gaussian'], function(x) summary(x)$dispersion )
    sigmasq.start = as.array(sigmasq.start)
    
    ## approximate MLE of correlation matrix
    betalist = lapply(fitlist, function(f) coef(f))
    Un = matrix(nrow = N, ncol = Jn)
    Ub = matrix(nrow = N, ncol = Jb)
    Up = matrix(nrow = N, ncol = Jp)
    j = 1
    jn = 1
    jb = 1
    jp = 1
    while ( jn <= Jn ) {
      Un[, jn] = pnorm( Yn[, j], mean = family.list[[j]]$linkinv(Xlist[[j]] %*% betalist[[j]]), sd = sqrt(summary(fitlist[[j]])$dispersion) )
      jn = jn + 1
      j  = j  + 1
    }
    while ( jb <= Jb ) {
      mu = family.list[[j]]$linkinv(Xlist[[j]] %*% betalist[[j]])
      Ub[, jb] = rowMeans(
        cbind(
          min = pbinom(Yb[, jb] - 1, size = 1, prob = mu),
          max = pbinom(Yb[, jb], size = 1, prob = mu)
        )
      )
      jb = jb + 1
      j  = j  + 1
    }
    while ( jp <= Jp ) {
      mu = family.list[[j]]$linkinv(Xlist[[j]] %*% betalist[[j]])
      Up[, jp] = rowMeans(
        cbind(
          min = ppois(Yp[, jp] - 1, lambda = mu),
          max = ppois(Yp[, jp], lambda = mu)
        )
      )
      jp = jp + 1
      j  = j  + 1
    }
    Gamma.start = cor(cbind(Un, Ub, Up))
    L = t(chol(Gamma.start))
    
    init = lapply(1:chains, function(i) 
      list(
        'beta' = unlist(beta.start), 
        'sigmasq' = sigmasq.start, 
        'L' = L,
        'uraw_b' = Ub,
        'uraw_p' = Up
      )
    )
    remove(list = 'Un', 'Ub', 'Up', 'Gamma.start', 'L', 'beta.start')
  }
  
  ## obtain matrix giving start and end values of x matrix per endpoint
  Xindx  = matrix(nrow = J, ncol = 2)
  xstart = 1
  for ( j in 1:J ) {
    xend = xstart + Kj[j] - 1
    Xindx[j, ] = c(xstart, xend)
    xstart = xend + 1
  }
  
  ## Set default value of mu0
  if ( is.null(mu0) ) {
    mu0 = matrix(0, nrow = N, ncol = J)
    for ( j in 1:J ) {
      if(family.list[[j]]$family == 'binomial')
        mu0[, j] = 0.5
      if(family.list[[j]]$family == 'poisson')
        mu0[, j] = 1
    }
  }
  
  ## Set default value of lambda
  if ( is.null(lambda) )
    lambda = rep(0.01, J)
  if ( length(lambda) == 1 )
    lambda = rep(lambda, J)
  if ( length(lambda) != J )
    stop("lambda0 must be a scalar or vector of dimension length(formula.list)")
  
  ## construct stan data
  standat = list(
    'N'      = N,
    'J'      = J,
    'Jn'     = Jn,
    'Jb'     = Jb,
    'Jp'     = Jp,
    'K'      = K,
    'Yn'     = Yn,
    'Yb'     = Yb,
    'Yp'     = Yp,
    'X'      = Xbig,
    'Xindx'  = Xindx,
    'Kj'     = Kj,
    'mu0'    = mu0,
    'lambda' = lambda
  )
  
  ## clean up--remove all un-needed variables to optimize memory usage
  suppressWarnings(
    remove(
      list = c(
        'xstart', 'xend', 'U', 'Y', 'n.indx', 'b.indx', 'p.indx',
        'Gamma.start', 'L', 'N', 'K', 'Yn', 'Yb',
        'Yp', 'Xbig', 'Xindx', 'Kj', 'mu0_Xj', 'lambda'
      )
    )
  )
  ## fit stan model
  fit = mod$sample(data = standat, init = init, ...)
  # return(fit)
  ## subset names of variables we care about
  draws = fit$draws()
  smpl.names = dimnames(draws)$variable
  keep.names = smpl.names[
    ! (grepl('uraw_', smpl.names) | grepl('lp__', smpl.names) | grepl('L\\[', smpl.names) )
  ]
  draws = draws[,,keep.names]
  
  
  ## rename betas for clarity about which equation each comes from, e.g., beta[j][0:(Kj-1)]
  beta.old.indx = which(grepl('beta', keep.names))
  beta.new      = character()
  ynames        = sapply(formula.list, function(x) all.vars(x)[1])
  xnames        = lapply(fitlist, function(x) names(x$coefficients) )
  xnames        = lapply(xnames, function(x) paste0('[', x, ']'))
  for( j in 1:J ) {
    beta.new = c(beta.new, paste0(ynames[j], xnames[[j]]))
  }
  dimnames(draws)$variable[beta.old.indx] = beta.new
  
  ## rename sigmasq for clarity
  sigmasq.indx = which(grepl('sigmasq', dimnames(draws)$variable))
  sigmasq.new = character(length = Jn)
  for ( j in 1:Jn ) {
    sigmasq.new[j] = paste0(all.vars(formula.list[[j]])[1], '[sigmasq]')
  }
  dimnames(draws)$variable[sigmasq.indx] = sigmasq.new
  
  corr.names = dimnames(draws)$variable[grepl('Gamma', dimnames(draws)$variable)]
  corr.indx  = grepl('Gamma', dimnames(draws)$variable)
  corr.names = gsub('Gamma', 'corr', dimnames(draws)$variable[corr.indx])
  dimnames(draws)$variable[corr.indx] = corr.names
  
  ## return list giving draws and names for easy indexing.
  list(
    'samples'      = as_draws_matrix(draws),
    'formula.list' = formula.list,
    'family.list'  = family.list,
    'coef.names'   = beta.new,
    'var.names'    = sigmasq.new,
    'corr.names'   = corr.names,
    'n.gaussian'   = Jn,
    'n.bernoulli'  = Jb,
    'n.poisson'    = Jp,
    'ynames'       = ynames
  )
}

#' Posterior prediction of Gaussian copula
#' 
#' Takes as input a result from calling \code{copula.glm} and a new data set
#' and outputs a data frame giving a posterior predictive sample
#' 
#' @param fit result of calling \code{copula.glm}
#' @param newdata a \code{data.frame} of new data for which predictions
#' are to be made
#' @param indx which sample from \code{fit} should be utilized to draw a 
#' prediction, defaults to first
#' 
#' @return a \code{data.frame} giving a posterior predictive sample (along with
#' covariates)
copula.glm.predict = function(fit, newdata, indx = 1) {
  N  = nrow(newdata)
  J  = length(fit$formula.list)
  Jn = fit$n.gaussian
  Jb = fit$n.bernoulli
  Jp = fit$n.poisson
  if(any(fit$ynames %in% names(newdata)))
    newdata = newdata[,-which(names(newdata) %in% fit$ynames)]
  
  ## Extract correlation matrix
  corr = fit$samples[indx, fit$corr.names]
  corr = matrix(corr, J, J)
  
  ## Generate Gaussian copula random variables
  Y = pnorm( matrix(rnorm(N * J), nrow = N, ncol = J) %*% chol(corr) )
  colnames(Y) = fit$ynames
  
  for ( j in 1:J ) {
    rhs.fmla  = fit$formula.list[[j]][-2]
    X         = model.matrix(rhs.fmla, newdata)
    parms     = fit$samples[indx, ]
    yname     = fit$ynames[j]
    parmnames = colnames(parms)
    parmnames = parmnames[grepl(paste0(yname, '\\['), colnames(parms))]
    parms     = as.vector( parms[,parmnames] )
    names(parms) = parmnames
    if ( fit$family.list[[j]]$family == 'gaussian' ) {
      sigmasq.indx = names(parms) == paste0(yname, '[sigmasq]')
      sigma        = sqrt(parms[sigmasq.indx])
      mean         = family.list[[j]]$linkinv( as.matrix(parms[!sigmasq.indx]) )
      Y[, j]       = qnorm(Y[, j], mean = mean, sd = sigma)
    } else if (fit$family.list[[j]]$family == 'binomial')
      Y[, j] = qbinom(Y[, j], size = 1, prob = family.list[[j]]$linkinv(X %*% as.matrix(parms)))
    else if (fit$family.list[[j]]$family == 'poisson')
      Y[, j] = qpois(Y[, j], lambda = family.list[[j]]$linkinv(X %*% as.matrix(parms)))
    else
      stop("Invalid family")
  }
  cbind(Y, newdata)
}




#' Generates from a Gaussian copula GLM model
#' 
#' Generate Gaussian copula data given means and dispersion parameters
#' 
#' @param mean a \code{matrix} of means
#' @param family.list a \code{list} giving which family each column of \code{mean} belongs to
#' @param disp  a vector giving dispersion parameters; ignored if family is binomial or poisson
#' 
#' @return a \code{matrix} of a draw from a Gaussian copula GLM.
copula.glm.generate = function(
  mean, corr, family.list, disp = rep(1, ncol(corr))
) {
  n = nrow(mean)
  J = ncol(corr)
  U = pnorm( mvtnorm::rmvnorm(n = n, sigma = corr) )
  Y = U
  for ( i in 1:J ) {
    if ( family.list[[i]]$family == 'gaussian' )
      Y[,i] = qnorm( U[, i], mean = mean[, i], sqrt(disp[i]) )
    else if ( family.list[[i]]$family == 'binomial' )
      Y[,i] = qbinom( U[, i], size = 1, prob = mean[, i] )
    else if ( family.list[[i]]$family == 'poisson' )
      Y[,i] = qpois(U[, i], mean[, i])
    else
      return(NA)
  }
  colnames(Y) = paste0('y', 1:J)
  Y
}


#' Compute posterior probabilities and correlation matrix of treatment effect
#' 
#' Given output from a copula.glm object, outputs the posterior probability
#' that each treamtment effect exceeds a user-specified threshold
#' 
#' @param fit output from calling copula.glm
#' @param varname name of covariate corresponding to treatment effect
#' @param nsmpl.normal number of multivariate normal samples to draw to approximate
#' posterior probability (if 0 or negative, uses MCMC draws)
#' @param thresh threshold for treatment success
#' 
#' @return a named vector giving posterior probabilities and correlations
#' between treatment effects
copula.glm.postprobs = function(fit, varname, nsmpl.normal = 0, thresh = 0) {
  J         = length(fit$formula.list)
  keeps     = colnames(fit$samples)[which(grepl(varname, colnames(fit$samples)))]
  smpl      = fit$samples[, keeps]
  smpl.cov  = cov(smpl)
  smpl.cor  = cov2cor(smpl.cov)
  cor.names = matrix(fit$corr.names, J, J)
  if ( nsmpl.normal > 0 ) {
    smpl = mvtnorm::rmvnorm(n = nsmpl.normal, mean = colMeans(smpl), sigma = smpl.cov)
    colnames(smpl) = keeps
  }
  postprobs = colMeans(smpl > thresh)
  smpl.cor  = as.vector( smpl.cor[lower.tri(smpl.cor)] )
  names(smpl.cor) = cor.names[lower.tri(cor.names)]
  c(postprobs, smpl.cor)
}


optfun = function(x, zmax, alpha = 0.05)
  (mean(zmax > x) - alpha)^2

#' Find threshold to declare clinical success
#' 
#' Based on posterior correlation matrix, uses simulation to find 
#' the optimal value to control type I error rate
#' 
#' @param cor correlation matrix
#' @param nsmpl number of samples from null to take
#' @param alpha desired type I error rate
#' @param ... optional arguments to pass onto \code{optimize}
#' 
#' @return vector giving cutoff based on maximum posterior probability given
#' correlation
threshold.gc = function(cor, nsmpl = 200000, alpha = .05, ...) {
  ## Get correlation matrix
  cormat = diag(3)
  cormat[lower.tri(cormat)] = as.numeric(cor)
  cormat[upper.tri(cormat)] = t(cormat)[upper.tri(cormat)]
  
  indx.list = list(1:3, 1:2, c(1, 3), 2:3)
  thresh = numeric(length(indx.list))
  
  ## Simulate from Gaussian Copula
  Z = data.frame( pnorm( mvtnorm::rmvnorm(n = nsmpl, sigma = cormat) ) )
  for ( i in seq_along(indx.list) ) {
    thresh[i] = optimize(
      optfun, interval = c(0.94, 1.00), 
      zmax = do.call(pmax, Z[, indx.list[[i]]]), 
      alpha = alpha
      , ...
    )$minimum
  }
  thresh
}






make_corr_exch = function(rho, J = 3, as.vector = T) {
  res = matrix(rho, J, J) + diag(1 - rho, J)
  if(as.vector)
    return(as.vector(res))
  res
}

corrvec2corrmat = function(corrvec) {
  J   = sqrt(length(corrvec))
  matrix(corrvec, nrow = J, ncol = J)
}

get.pvals = function(formula.list, family.list, data, trtname = 'TRT01PN') {
  fitlist = mapply(
    function(formula, family) { summary(glm(formula, family, data = data)) },
    formula = formula.list, family = family.list,
    SIMPLIFY = F
  )
  normal.indx = sapply(family.list, function(x) x$family == 'gaussian')
  zval = sapply(fitlist, function(x) x$coefficients['TRT01PN', 3])
  df   = sapply(fitlist, function(x) x$df.residual)
  pval = numeric(length(formula.list))
  pval[normal.indx] = pt(zval[normal.indx], df = df[normal.indx], lower.tail = F)
  pval[!normal.indx] = pnorm(zval[!normal.indx], lower.tail = F)
  names(pval) = names(formula.list)
  pval
}








## stan model
mod.2b1n = cmdstan_model('~/Projects/Paper2/twobernoulli_onenormal/glm_copula_opt.stan')

#' Bayesian copula GLM using stan
#' 
#' Samples from the posterior distribution of a Gaussian copula GLM.
#' Uses the conjugate prior of Chen and Ibrahim (2003) for regression
#' coefficients, a cauchy(20) prior for dispersion parameters, and
#' a LKJ(1.0) prior for the correlation matrix. Currently, only normal,
#' logistic, and Poisson (log link) models are supported.
#' 
#' @param formula.list list of formulas
#' @param family.list list of families
#' @param data the data set
#' @param mu0 a matrix giving a prior prediction parameter for conjugate prior. 
#' If NULL, default values are 0, 0.5, and 1 for normal, Bernoulli, and Poisson 
#' responses, respectively.
#' @param lambda prior precision parameter for conjugate prior.
#' @param init list of dimension number of chains each itself a named list 
#' giving initial values. Highly recommended to set as 'mle', which uses the 
#' maximum likelihood estimate for initial values of the regression coefficients.
#' @param ... other parameters to pass onto cmdstanr:::sample
#' 
#' @return a list of lists, first element of list is posterior samples.
copula.glm.2b1n = function(formula.list, family.list, data, lambda = .01, ...) {
  N  = nrow(data)
  J  = 3
  
  ## reorder variables if necessary
  n.indx    = sapply(family.list, function(x) x$family) == 'gaussian'
  b.indx    = sapply(family.list, function(x) x$family) == 'binomial'
  new.order = c(which(b.indx), which(n.indx))
  
  formula.list = formula.list[new.order]
  family.list  = family.list[new.order]
  
  standat        = lapply(formula.list, function(f) model.matrix(f, data = data))
  names(standat) = c('Xb1', 'Xb2', 'Xn')
  
  standat = c(
    list(
      'N'   = N,
      'Kb1' = ncol(standat[['Xb1']]),
      'Kb2' = ncol(standat[['Xb2']]),
      'Kn'  = ncol(standat[['Xn']]),
      'yb1' = data[, all.vars(formula.list[[1]])[1]],
      'yb2' = data[, all.vars(formula.list[[2]])[1]],
      'yn'  = data[, all.vars(formula.list[[3]])[1]]
    ),
    standat,
    list(
      'mu0' = t(replicate(N, c(0.5, 0.5, 0))),
      'lambda' = rep(lambda, 3)
    )
  )
  
  ## Get initial values to be MLE
  fitlist = mapply(
    function(formula, family) 
      glm(formula = formula, family = family, data = data),
    formula = formula.list, family = family.list,
    SIMPLIFY = F
  )
  init = list(
    'beta_b1'   = coef(fitlist[[1]]),
    'beta_b2'   = coef(fitlist[[2]]),
    'beta_n'    = coef(fitlist[[3]]),
    'sigmasq_n' = summary(fitlist[[3]])$dispersion
  )
  mu = predict(fitlist[[1]], type = 'response')
  init$uraw_b1 = rowMeans( cbind(
    pbinom(standat$yb1 - 1, size = 1, prob = mu),
    pbinom(standat$yb1    , size = 1, prob = mu)
  ) )
  mu = predict(fitlist[[2]], type = 'response')
  init$uraw_b2 = rowMeans( cbind(
    pbinom(standat$yb2 - 1, size = 1, prob = mu),
    pbinom(standat$yb2    , size = 1, prob = mu)
  ) )
  init$L = cbind(
    init$uraw_b1, init$uraw_b2, pnorm(standat$yn, mean = standat$Xn %*% init$beta_n, sd = sqrt(init$sigmasq_n))
  )
  init$L = t(chol(cor(qnorm(init$L))))
  
  args = list(...)
  chains = args$chains
  chains = ifelse(is.null(chains), 1, chains)
  init   = lapply(1:chains, function(i) init)
  fit    = mod.2b1n$sample(data = standat, init = init, ...)
  draws = fit$draws()
  smpl.names = dimnames(draws)$variable
  keep.names = smpl.names[
    ! (grepl('uraw_', smpl.names) | grepl('lp__', smpl.names) | grepl('L\\[', smpl.names) )
  ]
  draws = draws[,,keep.names]
  
  ## rename betas for clarity about which equation each comes from, e.g., beta[j][0:(Kj-1)]
  beta.old.indx = which(grepl('beta', keep.names))
  beta.new      = character()
  ynames        = sapply(formula.list, function(x) all.vars(x)[1])
  xnames        = lapply(fitlist, function(x) names(x$coefficients) )
  xnames        = lapply(xnames, function(x) paste0('[', x, ']'))
  for( j in 1:J ) {
    beta.new = c(beta.new, paste0(ynames[j], xnames[[j]]))
  }
  dimnames(draws)$variable[beta.old.indx] = beta.new
  
  ## rename sigmasq for clarity
  sigmasq.indx = which(grepl('sigmasq', dimnames(draws)$variable))
  sigmasq.new  = paste0(all.vars(formula.list[[3]])[1], '[sigmasq]')
  dimnames(draws)$variable[sigmasq.indx] = sigmasq.new
  return(
    list(
      'samples'      = as_draws_matrix(draws),
      'formula.list' = formula.list,
      'family.list'  = family.list,
      'coef.names'   = beta.new,
      'var.names'    = sigmasq.new,
      'corr.names'   = keep.names[grepl('corr\\[', keep.names)],
      'ynames'       = ynames
    )
  )
}


#' Bayesian copula GLM using stan
#' 
#' Samples from the posterior distribution of a Gaussian copula GLM.
#' Uses the conjugate prior of Chen and Ibrahim (2003) for regression
#' coefficients, a cauchy(20) prior for dispersion parameters, and
#' a LKJ(1.0) prior for the correlation matrix. Currently, only normal,
#' logistic, and Poisson (log link) models are supported.
#' 
#' @param formula.list list of formulas
#' @param family.list list of families
#' @param data the data set
#' @param mu0 a matrix giving a prior prediction parameter for conjugate prior. 
#' If NULL, default values are 0, 0.5, and 1 for normal, Bernoulli, and Poisson 
#' responses, respectively.
#' @param lambda prior precision parameter for conjugate prior.
#' @param init list of dimension number of chains each itself a named list 
#' giving initial values. Highly recommended to set as 'mle', which uses the 
#' maximum likelihood estimate for initial values of the regression coefficients.
#' @param ... other parameters to pass onto cmdstanr:::sample
#' 
#' @return a list of lists, first element of list is posterior samples.
copula.glm.2b1n.postprobs = function(formula.list, family.list, data, lambda = .01, trt.indx = 2, method = 'mcmc', clt.nsmpl = 10000, ...) {
  N  = nrow(data)
  J  = 3
  
  ## reorder variables if necessary
  n.indx    = sapply(family.list, function(x) x$family) == 'gaussian'
  b.indx    = sapply(family.list, function(x) x$family) == 'binomial'
  new.order = c(which(b.indx), which(n.indx))
  
  formula.list = formula.list[new.order]
  family.list  = family.list[new.order]
  
  standat        = lapply(formula.list, function(f) model.matrix(f, data = data))
  names(standat) = c('Xb1', 'Xb2', 'Xn')
  
  standat = c(
    list(
      'N'   = N,
      'Kb1' = ncol(standat[['Xb1']]),
      'Kb2' = ncol(standat[['Xb2']]),
      'Kn'  = ncol(standat[['Xn']]),
      'yb1' = data[, all.vars(formula.list[[1]])[1]],
      'yb2' = data[, all.vars(formula.list[[2]])[1]],
      'yn'  = data[, all.vars(formula.list[[3]])[1]]
    ),
    standat,
    list(
      'mu0' = t(replicate(N, c(0.5, 0.5, 0))),
      'lambda' = rep(lambda, 3)
    )
  )
  
  ## Get initial values to be MLE
  fitlist = mapply(
    function(formula, family) 
      glm(formula = formula, family = family, data = data),
    formula = formula.list, family = family.list,
    SIMPLIFY = F
  )
  init = list(
    'beta_b1'   = coef(fitlist[[1]]),
    'beta_b2'   = coef(fitlist[[2]]),
    'beta_n'    = coef(fitlist[[3]]),
    'sigmasq_n' = summary(fitlist[[3]])$dispersion
  )
  mu = predict(fitlist[[1]], type = 'response')
  init$uraw_b1 = rowMeans( cbind(
    pbinom(standat$yb1 - 1, size = 1, prob = mu),
    pbinom(standat$yb1    , size = 1, prob = mu)
  ) )
  mu = predict(fitlist[[2]], type = 'response')
  init$uraw_b2 = rowMeans( cbind(
    pbinom(standat$yb2 - 1, size = 1, prob = mu),
    pbinom(standat$yb2    , size = 1, prob = mu)
  ) )
  init$L = cbind(
    init$uraw_b1, init$uraw_b2, pnorm(standat$yn, mean = standat$Xn %*% init$beta_n, sd = sqrt(init$sigmasq_n))
  )
  init$L = t(chol(cor(qnorm(init$L))))
  
  args = list(...)
  chains = args$chains
  chains = ifelse(is.null(chains), 1, chains)
  init   = lapply(1:chains, function(i) init)
  fit    = mod.2b1n$sample(data = standat, init = init, ...)
  
  ## extract treatment as matrix
  fit = as_draws_matrix( fit$draws( paste0('beta_', c('b1', 'b2', 'n'), '[', trt.indx, ']') ) )
  cor = cor(fit)
  if ( method == 'clt')
    fit = mvtnorm::rmvnorm(n = clt.nsmpl, mean = colMeans(fit), sigma = cov(fit))
  post.probs = colMeans(fit > 0)
  names(post.probs) = sapply(formula.list, function(f) all.vars(f)[1])
  cor = cor(fit)
  cor = cor[lower.tri(cor)]
  names(cor) = c('cor_12', 'cor_13', 'cor_23')
  c(post.probs, cor)
}










