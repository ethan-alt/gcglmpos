remove(list = ls())

library(cmdstanr)
library(coda)
library(rstan)
ncores = parallel::detectCores() - 1
options(mc.cores = ncores)
rstan_options(auto_write = TRUE)


setwd('C:/users/ealt/Desktop')
mod = cmdstan_model('copula_glm_general_cmdstanr.stan')
# mod2 = rstan::stan_model('copula_glm_general_cmdstanr.stan')

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
#' giving inital values. Highly recommended to set as 'mle', which uses the 
#' maximum likelihood estimate for inital values of the regression coefficients.
#' @param stanfit should a stanfit object be returned
#' 
#' @return a list of lists, first element of list is posterior samples.
copula.glm = function(formula.list, family.list, data, mu0 = NULL, lambda = .01, init = 'mle', stanfit = T, ...) {
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
        'Gamma.start', 'L', 'N', 'Jb', 'Jp', 'K', 'Yn', 'Yb',
        'Yp', 'Xbig', 'Xindx', 'Kj', 'mu0_Xj', 'lambda'
      )
    )
  )
  
  ##
  ## OBTAIN NAMES OF PARAMETERS FOR EASY INDEXING
  ##
  
  
  
  ## fit stan model
  fit = mod$sample(data = standat, init = init, ...)
  if(!stanfit)
    return(
      list(
        'samples' = fit,
      )
    )
  
  ## convert to rstan object
  fit = rstan::read_stan_csv(fit$output_files())
  
  ## rename betas for clarity about which equation each comes from, e.g., beta[j][0:(Kj-1)]
  beta.old.indx = which(grepl('beta', names(fit)))
  beta.new      = character()
  ynames        = sapply(formula.list, function(x) all.vars(x)[1])
  xnames        = lapply(fitlist, function(x) names(x$coefficients) )
  xnames        = lapply(xnames, function(x) paste0('[', x, ']'))
  for( j in 1:J ) {
    beta.new = c(beta.new, paste0(ynames[j], xnames[[j]]))
  }
  names(fit)[beta.old.indx] = beta.new
  
  ## rename sigmasq for clarity
  sigmasq.indx = which(grepl('sigmasq', names(fit)))
  sigmasq.new = character()
  for ( j in 1:Jn ) {
    sigmasq.new = paste0(all.vars(formula.list[[j]])[1], '[sigmasq]')
  }
  names(fit)[sigmasq.indx] = sigmasq.new
  
  corr.names = names(fit)[grepl('Gamma', names(fit))]
  corr.names = matrix(corr.names, nrow = 3, ncol = 3)
  corr.names = corr.names[lower.tri(corr.names)]
  
  ## return stan object
  list(
    'samples' = fit,
    'formula.list' = formula.list,
    'family.list'  = family.list,
    'coef.names'   = beta.new,
    'var.names'    = sigmasq.new,
    'corr.names'   = corr.names
  )
}







gc_generate = function(
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

n = 200
J = 3
p = 3
repeat { 
  corr = randcorr::randcorr(J)
  LT = abs( corr[lower.tri(corr)] )
  if ( all(LT <= 0.60) )
    break
}

# expose_stan_functions(mod)

data = matrix(rnorm(n*p, sd = 2), nrow = n, ncol = p)
colnames(data) = paste0('x', 1:p)
data = data.frame(data)
formula.list = list(
  y1 ~ x1,
  y2 ~ x1 + x2,
  y3 ~ x1 + x2 + x3
)
family.list = list(
  gaussian(), binomial(), poisson()
)
rhs.formula.list = lapply(formula.list, function(f) f[-2])
Xlist            = lapply(rhs.formula.list, model.matrix, data = data)
betalist = list(
  c(-1, 1),
  c(0, -0.5, 0.5),
  c(1, -.25, .25, 0)
)
eta = mapply(function(X, beta) X %*% beta, X = Xlist, beta = betalist)
mu  = mapply(function(eta, family) family$linkinv(eta), eta = data.frame(eta), family = family.list )

data = cbind(gc_generate(mu, corr, family.list, disp = c(2, 1, 1)), data)




fitlist = mapply(
  function(formula, family) glm(formula, family, data = data), 
  formula = formula.list, family = family.list,
  SIMPLIFY = F
)

fit = copula.glm(formula.list, family.list, data, lambda = 0.01, chains = 1, iter_warmup = 1000, iter_sampling = 1000)
subset = as.mcmc( as.matrix(fit$samples)[, c(fit$coef.names, fit$var.names, fit$corr.names)] )
# 
# 
cbind(
  'post.mean' = colMeans(subset[, c(fit$coef.names, fit$var.names, fit$corr.names)]),
  'mle'       = c(unlist(sapply(fitlist, coef)), summary(fitlist[[1]])$dispersion, rep(NA, choose(J, 2))),
  'truth'     = c( unlist(betalist), 4, corr[lower.tri(corr)])
)


