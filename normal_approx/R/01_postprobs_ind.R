library(cmdstanr)
library(tidyverse)

n         = 100
warmup    = 1000
nsmpl     = 15000
ndatasets = 800
n.clt     = 2000
corr      = 'independent'


output.dir = '/pine/scr/e/t/ethanalt/Projects/Paper2/cmdstanr'


datapath = paste0("~/Projects/Paper2/vprior/power/", corr)

id = as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if(is.na(id))
  id = 530

data = readRDS(file.path(datapath, paste0('data_', id, '.rds')))
formula.list = data$formula.list
family.list  = data$family.list
data         = data$data[1:n, ]

## STAN DATA AND FUNCTION
## stan model
mod.2b1n = cmdstan_model('~/Projects/Paper2/twobernoulli_onenormal/glm_copula_opt.stan')

copula.glm.2b1n = function(formula.list, family.list, data, lambda = .01, method = 'mcmc', ...) {
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
  mod.2b1n$sample(data = standat, init = init, ...)
}

## Fit model
fit = copula.glm.2b1n(
  formula.list, family.list, data = data, 
  iter_warmup = warmup, iter_sampling = nsmpl, chains = 1, 
  output_dir = output.dir
)

## Compute posterior probabilities (MCMC)
coefnames = c('beta_b1[2]', 'beta_b2[2]', 'beta_n[2]')
postprob.mcmc        = fit$summary(coefnames, postprob = ~mean(.x > 0))[, 2]
postprob.mcmc        = as.numeric(t(postprob.mcmc))
names(postprob.mcmc) = paste0('mcmc_', 1:3)

## Compute posterior probabilities with normal approximation
clt                 = fit$summary(coefnames, postmean = ~mean(.x[1:n.clt]), postsd = ~sd(.x[1:n.clt]))
postprob.clt        = pnorm(0, mean = clt$postmean, sd = clt$postsd, lower.tail = F)
names(postprob.clt) = paste0('clt_', 1:3)

res = data.frame(t(postprob.mcmc), t(postprob.clt) )
rownames(res) = id


saveRDS(
  res,
  file = paste0('~/Projects/Paper2/normal_approx/Results/res_', id ,'.rds')
)


