remove(list = ls())
nchains = max(1, parallel::detectCores() - 1)
options(mc.cores = nchains)
## source gaussian copula R file to make functions available
source('~/Projects/Paper2/cmdstanr/glm_copula_rfuns.R')

n = 200
J = 3
p = 3
repeat { 
  corr = randcorr::randcorr(J)
  LT = abs( corr[lower.tri(corr)] )
  if ( all(LT <= 0.60) )
    break
}
data = matrix(rnorm(n*p, sd = 2), nrow = n, ncol = p)
colnames(data) = paste0('x', 1:p)
data = data.frame(data)
formula.list = list(
  y1 ~ x1,
  y2 ~ x1 + x2,
  y3 ~ x1 + x2 + x3
)
family.list = list(
  gaussian(), binomial(), binomial()  ## can also be poisson
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

data = cbind(copula.glm.generate(mu, corr, family.list, disp = c(2, 1, 1)), data)


## frequentist fit
fitlist = mapply(
  function(formula, family) glm(formula, family, data = data), 
  formula = formula.list, family = family.list,
  SIMPLIFY = F
)

## Bayes copula fit and posterior prediction
# fit      = copula.glm(formula.list, family.list, data, lambda = 0.001, chains = nchains, iter_warmup = 1000, iter_sampling = 3000)
# fit.2b1n = copula.glm.2b1n(formula.list, family.list, data, lambda = 0.001, chains = nchains, iter_warmup = 1000, iter_sampling = 3000)
# summ.fit = summary(fit$samples)
# summ.2b1n = summary(fit.2b1n$samples)

# rbind(
#   summ.fit[grepl('y3', summ.fit$variable), ],
#   summ.2b1n[grepl('y3', summ.2b1n$variable), ],
# )
# summary(fitlist[[3]])
# pred  = copula.glm.predict(fit, data, indx = 1)
# pp    = copula.glm.postprobs(fit, 'x1', nsmpl.normal = 10000, thresh = 0)  ## using normal approximation
# pp2   = copula.glm.postprobs(fit, 'x1', nsmpl.normal = 0,     thresh = 0)  ## using samples

fit = copula.glm.2b1n.postprobs(formula.list, family.list, data, lambda = 0.001, chains = nchains, iter_warmup = 1000, iter_sampling = 500, method = 'clt')

