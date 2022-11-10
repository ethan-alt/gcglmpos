library(cmdstanr)
library(posterior)
## c('pos', 'mod_neg', 'low_neg', 'independent', 'low_pos', 'mod_pos', 'high_pos')
# corr = 'pos'
# corr = 'mod_neg'
# corr = 'low_neg'
# corr = 'independent'
# corr = 'low_pos'
# corr = 'med_pos'
corr = 'high_pos'


## sim / cluster arguments
nvec          = seq(300, 500, by = 25)
lambda        = c(0.5, 0.75, 1.0)
ndatasets     = 10000
each.cl       = 250  ## number of data sets per array
output.dir    = '/pine/scr/e/t/ethanalt/Projects/Paper2/cmdstanr'
if(!dir.exists(output.dir))
  dir.create(output.dir, recursive = T)

## MCMC sampling arguments
chains   = 1
warmup   = 1000
nsmpl    = 2000
clt.smpl = 10000

# ## FOR TESTING--COMMENT OUT
# ndatasets  = 2
# each.cl    = 2
# warmup = 500
# nsmpl = 1000

save.every = min(5, each.cl)










hyp = c('type1error', 'power')
if (corr == 'pos')
  hyp = 'pos'

grid       = expand.grid('n' = nvec, 'hyp' = hyp, 'corr' = corr, stringsAsFactors = F)
grid       = grid[rep(seq_len(nrow(grid)), each = ndatasets / each.cl), ]
grid$end   = seq(from = each.cl, to = ndatasets, by = each.cl)
grid$start = grid$end - (each.cl - 1)


id = as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if(is.na(id)) {
  output.dir = NULL
  id = nrow(grid)
}


n.id        = grid[id, 'n']
hyp.id      = grid[id, 'hyp']
corr.id     = grid[id, 'corr']
start.id    = grid[id, 'start']
end.id      = grid[id, 'end']
indx.id     = start.id:end.id

folder = ifelse(corr.id == 'pos', 'pos', file.path(corr.id,hyp.id))
if ( corr.id != 'pos' )
  folder = file.path(corr.id,hyp.id)

savepath = file.path('~/Projects/Paper2/Results', folder)
if(!dir.exists(savepath))
  dir.create(savepath, recursive = T)

folder = ifelse(corr.id == 'pos', 'pos', file.path(corr,hyp.id))
folder.vprior = ifelse(corr.id == 'pos', 'pos', file.path(hyp.id,corr))
if ( corr.id != 'pos' )
  folder = file.path(corr,hyp.id)

## Create save path and get path of data from vprior
savepath = file.path('~/Projects/Paper2/Results', folder)
datapath = file.path('~/Projects/Paper2/vprior', folder.vprior)

## get formula / family
temp = readRDS(file.path(datapath, paste0('data_', 1, '.rds')))
formula.list = temp$formula.list
family.list  = temp$family.list
remove('temp')





## STAN DATA AND FUNCTION
## stan model
mod.2b1n = cmdstan_model('~/Projects/Paper2/twobernoulli_onenormal/glm_copula_opt.stan')

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





## loop through data sets and compute posterior probability / correlation
res = matrix(nrow = each.cl, ncol = 6)
rownames(res) = indx.id
for ( i in 1:each.cl ) {
  data = readRDS(file.path(datapath, paste0('data_', i, '.rds')))$data
  res[i, ] = tryCatch(
    copula.glm.2b1n.postprobs(
      formula.list, family.list, data, method = 'clt', clt.nsmpl = clt.smpl, 
      chains = chains, iter_warmup = warmup, iter_sampling = nsmpl, lambda = .01,
      output_dir = output.dir, output_basename = paste0("smpl_", id, '_', indx.id[i] )
      ),
    error = function(e) rep(NA, times = 6)
  )
  if ( i %% save.every == 0)
    saveRDS(
      list(
        'post.probs' = res,
        'n'          = n.id,
        'data.indx'  = indx.id,
        'corr'       = corr.id,
        'hyp'        = hyp.id
      ), file = file.path(savepath, paste0('postprobs_', id, '.rds'))
    )
}

saveRDS(
  list(
    'post.probs' = res,
    'n'          = n.id,
    'data.indx'  = indx.id,
    'corr'       = corr.id,
    'hyp'        = hyp.id
  ), file = file.path(savepath, paste0('postprobs_', id, '.rds'))
)