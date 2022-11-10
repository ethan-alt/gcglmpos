remove(list = ls())

library(truncnorm)
library(mvtnorm)
library(tidyverse)
library(harmonicmeanp)

corrvec = c(-0.4, -0.2, 0.0, 0.2, 0.4, 0.8)
nvec    = seq(10, 1600, by = 20)
hyp     = c('h0', 'h1')
nsims   = 100000
nsims = 10000
## get function to find threshold
source('~/Projects/Paper2/R/funs/XX_predict_thresh.R')

grid = expand.grid('corr' = corrvec, 'hyp' = hyp, 'n' = nvec, stringsAsFactors = F)
grid = expand.grid(
  'n' = nvec, 'hyp' = hyp, 'corr' = corrvec,
  stringsAsFactors = F
)

corvec2mat = function(x, J = 3) {
  matrix(x, J, J) + diag(1 - x, J)
}

id = as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if(is.na(id))
  id = 125

rho     = grid[id, 'corr']
hyp     = grid[id, 'hyp']
n       = grid[id, 'n']

cor = corvec2mat(rho)


postprobs = matrix(nrow = nsims, ncol = 3)
colnames(postprobs) = c(
  paste0('postProb_', 1:3)
)

# thresh = matrix(nrow = nsims, ncol = 4)
# colnames(thresh) = paste0('thresh_', c('123', '12', '13', '23'))
postprobs.mat = matrix(nrow = nsims, ncol = 3)
thresh.mat    = matrix(nrow = nsims, ncol = 4)
ind.all       = matrix(nrow = nsims, ncol = 2)
for ( i in 1:nsims ) {
  if ( hyp == 'h0' ) {
    mu = numeric(3)
  } else {
    mu = rtruncnorm(n = 3, a = 0, mean = 0, sd = 0.1)
  }
  Y = rmvnorm(n, mean = mu, sigma = cor)
  postMean  = colMeans(Y)
  postSigma = cov(Y) / n
  postSD    = sqrt(diag(postSigma))
  
  ## Compute posterior probability; order them
  postProbs      = pt(postMean / postSD, df = n - 1)
  ordIndx        = order(postProbs, decreasing = T)
  postProbs.ord  = postProbs[ordIndx]
  postprobs.mat[i, ] = postProbs
  
  ## Obtian posterior correlation--and get threshold
  postCor         = cov2cor(postSigma)
  postCor         = postCor[lower.tri(postCor)]
  thresh          = predict_thresh(postCor)
  thresh.mat[i, ] = thresh
  thresh          = c(thresh[1], thresh[ paste0('thresh_', paste0( sort(ordIndx[-1]), collapse = '') ) ] )
  ind.all[i, ]    = c(
    (postProbs.ord[1] > thresh[1]) & (postProbs.ord[2] > thresh[2]) & (postProbs.ord[3] > 0.95),
    (postProbs.ord[1] > (1 - .05/3)) & (postProbs.ord[2] > (1 - .05/2)) & (postProbs.ord[3] > 0.95)
  )
}

colnames(postprobs.mat) = paste0('postProb_', 1:3)
colnames(thresh.mat)    = paste0('thresh_', c('123', '12', '13', '23'))
colnames(ind.all)       = c('bayes_all', 'freq_all')

df = data.frame(postprobs.mat, thresh.mat)

ind.succ = df %>%
  summarize(
      'bayes_1'   = postProb_1 > 0.95
    , 'bayes_2'   = postProb_2 > 0.95
    , 'bayes_3'   = postProb_3 > 0.95
    , 'bayes_123' = pmax(postProb_1, postProb_2, postProb_3) > thresh_123
    , 'freq_123'  = pmax(postProb_1, postProb_2, postProb_3) > (1 - 0.05 / 3)
    , 'bayes_12'  = pmax(postProb_1, postProb_2) > thresh_12
    , 'freq_12'   = pmax(postProb_1, postProb_2) > (1 - .05 / 2)
    , 'bayes_13'  = pmax(postProb_1, postProb_3) > thresh_13
    , 'freq_13'   = pmax(postProb_1, postProb_3) > (1 - .05 / 2)
    , 'bayes_23'  = pmax(postProb_2, postProb_3) > thresh_23
    , 'freq_23'   = pmax(postProb_2, postProb_3) > (1 - .05 / 2)
  )
lst <- list(
  c('postProb_1', 'postProb_2')
  , c('postProb_1', 'postProb_3')
  , c('postProb_2', 'postProb_3')
  , c('postProb_1', 'postProb_2', 'postProb_3')
)


## Do averaged p-value approach
df1m           <- df %>%  select(paste0('postProb_', 1:3)) %>% mutate_all(function(x) 1 - x)
df1m <- lapply(df1m, function(x) ifelse(x == 0, 1e-20, x)) %>% data.frame
hmp            <- sapply(1:length(lst), function(i) { p <- df1m %>% select(lst[[i]]); apply(p, 1, p.hmp, L = ncol(p) ) } )
hmp.ind        <- data.frame(hmp < 0.05)
names(hmp.ind) <- c('hmp_12', 'hmp_13', 'hmp_23', 'hmp_123')
ind.succ       <- cbind(ind.succ, data.frame(hmp.ind))

ind.succ$nbayes = ind.succ %>%
  select('bayes_123', 'bayes_12', 'bayes_13', 'bayes_23') %>%
  rowSums

ind.succ$nfreq  = ind.succ %>%
  select('freq_123', 'freq_12', 'freq_13', 'freq_23') %>%
  rowSums

ind.succ$nhmp = ind.succ %>%
  select('hmp_12', 'hmp_13', 'hmp_23', 'hmp_123') %>%
  rowSums

res = data.frame(hyp = hyp, n = n, corr = rho, t(colMeans(ind.succ)))

saveRDS(
  res, file = file.path( '~/Projects/Paper2/normal_sims/Results/', paste0('res_', id, '.rds') )
)


