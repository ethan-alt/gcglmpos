library(tidyverse)
library(tools)


setwd('~/Projects/Paper2/Results')
corrvec   = c('pos', 'mod_neg', 'low_neg', 'independent', 'low_pos', 'mod_pos', 'high_pos')


for ( i in seq_along(corrvec) ) {
  corr = corrvec[i]
  file = readRDS(paste0('results_', corr, '.rds'))
  if ( i == 1 ) {
    res = file
    next
  }
  res = rbind(res, file)
}

## Clean up names
corr.new = c('POS', 'Moderately negative', 'Low negative', 'Independent', 'Low positive', 'Moderately positive', 'High positive')
names(corr.new) = corrvec
res$corr = res$corr %>% recode(!!!corr.new)

old = c('power', 'type1error', 'pos')
new = c('BCEP', 'Type I error', 'POS')
names(new) = old
res$hyp = res$hyp %>% recode(!!!new)

long = res %>%
  pivot_longer(
    cols = -c('corr', 'hyp', 'n', 'ndata'),
    names_to = 'Method',
    values_to = 'POS'
  )


split = str_split_fixed(long$Method, pattern = '_', n = 2)
method = split[, 1]
endpt  = split[, 2]

old = c('bayes', 'freq')
new = c('Bayesian', 'Frequentist')
names(new) = old
method = method %>% recode(!!!new)

long$Method   = method
long$Endpoint = endpt





setwd('~/Projects/Paper2/Results')
saveRDS(
  list('wide' = res, 'long' = long),
  file = 'sim_results.rds'
)

