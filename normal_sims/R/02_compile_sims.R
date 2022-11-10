remove(list = ls())
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(latex2exp)

setwd('~/Projects/Paper2/normal_sims/Results')
nfiles = length(list.files(pattern = '.rds'))


height = 9
width  = 7

for ( i in 1:nfiles ) {
  file = readRDS(paste0('res_', i, '.rds'))
  if ( i == 1 ) {
    res = file
    next
  }
  res = rbind(res, file)
}

long = res %>% pivot_longer(
  cols = !(c('hyp', 'n', 'corr')),
  names_to = 'method',
  values_to = 'pos'
)

# long$cor.latex  = paste0('$\\rho = ', long$corr, '$')
# corr.latex.lvls = paste0('$\\rho = ', sort(unique(long$corr)), '$')

long$cor.string <- formatC(long$corr, digits = 1, format = 'f')
long$cor.latex  <- factor(
  long$cor.string, levels = sort(unique(long$cor.string))
  , ordered = TRUE
)
levels(long$cor.latex) <- sapply(levels(long$cor.latex), function(x) {
  TeX( paste0('$\\rho = ', x, '$' ) )
})

method = str_split_fixed(long$method, '_', n = 2)
long$method = method[, 1]
long$endpt  = method[, 2]

old = c('bayes', 'freq', 'hmp', 'nbayes', 'nfreq', 'nhmp')
new = c('Bayesian', 'Frequentist', 'HMP', 'Bayesian', 'Frequentist', 'HMP')
names(new) = old
long$method = long$method %>% recode(!!!new)

long$endpt = ifelse(long$endpt == '', 'nrejections', long$endpt)

wide = long %>% pivot_wider(
  id_cols = c(hyp, n, corr),
  values_from = pos,
  names_from = c(method, endpt)
)


saveRDS(
  list(wide = wide, long = long),
  file = '~/Projects/Paper2/normal_sims/normal_sims_compiled.rds'
)


## Plot two-way
ptwoway.h0 <- ggplot(
  data = long %>%
    filter(endpt == '12', hyp == 'h0')
  , aes(x = n, y = pos, color = method, lty = method)
) + 
  geom_smooth(size = 0.25, se = F) +
  facet_wrap(~cor.latex, ncol = 1, labeller = label_parsed) + 
  theme(legend.title = element_blank()) + 
  scale_color_tableau() +
  ylim(0, 0.060) + 
  ylab('Family-wise error rate') + 
  xlab('Sample size') + 
  ggtitle('Type I error')

ptwoway.h1 <- ggplot(
  data = long %>%
    filter(endpt == '12', hyp == 'h1')
  , aes(x = n, y = pos, color = method, lty = method)
) + 
  geom_smooth(size = 0.25, se = F) +
  facet_wrap(~cor.latex, ncol = 1, labeller = label_parsed) + 
  theme(legend.title = element_blank()) +
  scale_color_tableau() + 
  ylab('Power') + 
  xlab('Sample size') + 
  ggtitle('Power')

ggarrange(ptwoway.h0, ptwoway.h1, common.legend = TRUE, legend = 'bottom')
ggsave('~/Projects/Paper2/normal_sims/sims_normal_2way.pdf', height = height, width = width, units = 'in', dpi = 500)







## Plot three-way
pthreeway.h0 <- ggplot(
  data = long %>%
    filter(endpt == '123', hyp == 'h0')
  , aes(x = n, y = pos, color = method, lty = method)
) + 
  geom_smooth(size = 0.25, se = F) +
  facet_wrap(~cor.latex, ncol = 1, labeller = label_parsed) + 
  theme(legend.title = element_blank()) + 
  scale_color_tableau() +
  ylim(0, 0.060) + 
  ylab('Family-wise error rate') + 
  xlab('Sample size') + 
  ggtitle('Type I error')

pthreeway.h1 <- ggplot(
  data = long %>%
    filter(endpt == '123', hyp == 'h1')
  , aes(x = n, y = pos, color = method, lty = method)
) + 
  geom_smooth(size = 0.25, se = F) +
  facet_wrap(~cor.latex, ncol = 1, labeller = label_parsed) + 
  theme(legend.title = element_blank()) +
  scale_color_tableau() + 
  ylab('Power') + 
  xlab('Sample size') + 
  ggtitle('Power')

ggarrange(pthreeway.h0, pthreeway.h1, common.legend = TRUE, legend = 'bottom')
ggsave('~/Projects/Paper2/normal_sims/sims_normal_3way.pdf', height = height, width = width, units = 'in', dpi = 500)



ggsave('~/Projects/Paper2/normal_sims/sims_normal_3way_type1error.pdf'
       , plot = pthreeway.h0
       ,height = height, width = width, units = 'in', dpi = 500)

ggsave('~/Projects/Paper2/normal_sims/sims_normal_3way_power.pdf'
       , plot = pthreeway.h1
       ,height = height, width = width, units = 'in', dpi = 500)



## Expected number of rejected hypotheses
pnrej.3way <- ggplot(
  data = long %>%
    filter(endpt == 'nrejections', hyp == 'h1')
  , aes(x = n, y = pos, color = method, lty = method)
) + 
  geom_smooth(size = 0.25, se = F) +
  facet_wrap(~cor.latex, ncol = 1, labeller = label_parsed) + 
  theme(legend.title = element_blank()) +
  scale_color_tableau() + 
  ylab('Power') + 
  xlab('Sample size') + 
  ggtitle('Power')
