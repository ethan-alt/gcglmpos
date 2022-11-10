remove(list = ls())

library(tidyverse)
source('~/Projects/Paper2/cmdstanr/glm_copula_rfuns.R')

## get function to find threshold
source('~/Projects/Paper2/R/funs/XX_predict_thresh.R')

alpha = .05
gamma = 1 - alpha

id = as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if (is.na(id) )
  id = 7

corrvec   = c('mod_neg', 'low_neg', 'independent', 'low_pos', 'mod_pos', 'high_pos', 'pos')
corr      = corrvec[id]
dirs      = file.path('~/Projects/Paper2/Results', corr)
suppressMessages(
  for ( j in seq_along(dirs) ) {
    setwd(dirs[j])
    nfiles = length( list.files(pattern = '.rds') )
    print(paste('Beginning', corr[j]))
    for ( i in 1:nfiles ) {
      if ( i %% 100 == 0 )
        print( paste( round(i/nfiles * 100, 1), '% finished' ) )
      filename = paste0('postprobs_', i, '.rds')
      if ( !(file.exists(filename) ) )
        next
      file = readRDS(paste0('postprobs_', i, '.rds'))
      if ( any( is.na( file$post.probs ) ) )
        colnames(file$post.probs) = c('postprob_MSWS', 'postprob_TUG', 'postprob_BBS', 'cor_12', 'cor_13', 'cor_23', 'pval_MSWS', 'pval_TUG', 'pval_BBS')
      temp = data.frame('corr' = file$corr, 'hyp' = file$hyp, 'n' = file$n, file$post.probs)
      rownames(temp) = file$data.indx
      ## Remove sims resulting in error or incompleted sims
      indx.drop = !(complete.cases(temp)) | (rowSums(temp[, c('pval_MSWS', 'pval_TUG', 'pval_BBS')]) == 0)
      temp = temp[!indx.drop, ]
      
      ## Extract correlations and get thresholds
      corrs  = temp[, grepl('cor_', names(temp))]
      thresh = t(apply(corrs, 1, predict_thresh, nclose = 3))
      
      temp = temp %>%
        group_by(corr, hyp, n) %>%
        summarize(
          'ndatasets'           = sum(complete.cases(temp))
          , 'pos_MSWS'          = mean(postprob_MSWS > gamma)
          , 'freq_MSWS'         = mean(pval_MSWS < alpha)
          , 'pos_TUG'           = mean(postprob_TUG > gamma)
          , 'freq_TUG'          = mean(pval_TUG < alpha)
          , 'pos_BBS'           = mean(postprob_BBS > gamma)
          , 'freq_BBS'          = mean(pval_BBS < alpha)
          , 'pos_MSWS_TUG_BBS'  = mean( pmax(postprob_MSWS, postprob_BBS, postprob_TUG) > thresh[, 'thresh_123'] )
          , 'freq_MSWS_TUG_BBS' = mean( pmin( pval_MSWS, pval_BBS, pval_TUG ) < alpha / 3 )
          , 'pos_MSWS_TUG'      = mean( pmax(postprob_MSWS, postprob_TUG) > thresh[, 'thresh_12'] )
          , 'freq_MSWS_TUG'     = mean( pmin( pval_MSWS, pval_TUG ) < alpha / 2 )
          , 'pos_MSWS_BBS'      = mean( pmax(postprob_MSWS, postprob_BBS) > thresh[, 'thresh_13'] )
          , 'freq_MSWS_BBS'     = mean( pmin( pval_MSWS, pval_BBS ) < alpha / 2 )
          , 'pos_TUG_BBS'       = mean( pmax(postprob_TUG, postprob_BBS) > thresh[, 'thresh_23'] )
          , 'freq_TUG_BBS'      = mean( pmin( pval_TUG, pval_BBS ) < alpha / 2 )
          , 'pos_MSWS.TUG_BBS'  = mean( postprob_MSWS > gamma & ( pmax(postprob_TUG, postprob_BBS) > thresh[, 'thresh_23'] ) )
          , 'freq_MSWS.TUG_BBS' = mean( (pval_MSWS < alpha) & ( pmin(pval_TUG, pval_BBS) < alpha/2 ) )
        )
      if ( i == 1 ) {
        res = temp
        next
      }
      res = rbind(res, temp)
    }
    summ = res %>% 
      group_by(corr, hyp, n) %>%
      summarize(
        'ndata'                = sum(ndatasets)
        , 'bayes_MSWS'         = weighted.mean(pos_MSWS,          w = ndatasets)
        , 'freq_MSWS'          = weighted.mean(freq_MSWS,         w = ndatasets)
        , 'bayes_TUG'          = weighted.mean(pos_TUG,           w = ndatasets)
        , 'freq_TUG'           = weighted.mean(freq_TUG,          w = ndatasets)
        , 'bayes_BBS'          = weighted.mean(pos_BBS,           w = ndatasets)
        , 'freq_BBS'           = weighted.mean(freq_BBS,          w = ndatasets)
        , 'bayes_MSWS_TUG_BBS' = weighted.mean(pos_MSWS_TUG_BBS,  w = ndatasets)
        , 'freq_MSWS_TUG_BBS'  = weighted.mean(freq_MSWS_TUG_BBS, w = ndatasets)
        , 'bayes_MSWS_TUG'     = weighted.mean(pos_MSWS_TUG,      w = ndatasets)
        , 'freq_MSWS_TUG'      = weighted.mean(freq_MSWS_TUG,     w = ndatasets)
        , 'bayes_MSWS_BBS'     = weighted.mean(pos_MSWS_BBS,      w = ndatasets)
        , 'freq_MSWS_BBS'      = weighted.mean(freq_MSWS_BBS,     w = ndatasets)
        , 'bayes_TUG_BBS'      = weighted.mean(pos_TUG_BBS,       w = ndatasets)
        , 'freq_TUG_BBS'       = weighted.mean(freq_TUG_BBS,      w = ndatasets)
        , 'bayes_MSWS.TUG_BBS' = weighted.mean(pos_MSWS.TUG_BBS,  w = ndatasets)
        , 'freq_MSWS.TUG_BBS'  = weighted.mean(freq_MSWS.TUG_BBS, w = ndatasets)
      )
    saveRDS(summ, file = file.path('~/Projects/Paper2/Results', paste0('results_', corr, '.rds')))
  }
)
