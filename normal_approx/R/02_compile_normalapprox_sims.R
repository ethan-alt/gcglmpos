


dir = '~/Projects/Paper2/normal_approx'

resdir = file.path(dir, 'Results')
nfiles = length( list.files(resdir, pattern = '.rds') )


for ( i in 1:nfiles ) {
  file = readRDS( file.path(resdir, paste0('res_', i, '.rds') ) )
  if ( i == 1 ) {
    res = file
    next
  }
  res = rbind(res, file)
}


saveRDS(
  res, file = file.path(dir, 'normalapprox_compiled_sims.rds')
)


