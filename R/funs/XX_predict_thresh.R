thresh.data = readRDS('~/Projects/Paper2/Results/corr_threshold.rds')

cors.data   = thresh.data$thresh[, 1:3]
thresh.data = thresh.data$thresh[, 4:7]

predict_thresh = function(corvec, nclose = 3) {
  dist = sqrt( colSums((t(cors.data) - corvec)^2) )
  if (nclose == 1) {
    return( thresh.data[which.min(dist),] )
  } else {
    indx = order(dist)[1:nclose]
  }
  wt  = 1 / dist[indx]
  val = thresh.data[indx, ]
  colSums(wt * val) / sum(wt)
}



