RT.thres <-
function(R.thres) {
  is.even <- function(x) x %% 2 == 0
  cat = 1:length(R.thres)
  n.cat = length(R.thres)
  cat1 = cat[!is.even(cat)]
  cat2 = cat[is.even(cat)]
  b = cbind(R.thres[cat1], R.thres[cat2])
  b1.r = b[,1]
  b2.r = b[,2]
  b1.rt = b2.r + log((-1 + sqrt(1 + 4 * exp(b1.r-b2.r)))/2)
  b2.rt = b1.r + b2.r - b1.rt
  b12.rt = c(rbind(b1.rt, b2.rt))
  return(list(b1.rt = b1.rt, b2.rt = b2.rt, b12.rt = b12.rt))
}
