wle.a.fit1 <- function(x, .data, .beta, .w) {
  XX = .data
  wt=.w
  .rawscores = rv = rowSums(.data)
  k = ncol(.data)
  ml.select <- which(!((rv == 0) | (rv == k)))
  .data <- XX[ml.select,]
  .rawscores <- rv[ml.select]
  .w = wt[ml.select]
  teil1 <- .rawscores*x[.rawscores]
  teil2a <- matrix(rep(NA, length(.rawscores)*k), ncol=k)
  for(j in 1:k) 
    teil2a[,j] <- x[.rawscores]-.beta[j] 
  teil2 <- rowSums(log(1+exp(teil2a)))
  teil3a <- .data
  for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
  teil3 <- rowSums(teil3a)
  LL <- sum(.w * (teil1 - teil2 - teil3))  
  return(-LL)
}

# Likelihood: For extreme raw-scores: http://www.rasch.org/rmt/rmt122h.htm
wle.a.fit.extr1 <- function(x, .data, .beta, .w, d = NULL) {
  l.d = length(d)
  extr=1:l.d
  i = extr
  XX = .data
  .rawscores = rv = rowSums(.data)
  k = ncol(.data)
  sel.zeros = which(.rawscores == 0)
  sel.ks = which(.rawscores == k)
  if(l.d == 2){
    .rawscores.zeros = .rawscores[sel.zeros] = 1 - d[1]
    .rawscores.ks = .rawscores[sel.ks] = d[2]
  } else if(l.d == 3){
    if(d[2] < 1){
      .rawscores.zeros1 = .rawscores[sel.zeros] = 1 - d[1]
      .rawscores.zeros2 = .rawscores[sel.zeros] = 1 - d[2]
      .rawscores.ks = .rawscores[sel.ks] = d[3]   
    } else {
      .rawscores.zeros = .rawscores[sel.zeros] = 1 - d[1]
      .rawscores.ks1 = .rawscores[sel.ks] = d[2]
      .rawscores.ks2 = .rawscores[sel.ks] = d[3]  
    }
  } 
  else if(l.d == 4){
    .rawscores.zeros1 = .rawscores[sel.zeros] = 1 - d[1]
    .rawscores.zeros2 = .rawscores[sel.zeros] = 1 - d[2]
    .rawscores.ks1 = .rawscores[sel.ks] = d[3]
    .rawscores.ks2 = .rawscores[sel.ks] = d[4]
  } 
  ### For the zeros
  if(l.d == 2){
    if(i == 1){
      teil1.zeros <- .rawscores.zeros*x
      teil2a.zeros <- matrix(rep(NA, length(.rawscores.zeros)*k), ncol=k)
      for(j in 1:k) 
        teil2a.zeros[,j] <- x-.beta[j] 
      teil2.zeros <- rowSums(log(1+exp(teil2a.zeros)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.zeros <- sum(.w * (teil1.zeros - teil2.zeros - teil3))
      LL = LL.zeros
    } else  {
      ### For the ks
      teil1.ks <- .rawscores.ks * x
      teil2a.ks <- matrix(rep(NA, length(.rawscores.ks)*k), ncol=k)
      for(j in 1:k) 
        teil2a.ks[,j] <- (x -.beta[j]) 
      teil2.ks <- rowSums(log(1+exp(teil2a.ks)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.ks <- sum(.w * (teil1.ks - teil2.ks - teil3))
      LL = LL.ks
    } 
  } else if(l.d == 3) {
    if(d[2] < 1){
      if(i == 1){
        teil1.zeros <- .rawscores.zeros1*x
        teil2a.zeros <- matrix(rep(NA, length(.rawscores.zeros1)*k), ncol=k)
        for(j in 1:k) 
          teil2a.zeros[,j] <- x-.beta[j] 
        teil2.zeros <- rowSums(log(1+exp(teil2a.zeros)))
        teil3a <- .data
        for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
        teil3 <- rowSums(teil3a, na.rm = T)
        LL.zeros <- sum(.w * (teil1.zeros - teil2.zeros - teil3))
        LL = LL.zeros
      } else if(i == 2) {
        teil1.zeros <- .rawscores.zeros2*x
        teil2a.zeros <- matrix(rep(NA, length(.rawscores.zeros2)*k), ncol=k)
        for(j in 1:k) 
          teil2a.zeros[,j] <- x-.beta[j] 
        teil2.zeros <- rowSums(log(1+exp(teil2a.zeros)))
        teil3a <- .data
        for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
        teil3 <- rowSums(teil3a, na.rm = T)
        LL.zeros <- sum(.w * (teil1.zeros - teil2.zeros - teil3))
        LL = LL.zeros
      } else if (i == 3) {
        teil1.ks <- .rawscores.ks * x
        teil2a.ks <- matrix(rep(NA, length(.rawscores.ks)*k), ncol=k)
        for(j in 1:k) 
          teil2a.ks[,j] <- (x -.beta[j]) 
        teil2.ks <- rowSums(log(1+exp(teil2a.ks)))
        teil3a <- .data
        for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
        teil3 <- rowSums(teil3a, na.rm = T)
        LL.ks <- sum(.w * (teil1.ks - teil2.ks - teil3))
        LL = LL.ks
      }
    } else if(i == 1){
      teil1.zeros <- .rawscores.zeros*x
      teil2a.zeros <- matrix(rep(NA, length(.rawscores.zeros)*k), ncol=k)
      for(j in 1:k) 
        teil2a.zeros[,j] <- x-.beta[j] 
      teil2.zeros <- rowSums(log(1+exp(teil2a.zeros)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.zeros <- sum(.w * (teil1.zeros - teil2.zeros - teil3))
      LL = LL.zeros
    } else if(i == 2) {
      teil1.ks <- .rawscores.ks1 * x
      teil2a.ks <- matrix(rep(NA, length(.rawscores.ks1)*k), ncol=k)
      for(j in 1:k) 
        teil2a.ks[,j] <- (x -.beta[j]) 
      teil2.ks <- rowSums(log(1+exp(teil2a.ks)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.ks <- sum(.w * (teil1.ks - teil2.ks - teil3))
      LL = LL.ks
    } else if (i == 3) {
      teil1.ks <- .rawscores.ks2 * x
      teil2a.ks <- matrix(rep(NA, length(.rawscores.ks2)*k), ncol=k)
      for(j in 1:k) 
        teil2a.ks[,j] <- (x -.beta[j]) 
      teil2.ks <- rowSums(log(1+exp(teil2a.ks)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.ks <- sum(.w * (teil1.ks - teil2.ks - teil3))
      LL = LL.ks
    }
  } else if(l.d == 4) {
    if(i == 1){
      teil1.zeros <- .rawscores.zeros1*x
      teil2a.zeros <- matrix(rep(NA, length(.rawscores.zeros1)*k), ncol=k)
      for(j in 1:k) 
        teil2a.zeros[,j] <- x-.beta[j] 
      teil2.zeros <- rowSums(log(1+exp(teil2a.zeros)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.zeros <- sum(.w * (teil1.zeros - teil2.zeros - teil3))
      LL = LL.zeros
    } else if(i == 2) {
      teil1.zeros <- .rawscores.zeros2*x
      teil2a.zeros <- matrix(rep(NA, length(.rawscores.zeros2)*k), ncol=k)
      for(j in 1:k) 
        teil2a.zeros[,j] <- x-.beta[j] 
      teil2.zeros <- rowSums(log(1+exp(teil2a.zeros)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.zeros <- sum(.w * (teil1.zeros - teil2.zeros - teil3))
      LL = LL.zeros
    } else if (i == 3) {
      teil1.ks <- .rawscores.ks1 * x
      teil2a.ks <- matrix(rep(NA, length(.rawscores.ks1)*k), ncol=k)
      for(j in 1:k) 
        teil2a.ks[,j] <- (x -.beta[j]) 
      teil2.ks <- rowSums(log(1+exp(teil2a.ks)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.ks <- sum(.w * (teil1.ks - teil2.ks - teil3))
      LL = LL.ks
    } else if (i == 4) {
      teil1.ks <- .rawscores.ks2 * x
      teil2a.ks <- matrix(rep(NA, length(.rawscores.ks2)*k), ncol=k)
      for(j in 1:k) 
        teil2a.ks[,j] <- (x -.beta[j]) 
      teil2.ks <- rowSums(log(1+exp(teil2a.ks)))
      teil3a <- .data
      for(j in 1:k)  teil3a[,j] <- apply(.data,2,sum)[j]*.beta[j]
      teil3 <- rowSums(teil3a, na.rm = T)
      LL.ks <- sum(.w * (teil1.ks - teil2.ks - teil3))
      LL = LL.ks
    }
  }
  return(-LL)
}
# Probability for the SE estimation of the person parameters
P.i <- function(x, .beta) {
#   .beta=b.w
  k=length(.beta)
  teil2a <- matrix(NA, k-1, k)
  for (j in 1:k)  teil2a[,j] <- x - .beta[j] 
  teil2 <- (1+exp(teil2a))
  teil3 <- exp(teil2a)
  LL <- teil3/teil2  
  return(LL)
}