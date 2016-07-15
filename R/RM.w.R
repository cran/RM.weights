RM.w <-
function(.data, .w = NULL, .d=NULL, country=NULL, se.control = T,
                quantile.seq = NULL, write.file = F) {
  if(is.null(country)) country = "country"
  if(is.null(.w)) .w = rep(1, nrow(.data))
	rv = rowSums(.data)
	k = ncol(.data)
	XX = .data
	wt = .w
  if(is.null(.d)) .d = c(0.5,(k-1)+0.5)
  d = .d
  l.d = length(d)
	wle.fit <- function(x, .data, .rawscores, .w) {
		.data <- XX
		.rawscores = rv
		ml.select <- which(!((.rawscores == 0) | (.rawscores == k)))
		.data = .data[ml.select,]
		.rawscores = .rawscores[ml.select]
		.w = wt[ml.select]
    gamma.r.v = elementary_symmetric_functions(x)$"0"
		teil1a <- - as.matrix(.data) %*% x
		LL <- sum(.w * (teil1a - log(gamma.r.v[.rawscores+1])))
		return(-LL)
	}
	P.i.b <- function(x) {
		.data <- XX
		.rawscores = rv
		ml.select <- which(!((.rawscores == 0) | (.rawscores == k)))
		.data = .data[ml.select, ]
		.rawscores = .rawscores[ml.select]
		.w = wt[ml.select]
    .w = (.w) * length(ml.select)/sum(.w)
		gamma.r.v = elementary_symmetric_functions(x)$"0"
		p = LL = matrix(NA, ncol=k, nrow=nrow(.data))
		for(i in 1:k){ 
		gamma.r.v1 = elementary_symmetric_functions(x[(-i)],order=0)$"0"
		teil1a <- exp(- x[i]) * gamma.r.v1[.rawscores]
		p[,i] = teil1a / gamma.r.v[.rawscores+1] 
		LL[,i] <- (p[,i] * (1 - p[,i]))*.w
		} 
		LL = rowsum(LL, group = .rawscores)	
		return(LL)
	}
  # Estimation: item
	opt.w <- optim(seq(-3,3,length.out=k), wle.fit, method = "BFGS", hessian = T)
  b.w = opt.w$par
  names(b.w) = colnames(XX)
	p.i.b = colSums( P.i.b(b.w) )
	se.b.w = sqrt(1/p.i.b)
  # Likelihood for the post-hoc estimation of the person parameters
	wle.a.fit <- function(x, .data, .beta, .rawscores, .w) {
    .data = XX
    k = ncol(.data)
		ml.select <- which(!((rv == 0) | (rv == k)))
		.data <- XX[ml.select,]
		.beta <- b.w
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
  wle.a.fit.extr <- function(x, .data, .beta, .rawscores, .w, d = NULL, 
                             extr=1:l.d) {
    i = extr
    if(is.null(d)) d = .d
    .data = XX
    .w = wt
    .rawscores = rv
    .beta <- b.w
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
    k=length(.beta)
		teil2a <- matrix(NA, k-1, k)
		for (j in 1:k)  teil2a[,j] <- x - .beta[j] 
		teil2 <- (1+exp(teil2a))
		teil3 <- exp(teil2a)
		LL <- teil3/teil2  
		return(LL)
	}
  # Estimation: person parameters
	opt.a = optim(seq(-3,3,length.out=k-1),wle.a.fit,method="BFGS",hessian=T)
  # Extremes
  if(l.d == 2){
    opt.a.zeros = optim(c(-3), function(i) wle.a.fit.extr(i, extr=1), 
                        method="BFGS")
    opt.a.ks = optim(c(3), function(i) wle.a.fit.extr(i, extr = 2), 
                     method="BFGS")
  } else if(l.d == 3){
    if(d[2]<1){
      opt.a.zeros1 = optim(c(-3), function(i) wle.a.fit.extr(i, extr=1), 
                          method="BFGS")
      opt.a.zeros2 = optim(c(-3), function(i) wle.a.fit.extr(i, extr=2), 
                           method="BFGS")
      opt.a.ks = optim(c(3), function(i) wle.a.fit.extr(i, extr = 3), 
                       method="BFGS")
    } else {
      opt.a.zeros = optim(c(-3), function(i) wle.a.fit.extr(i, extr=1), 
                           method="BFGS")
      opt.a.ks1 = optim(c(3), function(i) wle.a.fit.extr(i, extr = 2), 
                       method="BFGS")
      opt.a.ks2 = optim(c(3), function(i) wle.a.fit.extr(i, extr = 3), 
                        method="BFGS")
    }

  } else if(l.d == 4) {
    opt.a.zeros1 = optim(c(-3), function(i) wle.a.fit.extr(i, extr=1), 
                        method="BFGS")
    opt.a.zeros2 = optim(c(-3), function(i) wle.a.fit.extr(i, extr=2), 
                        method="BFGS")
    opt.a.ks1 = optim(c(3), function(i) wle.a.fit.extr(i, extr = 3), 
                     method="BFGS")
    opt.a.ks2 = optim(c(3), function(i) wle.a.fit.extr(i, extr = 4), 
                     method="BFGS")
  }
  # 0.5 distance - to get extremes standard error
  opt.a.mid = optim(c(-3), function(i) wle.a.fit.extr(i, d = c(0.5,0.5), extr = 1), 
                      method="BFGS")
	a.w <- opt.a$par
	se.a.w = sqrt(1/rowSums(P.i(a.w,b.w)*(1-P.i(a.w,b.w))))
  hessian.a = P.i(a.w,b.w)*(1-P.i(a.w,b.w))
  if(l.d == 2) {
    a.zero = opt.a.zeros$par
    a.k = opt.a.ks$par
  } else if(l.d == 3){
    if(d[2]<1){
      a.zero1 = opt.a.zeros1$par
      a.zero2 = opt.a.zeros2$par
      a.k = opt.a.ks$par
    } else {
      a.zero = opt.a.zeros$par
      a.k1 = opt.a.ks1$par
      a.k2 = opt.a.ks2$par
    }
  } else {
    a.zero1 = opt.a.zeros1$par
    a.zero2 = opt.a.zeros2$par
    a.k1 = opt.a.ks1$par
    a.k2 = opt.a.ks2$par
  }
  # SE Extremes
  if(l.d == 2){
    if(se.control){
      a.mid = opt.a.mid$par
      se.a.zero = se.a.k = 
        sqrt(1/rowSums(P.i(a.mid,b.w)*(1-P.i(a.mid,b.w))))[1]
    } else {
      se.a.zero = sqrt(1/rowSums(P.i(a.zero,b.w)*(1-P.i(a.zero,b.w))))[1]
      se.a.k = sqrt(1/rowSums(P.i(a.k,b.w)*(1-P.i(a.k,b.w))))[1]
    }  
  }
  if(l.d == 3){
    if(d[2]<1){
      se.a.zero1 = sqrt(1/rowSums(P.i(a.zero1,b.w)*(1-P.i(a.zero1,b.w))))[1]
      se.a.zero2 = sqrt(1/rowSums(P.i(a.zero2,b.w)*(1-P.i(a.zero2,b.w))))[1]
      se.a.k = sqrt(1/rowSums(P.i(a.k,b.w)*(1-P.i(a.k,b.w))))[1]
    } else{
      se.a.zero = sqrt(1/rowSums(P.i(a.zero,b.w)*(1-P.i(a.zero,b.w))))[1]
      se.a.k1 = sqrt(1/rowSums(P.i(a.k1,b.w)*(1-P.i(a.k1,b.w))))[1]
      se.a.k2 = sqrt(1/rowSums(P.i(a.k2,b.w)*(1-P.i(a.k2,b.w))))[1]
    }
  }
  if(l.d == 4){
    se.a.zero1 = sqrt(1/rowSums(P.i(a.zero1,b.w)*(1-P.i(a.zero1,b.w))))[1]
    se.a.zero2 = sqrt(1/rowSums(P.i(a.zero2,b.w)*(1-P.i(a.zero2,b.w))))[1]
    se.a.k1 = sqrt(1/rowSums(P.i(a.k1,b.w)*(1-P.i(a.k1,b.w))))[1]
    se.a.k2 = sqrt(1/rowSums(P.i(a.k2,b.w)*(1-P.i(a.k2,b.w))))[1]
  }
  if(se.control){
    a.mid = opt.a.mid$par
    if(l.d == 3){
      if(d[2]<1){
        se.a.zero1 = se.a.zero2 = se.a.k = 
          sqrt(1/rowSums(P.i(a.mid,b.w)*(1-P.i(a.mid,b.w))))[1]
      } else         
        se.a.zero = se.a.k1 = se.a.k2 =
        sqrt(1/rowSums(P.i(a.mid,b.w)*(1-P.i(a.mid,b.w))))[1]
    } else if(l.d == 4){
      se.a.zero1 = se.a.zero2 = se.a.k1 = se.a.k2 =
        sqrt(1/rowSums(P.i(a.mid,b.w)*(1-P.i(a.mid,b.w))))[1]
    }
  }
  # SE global
  if(l.d == 2){
    a.w = c(a.zero, a.w, a.k)
    se.a.w = c(se.a.zero, se.a.w, se.a.k) 
  }
  if(l.d == 3){
    if(d[2] < 1){
      a.w = c(a.zero1, a.zero2, a.w, a.k)
      se.a.w = c(se.a.zero1, se.a.zero2, se.a.w, se.a.k) 
    } else {
      a.w = c(a.zero, a.w, a.k1, a.k2)
      se.a.w = c(se.a.zero, se.a.w, se.a.k1, se.a.k2) 
    }
  }
  if(l.d == 4){
    a.w = c(a.zero1,a.zero2, a.w, a.k1, a.k2)
    se.a.w = c(se.a.zero1, se.a.zero2, se.a.w, se.a.k1, se.a.k2) 
  }
  # Infit and outfit statistics
	ml.select <- which(!((rv == 0) | (rv == k) |  is.na(rv)))
	rv2 = rv[ml.select]
	XX2 = XX[ml.select,]
	wt2 = wt[ml.select]
  wt2 = (wt2) * nrow(XX2)/sum(wt2)
	P.E <- function(x, .data = XX) {
	  .rawscores = rowSums(.data)
		ml.select <- which(!((.rawscores == 0) | (.rawscores == k)))
		.data = .data[ml.select, ]
		.rawscores = .rawscores[ml.select]
		.w = wt[ml.select]
		gamma.r.v = elementary_symmetric_functions(x)$"0"
		p = LL = matrix(NA, ncol=k, nrow=nrow(.data))
		for(i in 1:k){ 
		gamma.r.v1 = elementary_symmetric_functions(x[(-i)],order=0)$"0"
		teil1a <- exp(- x[i]) * gamma.r.v1[.rawscores]
		p[,i] = teil1a / gamma.r.v[.rawscores+1] 
		LL[,i] <- p[,i]
		} 
		return(LL)
	}
	E = P.E(b.w) 
	V = P.E(b.w) * (1-P.E(b.w))
	z = (XX2 - E )/sqrt(V) 
	outfit.w = apply(z^2 * wt2, 2, sum)/(sum(wt2))
	infit.w = apply(V * z^2* wt2,2,sum)/ apply(V* wt2,2,sum)
  # Residuals
    P.i2 <- function(x, .beta = 1){
      x = x[rv2]
      teil2a = XX2
      for(i in 1:ncol(XX2)){
        teil2a[,i] <- x - .beta[i] 
      }
      teil2 <- (1+exp(teil2a))
      teil3 <- exp(teil2a)
      LL <- teil3/teil2  
      return(LL)
    }
    E.i2 = P.i2(a.w, b.w)
    mat.res = ((XX2 - E.i2 ))
  # Reliability
  wt.rv2 = sapply(1:(k-1), function(i) sum(wt2[rv2==i], na.rm = T))
  wt.rv2.tot = sapply(1:(k-1), function(i) wt2[rv2==i])
  obs.mean = sum(a.w[2:k]*wt.rv2)/sum(wt.rv2) #observed distribution mean
  model.var=(a.w[2:k]-obs.mean)^2*wt.rv2
  model.var.tot =sum(model.var)
  resid.var=se.a.w[2:k]^2*wt.rv2
  resid.var.tot = sum(resid.var)
  reliab =model.var.tot/(model.var.tot+resid.var.tot)
  mean.flat=mean(a.w[2:k]) #flat distribution mean
  model.var.fl=(a.w[2:k]-mean.flat)^2
  model.var.fl.tot=sum(model.var.fl)
  resid.var.fl=se.a.w[2:k]^2
  resid.var.fl.tot = sum(resid.var.fl)
  reliab.fl =model.var.fl.tot/(model.var.fl.tot+resid.var.fl.tot)
  # Conditional correlation
  nitem = k
  ri_in = XX2
  r_wt = wt2
  nraw=k-1
  npat = 2^k-2 
  totwt = sum(r_wt)
  i_pyeswt=rep(0,nitem) #pos resp proportion
  r_nyes=rowSums(ri_in)
  rri_nyeswt= mat.or.vec(nraw,nitem) #wt sum of positive responses
  i_yeswt=rep(0,nitem) #wt sum of pos responses  
  for (rr in 1:nraw) {totrx=ri_in[,1:nitem]*(r_nyes==rr)*r_wt
                      for (c in 1:nitem) {rri_nyeswt[rr,c]=sum(totrx[,c])}}
  for(i in 1:nitem) {i_yeswt[i]=sum(rri_nyeswt[,i])} 
  i_pyeswt=i_yeswt/totwt
  rr_totwt=rep(0,nraw) #wt sum of cases
  for (rr in 1:nraw) {rr_totwt[rr]=sum(r_wt*(r_nyes==rr))}
  ii_predcorr=mat.or.vec(nitem,nitem) #predicted correlation
  ii_predyesprop=mat.or.vec(nitem,nitem) #predicted wt prop resp in 1,1
  ii_predyes=mat.or.vec(nitem,nitem) #predicted wt resp in 1,1 this case
  ii_predyestot=mat.or.vec(nitem,nitem) #predicted total wt resp in 1,1
  ii_residcorr=mat.or.vec(nitem,nitem) #residual correlation
  colnames(ii_residcorr) = rownames(ii_residcorr) = colnames(XX)
  pi_logpresp=mat.or.vec(npat,nitem) #log probability of item response if theta=0
  pi_resp=mat.or.vec(npat,nitem) #item responses
  for (p in 1:npat) {trem=p;
                     for(i in 1:nitem) {if(trem>=2^(nitem-i)) {pi_resp[p,i]=1; 
                                                               trem=trem-2^(nitem-i)} 
                                        else pi_resp[p,i]=0}}
  p_raw=rowSums(pi_resp)
  p_pinrr=rep(0,npat) #probability of pattern given raw score
  p_presp=rep(0,npat) #probability of response pattern
  rr_totp=rep(0,nraw) #sum of probabilities of patterns with raw score rr
  i_pyeswt=i_yeswt/totwt
  i_parm=rep(0,nitem) #item parameters
  i_parm2=rep(0,nitem) #adjusted item parameters
  i_pyes=1/(1+1/exp(0-b.w)) #probability of affirmative response given theta=0
  for(p in 1:npat) {for (i in 1:nitem) {
    if (pi_resp[p,i]==1) pi_logpresp[p,i]=log(i_pyes[i]) else pi_logpresp[p,i]=log(1-i_pyes[i])}}
  for (rr in 1:nraw) {rr_totp[rr]=sum(p_presp*(p_raw==rr))}
  for(p in 1:npat) p_pinrr[p]=p_presp[p]/rr_totp[p_raw[p]]
  p_presp=exp(rowSums(pi_logpresp))
  for (rr in 1:nraw) {rr_totp[rr]=sum(p_presp*(p_raw==rr))}
  for(p in 1:npat) p_pinrr[p]=p_presp[p]/rr_totp[p_raw[p]] 
  for (p in 1:npat){for (i in 1:nitem){for (ii in 1:nitem){
    ii_predyes[i,ii]=pi_resp[p,i]*pi_resp[p,ii]*p_pinrr[p]*rr_totwt[p_raw[p]]}};
                    ii_predyestot=ii_predyestot+ii_predyes}
  ii_predyesprop=ii_predyestot/totwt
  rri_outobs=rri_nyeswt/rr_totwt # Observed response proportion
  obs.prop = rri_outobs
#   lev.rv2 = as.numeric(levels(as.factor(rv2)))
#   pred.count.tot = sapply(lev.rv2, function(i) (P.E(b.w)[rv2==i,]*wt.rv2.tot[[i]]))
#   pred.count = sapply(lev.rv2, function(i) apply(pred.count.tot[[i]], 2, sum))
#   pred.prop = t(pred.count/wt.rv2)
#   colnames(pred.prop) = colnames(obs.prop) = colnames(XX)
  wt.item.tot = sapply(1:k, function(i) sum(wt2[XX2[,i] == 1]))/sum(wt2)
  cov1.obs.wt = sapply(1:k, function(j)
    sapply(1:k, function(i) sum(wt2[XX2[,j] == 1 & XX2[,i] == 1])/sum(wt2))) 
  cov2.obs.wt = sapply(1:k, function(j)
    sapply(1:k, function(i) sum(wt2[XX2[,j] == 1])/sum(wt2) * 
             sum(wt2[XX2[,i] == 1])/sum(wt2))) 
  den.obs.wt = sapply(1:k, function(j)
    sapply(1:k, function(i) 
      prod(c(sum(wt2[XX2[,j] == 1])/sum(wt2),sum(wt2[XX2[,j] == 0])/sum(wt2),
             sum(wt2[XX2[,i] == 1])/sum(wt2), sum(wt2[XX2[,i] == 0])/sum(wt2)))
    ))
  cov.obs.wt = cov1.obs.wt-cov2.obs.wt
  cor.obs.wt = cov.obs.wt/sqrt(den.obs.wt)
  ii_obscorr = cor.obs.wt
  for (i in 1:nitem){for (ii in 1:nitem){
    ii_predcorr[i,ii]=(ii_predyesprop[i,ii]-i_pyeswt[i]*i_pyeswt[ii])/
      sqrt(i_pyeswt[i]*(1-i_pyeswt[i])*i_pyeswt[ii]*(1-i_pyeswt[ii]));
    if (i==ii) ii_predcorr[i,ii]=1}}
  
  for (i in 1:nitem){for (ii in 1:nitem){
    if (i==ii) ii_residcorr[i,ii]=1 else ii_residcorr[i,ii]=
      (ii_obscorr[i,ii]-ii_predcorr[i,ii])/(1-abs(ii_predcorr[i,ii]))}}
  # Individual fit
  outfit.person = apply(z^2, 1, mean)
  infit.person = apply(V * z^2,1,sum)/ apply(V,1,sum)
  rri_pinrr = matrix(NA, nraw, nitem)
  for (rr in 1:nraw) {
    for(i in 1:nitem) {
      rri_pinrr[rr,i]=sum(p_pinrr*pi_resp[,i]*(p_raw==rr))
    }
  }
  colnames(rri_outobs) = colnames(rri_pinrr) = colnames(XX)
  p_wt=rep(0,npat) #weight calc as p_inrr x rr_totwt
  p_infit=rep(0,npat) #infit
  p_outfit=rep(0,npat) #outfit
  for(p in 1:npat) p_wt[p]=p_pinrr[p]*rr_totwt[p_raw[p]]
  for(p in 1:npat){errsq=0;eesq=0;outfit=0;
                   for(i in 1:nitem){px=rri_pinrr[p_raw[p],i];
                                     if (pi_resp[p,i]==1) errx=1-px else errx=px;
                                     errsq=errsq+errx^2; eesq=eesq+px*(1-px); outfit=outfit+errx^2/(px*(1-px))}
                   p_infit[p]=errsq/eesq; p_outfit[p]=outfit/nitem}
  infit.person.theor = p_infit
  outfit.person.theor = p_outfit
  if(is.null(quantile.seq))
    quantile.seq = seq(0,1,0.01) 
  c_ppoints = quantile.seq
  w.infit = p_wt[order(p_infit)]
  w.outfit = p_wt[order(p_outfit)]
  q.infit.theor = wtd.quantile(sort(p_infit), weights=w.infit, probs=c_ppoints)
  q.infit = wtd.quantile(infit.person, weights=wt2, probs=c_ppoints)
  q.outfit.theor = wtd.quantile(sort(p_outfit), weights=w.outfit, probs=c_ppoints)
  q.outfit = wtd.quantile(outfit.person, weights=wt2, probs=c_ppoints)
  i_sumpq=rep(0,nitem)
  i_sumpqsq=rep(0,nitem)
  i_sumpq=i_sumpq*0
  for (rr in 1:nraw){for (i in 1:nitem){
    i_sumpqsq[i]=i_sumpqsq[i]+
      (rri_pinrr[rr,i]*(1-rri_pinrr[rr,i]))^2*rr_totwt[rr]}}
  for (rr in 1:nraw) {
    for (i in 1:nitem) {i_sumpq[i]=
                          i_sumpq[i]+rri_pinrr[rr,i]*(1-rri_pinrr[rr,i])*rr_totwt[rr]}}
  i_seinfit=sqrt(i_sumpq-4*i_sumpqsq)/i_sumpq
  n.rv.tot = table(factor(rv, levels = 0:k))
  wt.rv.tot = sapply(1:(k+1), function(i) sum(wt[rv==i-1], na.rm = T))
#   wt.rv.tot = sapply(as.numeric(names(n.rv.tot)), function(i) 
#     sum(wt[rv==i], na.rm = T))
  # Weighted N of yes on the complete non-extreme sample
  wt.yes = sapply(1:k, function(i) sum(wt2[XX2[,i]==1], na.rm = T))
  valid.resp = nrow(XX2)
  valid.resp.w = sum(wt2)
  wt.yes.perc = sapply(1:k, function(i) sum(wt2[XX2[,i]==1], na.rm = T)/sum(wt2)
                       * 100)
  # Weighted N of yes on total sample
  wt.yes.tot = sapply(1:k, function(i) sum(wt[XX[,i]==1], na.rm = T))
  wt.yes.perc.tot = sapply(1:k, function(i) sum(wt[XX[,i]==1], na.rm = T)/sum(wt)
                       * 100)
  # N of yes on the complete non-extreme sample
  n.yes = apply(XX2,2,sum)
  perc.yes = n.yes/nrow(XX2)*100
  # N of yes on the total sample
  n.yes.tot = apply(XX,2, function(i) sum(i, na.rm=T))
  perc.yes.tot = n.yes.tot/nrow(XX)*100
  #missing.resp = sum(is.na(XX))
  missing.resp = sum(is.na(rowSums(XX)))
  miss.item.mat = apply(XX, 2,is.na)
  miss.item = apply(miss.item.mat, 2, sum)
  miss.item.w = sapply(1:k, function(i) sum(wt[miss.item.mat[,i]]))
  missing.resp.w = sum(wt[is.na(rowSums(XX))])
  # Prepare data for missing analysis
  #
  ri_in_all=XX
  rs = r_raw_all = rowSums(XX)
  r_nvalid = rowSums(!is.na(XX))
  ri_inz_all=as.matrix(ri_in_all)
  ri_inz_all[is.na(ri_inz_all)]=-1
  ncaseall=nrow(XX)
  r_wt_all=wt
  i_names=colnames(ri_in_all)
  #
  #calculate raw score and distribution by raw score for any valid and all valid
  all.na = apply(XX,1,function(i) sum(is.na(i)))
  # Any
  rs2=rowSums(XX,na.rm=T)
  w.anyv=wt[all.na!=ncol(XX)]
  rr_wt_anyv = tab.weight(as.factor(rs2[all.na!=ncol(XX)]), w.anyv)$tab.ext.w
  rr_ncase_anyv = table(rs2[all.na!=ncol(XX)])
  totwt_anyv = sum(w.anyv)
  perc_wt_anyv = round(totwt_anyv/nrow(XX)*100,1)
  rr_pctwt_anyv = rr_wt_anyv*100/totwt_anyv
  #All
  w.allv=wt[all.na!=ncol(XX)]
  rr_wt_allv = tab.weight(as.factor(rs[all.na!=ncol(XX)]), w.allv)$tab.ext.w
  rr_ncase_allv=unname(table(rs[all.na!=ncol(XX)]))
  totwt_allv=sum(rr_wt_allv)
  rr_pctwt_allv=rr_wt_allv*100/totwt_allv
  # Missing information
  nitem = k
  nvall=nitem+1
  vall_value=rep(0,nvall) #
  vall_value = 0:nitem
  i_max = apply(XX,2,function(x) max(x,na.rm=T))
  rawmax=sum(i_max)
  ri_in_all=XX
  nrawall=rawmax+1 #number of cases including 0
  rrall_value=0:rawmax 
  #calculate distribution of missing: 
  # Distribution of valid responses = da 0 ad 8 quante sono i valori validi
  tab.not.na = apply(XX, 1, function(i) sum(!is.na(i)))
  vall_ncaseinv = sapply(0:nitem, function(i) sum(tab.not.na==i)) #number of cases by number of valid items
  vall_wtcaseinv = sapply(0:nitem, function(i) sum(wt[tab.not.na==i])) #sum of weights of cases by number valid
  
  vall_pctncaseinv=vall_ncaseinv*100/ncaseall
  vall_pctwtcaseinv=vall_wtcaseinv*100/ncaseall #total weight is constrained to ncase
  
  ncaseanyv=ncaseall-vall_ncaseinv[1] #subtract cases w zero valid
  wtcaseanyv=ncaseall-vall_wtcaseinv[1] #subtract wt of cases w zero valid
  
  vall_pctncaseinv_anyv=vall_ncaseinv*100/ncaseanyv
  vall_pctwtcaseinv_anyv=vall_wtcaseinv*100/wtcaseanyv
  vall_pctncaseinv_anyv[1]=0
  vall_pctwtcaseinv_anyv[1]=0
  
  #calculate item missing if any item valid
  # Missing by item if any valid: quanti sono i missing per ogni item
  
  i_nmiss=rep(0,nitem) #number missing if any valid
  i_nvalid=rep(0,nitem) #number of cases w valid response to the item
  i_pctnmiss=rep(0,nitem) #pct missing if any valid
  i_wtmiss=rep(0,nitem) #wt number missing if any valid
  i_wtvalid=rep(0,nitem) #wt number w valid resp to the item
  i_pctwtmiss=rep(0,nitem) #wt pct missing if any valid
  
  i_nmiss=colSums(is.na(XX))
  i_nvalid=ncaseall-i_nmiss
  i_nmiss=i_nmiss-vall_ncaseinv[1]

  for (i in 1:nitem)i_wtmiss[i]=sum((ri_inz_all[,i]==-1)*r_wt_all)
  i_wtvalid=ncaseall-i_wtmiss
  i_wtmiss=i_wtmiss-vall_wtcaseinv[1]
  
  i_pctnmiss=i_nmiss*100/ncaseanyv
  i_pctwtmiss=i_wtmiss*100/wtcaseanyv
  
  #i_nvalid and i_wtvalid already calculated
  ic_numcase_all=mat.or.vec(nitem,4) #responses by value index is value+1
  ic_wtcase_all=mat.or.vec(nitem,4) #wt response by value index is valueC+1
  ic_pctwt_all=mat.or.vec(nitem,4) #wt pct of valid resp for item
  c_catval=c(0,1,2,3) #value of category
  
  for (r in 1:ncaseall) for (i in 1:nitem) if(!is.na(XX[r,i])){
    ic_numcase_all[i,(XX[r,i]+1)]=ic_numcase_all[i,(XX[r,i]+1)]+1
    ic_wtcase_all[i,(XX[r,i]+1)]=ic_wtcase_all[i,(XX[r,i]+1)]+r_wt_all[r]}
  
  for (c in 1:4) ic_pctwt_all[,c]=ic_wtcase_all[,c]*100/i_wtvalid
  
  ##5 subset to complete nonexreme with weight and recalculate wt
  ic_numcase_cnext=mat.or.vec(nitem,4) #responses by value index is value+1
  
  rdf_cnext=data.frame(ri_inz_all,r_wt_all,r_nvalid,r_raw_all)
  
  rdf_cnext=subset(rdf_cnext,(r_nvalid==nitem & r_raw_all>0 & 
                                r_raw_all<rawmax & r_wt_all>0))
  
  nccnext=nrow(rdf_cnext)
  
  #calculate weights at mean 1 for cnext cases
  r_wt=rep(0,nccnext)
  totwtcnext=sum(rdf_cnext$r_wt_all)
  r_wt=rdf_cnext$r_wt_all/totwtcnext*nccnext
  ri_in=as.matrix(rdf_cnext[,1:nitem])
  ri_in[is.na(ri_in)]=-1
  r_raw=as.matrix(rdf_cnext$r_raw_all)
  # r_raw=rs
  totwtcnext=nccnext #redefined totwtcnext for later use
  
  #calculate cnext cases by item by category
  
  for (r in 1:nccnext) for (i in 1:nitem) 
    ic_numcase_cnext[i,(ri_in[r,i]+1)]=ic_numcase_cnext[i,(ri_in[r,i]+1)]+1
  
  ##6 code response into riv matrix with 1 indicating at or above
  riv_in=rep(0,nccnext*nitem*3)
  dim(riv_in)=c(nccnext,nitem,3)
  rriv_wt=rep(0,(rawmax-1)*nitem*3) #sum of wt in or greater than v
  dim(rriv_wt)=c(rawmax-1,nitem,3)
  rriv_wtprop=rep(0,(rawmax-1)*nitem*3) #weighted proportion within raw score
  dim(rriv_wtprop)=c((rawmax-1),nitem,3)
  iv_wt=mat.or.vec(nitem,3) #control totals for solution
  rr_wt=rep(0,rawmax-1)
  rr_ncase=rep(0,rawmax-1)
  
  for (r in 1:nccnext) for (i in 1:nitem) for (v in 1:i_max[i])
    if(ri_in[r,i]>=v) riv_in[r,i,v]=1
  
  #accumulate to raw scores and to totals for cnext
  for (r in 1:nccnext) {rr_wt[r_raw[r]]=rr_wt[r_raw[r]]+r_wt[r];
                        rr_ncase[r_raw[r]]=rr_ncase[r_raw[r]]+1;
                        for (i in 1:nitem) for (v in 1:3){
                          rriv_wt[r_raw[r],i,v]=rriv_wt[r_raw[r],i,v]+r_wt[r]*riv_in[r,i,v];
                          iv_wt[i,v]=iv_wt[i,v]+r_wt[r]*riv_in[r,i,v]}}
  
  #calculate wt ge v as proportion of rr_wt
  rriv_wtprop=rriv_wt/rr_wt #missing if zero rr_wt but no problem

  # Save output
  if(write.file){
    print.title = paste(country, "input data from R datafile")
    write.table(print.title, paste("Output", country,".csv", sep=""), append = F, 
                sep = ",", eol = "\n", na = "NA", dec = ".", col.names=F, row.names=F)
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    valid.cases = rbind("N complete non-extreme" = valid.resp, 
                        "WN complete non-extreme" = valid.resp.w,
                        "N total" = nrow(XX),
                        "N Any Missing" = missing.resp, 
                        "WN Any missing" = missing.resp.w)
    write.table(valid.cases, paste("Output", country,".csv", sep=""), append = T, 
                sep = ",", eol = "\n", na = "NA", dec = ".", col.names = F)
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    item.data = cbind("Item"=colnames(XX),
                      "Severity" = b.w, "SE severity" = se.b.w,
                       "Infit" = infit.w, "SE infit" = i_seinfit,
                      "Outfit" = outfit.w, 
                      "N Yes on complete non-extreme sample" = n.yes, "Perc Yes on complete non-extreme sample"= perc.yes,
                      "WN Yes on complete non-extreme sample" = wt.yes, "WPerc Yes on complete non-extreme sample"= wt.yes.perc,
                      "N Yes on total sample" = n.yes.tot, "Perc Yes on total sample"= perc.yes.tot,
                      "WN Yes on total sample" = wt.yes.tot, "WPerc Yes on total sample"= wt.yes.perc.tot)
#                       "N missing" = miss.item, "W missing" = miss.item.w)
    suppressWarnings(
      write.table(item.data, paste("Output", country,".csv", sep=""), append = T, 
                sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F)
    )
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    if(l.d == 2){
      rs.data = cbind("Raw-score" = 0:k, "Severity" = a.w, "Error" = se.a.w,
                      "N cases" = n.rv.tot, "W cases" = wt.rv.tot)
    }
    if(l.d == 3 & d[2]<1){
      rs.data = cbind("Raw-score" = c(c("0_1","0_2"),1:k), "Severity" = a.w, "Error" = se.a.w,
                      "N cases" = c("",n.rv.tot), "W cases" = c("",wt.rv.tot))
    }
    if(l.d == 3 & d[2]>1){
      rs.data = cbind("Raw-score" = c(0:(k-1), c(paste(k,"_1",sep=""),paste(k,"_2",sep=""))), 
                      "Severity" = a.w, "Error" = se.a.w,
                      "N cases" = c(n.rv.tot, ""), "W cases" = c(wt.rv.tot,""))
    }
    if(l.d == 4){
      rs.data = cbind("Raw-score" = c(c("0_1","0_2"),
                                      1:(k-1), 
                                      c(paste(k,"_1",sep=""),paste(k,"_2",sep=""))), 
                      "Severity" = a.w, "Error" = se.a.w,
                      "N cases" = c("",n.rv.tot,""), "W cases" = c("",wt.rv.tot,""))
    }
    suppressWarnings(
      write.table(rs.data, paste("Output", country,".csv", sep=""), append = T, 
                sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F))
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    reliab.data = rbind("Reliab."=reliab, "Reliab. flat" = reliab.fl)
    suppressWarnings(
      write.table(reliab.data, paste("Output", country,".csv", sep=""), append = T, 
                  sep = ",", eol = "\n", na = "NA", dec = ".", col.names = F))
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)    
    #write missing cases analysis
    vall_value=0:nitem
    vall_df_missing=data.frame(vall_value,vall_ncaseinv,vall_pctncaseinv,vall_wtcaseinv,
                               vall_pctwtcaseinv,vall_pctncaseinv_anyv,vall_pctwtcaseinv_anyv)
    ttext1=paste("\nDistribution of valid responses")
    ttext2=paste("\nNum valid,Num cases,Pct cases,Wt cases,Wt pct,Pct if any valid,Wt pct if any valid")
    vall_dfheader_missing=paste(ttext1,ttext2)
    write(vall_dfheader_missing,file = paste("Output", country,".csv", sep=""), 
          append=T)
    suppressWarnings(write.table(vall_df_missing, file = paste("Output", country,".csv", sep=""), 
                     append = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", 
                     row.names = FALSE, col.names = FALSE))
    #write missing by item if any valid
    fn_output = paste("Output", country,".csv", sep="")
    i_df_missing=data.frame(i_names,i_nmiss,i_pctnmiss,i_wtmiss,i_pctwtmiss)
    i_dfheader_missing=paste("\nMissing by item if any valid\nItem,Num missing,Pct missing, Wt missing, Wt pct missing")
    write(i_dfheader_missing,file=fn_output,append=T)
    suppressWarnings(write.table(i_df_missing, file=fn_output, append = TRUE, 
                                   sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                                   col.names = FALSE))    
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("Residual correlation", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    corr.data = ii_residcorr
    corr.data[lower.tri(corr.data, diag = T)] = ""
    corr.data = cbind(colnames(XX), corr.data)
    corr.data = corr.data[,-2]
    corr.data = corr.data[-k,]
    suppressWarnings(
      write.table(corr.data, paste("Output", country,".csv", sep=""), append = T, 
                  sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F))
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("############# Detailed output#############",
        file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("Observed response proportion", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    add.col = c(1:(k-1))
    obs.prop.mat = cbind("Raw score" = add.col, rri_outobs)
    suppressWarnings(
      write.table(obs.prop.mat, paste("Output", country,".csv", sep=""), append = T, 
                  sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F))
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("Predicted response proportion", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    add.col = c(1:(k-1))
    pred.prop.mat = cbind("Raw score" = add.col, rri_pinrr)
    suppressWarnings(
      write.table(pred.prop.mat, paste("Output", country,".csv", sep=""), append = T, 
                  sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F))
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("Observed and expected respondent infit distribution (weighted)", 
        file = paste("Output", country,".csv", sep=""), append = TRUE)
    #
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    fit.cases.mat = cbind("Percentile" = quantile.seq*100, "Obs infit"=q.infit, 
                            "Pred infit"=q.infit.theor, 
                            "Obs outfit"=q.outfit, "Pred outfit"=q.outfit.theor)
    suppressWarnings(
      write.table(fit.cases.mat, paste("Output", country,".csv", sep=""), append = T, 
                  sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F))
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Output", country,".csv", sep=""), append = TRUE)
    log.lik = round(opt.w$value,3)
    cat(paste("Rasch log-lik", log.lik), 
        file = paste("Output", country,".csv", sep=""), append = TRUE)
  }
	return(list(country = country,
              b = b.w, a = a.w, se.b = se.b.w, se.a = se.a.w, infit = infit.w, 
              outfit = outfit.w, reliab = reliab, reliab.fl = reliab.fl, 
              infit.person = infit.person, infit.person.theor=infit.person.theor,
	            outfit.person = outfit.person,
	            outfit.person.theor = outfit.person.theor,
	            q.infit.theor = q.infit.theor, q.infit = q.infit,
	            q.outfit.theor = q.outfit.theor, q.outfit = q.outfit,
              res.corr = ii_residcorr, se.infit = i_seinfit, mat.res = mat.res,
              d = d, XX=XX, wt = wt, n.compl = valid.resp,
              wt.rs=wt.rv.tot) )
}
