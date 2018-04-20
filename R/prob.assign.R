prob.assign <-
function(rr = NULL,  rwthres = NULL, sthres = NULL, 
         eps.a = NULL, flex = list(a=NULL, se.a=NULL, d=NULL, XX=NULL, wt=NULL)) {
  if(!is.null(rr)){
  a.w = rr$a
  se.a.w = rr$se.a
  d = rr$d
  .data = rr$XX
  wt = rr$wt
  } else if(!is.null(flex)){
    a.w = flex$a
    se.a.w = flex$se.a
    d = flex$d
    .data = flex$XX
    wt = flex$wt
  }
  #
  XX = .data
  rv = rowSums(XX)
  n = nrow(XX)
  k = ncol(XX)
  #
  if (is.null(d)) d = c(0.5,(k-1)+0.5)
  if (is.null(eps.a)) eps.a = 0.001
#   if (is.null(rwthres)) rwthres = c(1,2)
#   if (is.null(sthres)) sthres = c(-1,0,1)
  # Select extremes
  l.a = length(a.w)
  if(l.a == (k+1)) {
    # Eliminating rs 0
    a.w = a.w[-1]
    se.a.w = se.a.w[-1]
    } 
  if(l.a == (k+2)){
    if(d[2] < 1) {
      # Eliminating rs 0
      a1 = c(a.w[1], a.w[3:(k+2)])
      a2 = c(a.w[2:(k+2)])
      #
      se.a1 = c(se.a.w[1], se.a.w[3:(k+2)])
      se.a2 = c(se.a.w[2:(k+2)])
      # Eliminating rs 0
      a1 = a1[-1]
      a2 = a2[-1]
      se.a1 = se.a1[-1]
      se.a2 = se.a2[-1]
      #
      a.list = list(a1, a2)
      se.a.list = list(se.a1, se.a2)
      #
      a.w = a.list
      se.a.w = se.a.list
    } else {
      a1 = c(a.w[1:k], a.w[(k+2)])
      a2 = c(a.w[1:(k+1)])
      #
      se.a1 = c(se.a.w[1:k], se.a.w[(k+2)])
      se.a2 = c(se.a.w[1:(k+1)])
      # Eliminating rs 0
      a1 = a1[-1]
      a2 = a2[-1]
      se.a1 = se.a1[-1]
      se.a2 = se.a2[-1]
      #
      a.list = list(a1, a2)
      se.a.list = list(se.a1, se.a2)
      #
      a.w = a.list
      se.a.w = se.a.list
    } 
  } else if(l.a == (k+3)) {
    a1 = c(a.w[2:(k+2)])
    a2 = c(a.w[1], a.w[3:(k+1)], a.w[(k+3)])
    a3 = c(a.w[1], a.w[3:(k+2)])
    a4 = c(a.w[2:(k+1)], a.w[(k+3)])
    #
    se.a1 = c(se.a.w[2:(k+2)])
    se.a2 = c(se.a.w[1], se.a.w[3:(k+1)], se.a.w[(k+3)])
    se.a3 = c(se.a.w[1], se.a.w[3:(k+2)])
    se.a4 = c(se.a.w[2:(k+1)], se.a.w[(k+3)])
    #
    # Eliminating rs 0
    a1 = a1[-1]
    a2 = a2[-1]
    a3 = a3[-1]
    a4 = a4[-1]
    se.a1 = se.a1[-1]
    se.a2 = se.a2[-1]
    se.a3 = se.a3[-1]
    se.a4 = se.a4[-1]
    #
    a.list = list(a1, a2, a3, a4)
    se.a.list = list(se.a1, se.a2, se.a3, se.a4)
    #
    a.w = a.list
    se.a.w = se.a.list
  }
  # Absolute frequency distribution for each raw score
  rvv = sort(rv)
  n_j = sapply(0:(k), function(i) sum(wt[rv == i], na.rm = T))
  # Discrete frequency distribution for each raw score 
  f_j = n_j/sum(n_j)
  f_j = f_j[-1]
  l_j = length(f_j)
  if(!is.null(rwthres)) n.rwthres = length(rwthres)
  if(!is.null(sthres)) n.sthres = length(sthres)
  if(!is.null(rwthres)){
  p = sapply(1:n.rwthres, function(i) sum(f_j[rwthres[i]:l_j]))
  #Probabilistic Assignment with thresholds defined in terms of raw scores
  assig = function(rwthres, a.w, se.a.w, f_j, eps = NULL) {	
    if (is.null(eps)) eps = 0.001
    if(l.a == (k+1)){
      prob = sapply(1:n.rwthres, 
                    function(i) (1-pnorm(rwthres[i],a.w,se.a.w))*f_j)
      diff = apply(prob,2,sum)-p
      while (abs(min(diff)) > eps) {
        prob = sapply(1:n.rwthres, 
                      function(i) (1-pnorm(rwthres[i],a.w,se.a.w))*f_j)
        diff = apply(prob,2,sum)-p		
        thres.new = c(rwthres+diff)		
        rwthres = thres.new
      } } else {
        prob = list()
        for(j in 1:length(a.w)){
          prob[[j]] = sapply(1:n.rwthres, 
                             function(i) (1-pnorm(rwthres[i],a.w[[j]],se.a.w[[j]]))*f_j) 
          diff = apply(prob[[j]],2,sum)-p
          while (abs(min(diff)) > eps) {
            prob[[j]] = sapply(1:n.rwthres, 
                               function(i) (1-pnorm(rwthres[i],a.w[[j]],se.a.w[[j]]))*f_j)
            diff = apply(prob[[j]],2,sum)-p  	
            thres.new = c(rwthres+diff)		
            rwthres = thres.new
          }
        }
      }
    return(list("thres"=rwthres, "prob" = prob))	
  }
  res = assig(rwthres, a.w, se.a.w, f_j)
  }
  if(!is.null(sthres)){
  #Probabilistic assignment with pre-defined thresholds
  if(l.a == (k+1)){
    sprob = colSums(sapply(1:n.sthres,
                           function(i) (1-pnorm(sthres[i],a.w,se.a.w))*f_j))
  } else {
    sprob = list()
    for(j in 1:length(a.w)){
      sprob[[j]] = colSums(sapply(1:n.sthres,
                                  function(i) (1-pnorm(sthres[i],a.w[[j]],se.a.w[[j]]))*f_j))
    }
  }}
  if(is.null(sthres)) sprob=NULL
  if(is.null(rwthres)) thres=NULL else thres=res$thres
  if(is.null(rwthres)) f=NULL else f=res$prob
  if(is.null(rwthres)) p=NULL 
  return(list("sprob" = sprob, "thres" = thres, "f" = f, "p" = p,
              f_j = f_j))
}
