equating.fun = function(rr1, st=NULL, tol = .35, spec.com1 = 1:8, 
                           spec.com2=1:8, thres = NULL, maxuniq=3, write.file=FALSE, 
                            plot=FALSE, iterative=TRUE, excl.prior1=NULL, excl.prior2=NULL, wt.spec=NULL){
  # spec.com1 is the set of a priori common items for the country
  # spec.com2 is the set of a priori common items for the standard
  tol2=tol
  b1=rr1$b
  se1=rr1$se.b
  # If null, the standard is set as the VoH 2014-2016 global standard
  if(is.null(st)) st = c(-1.2230564, -0.8471210, -1.1056616,  0.3509848, -0.3117999,  0.5065051,  0.7546138, 1.8755353)
  b.tot = st
  # Set thresholds as 5th and 8th item on the global standard if not specified
  if(is.null(thres)) {
    thres=b.tot[c(5,8)]
    truetresh=T
  } else{
    truetresh=F
  }
  # We start from the assumption that all items (specified in spec.comm1) 
  # are common
  common = rep(F,length(b1))
  common[spec.com1]=T
  plus=sum(!common)
  common2 = rep(F,length(b.tot))
  common2[spec.com2]=T
  plus2=sum(!common2)
  # Standardize everything to mean and SD of a priori common items on the glob.st. scale
  b.st1 = (b1-mean(b1[spec.com1]))/sd(b1[spec.com1])*sd(b.tot[spec.com2])+mean(b.tot[spec.com2])
  a = 1
  oldcommon = common
  oldcommon2 = common2
  if(iterative){
    while(a<=length(b.st1)){
      oldcommon = common
      oldcommon2 = common2
      diff = rep(100, length(b1))
      diff[spec.com1] = abs(b.st1[spec.com1] - b.tot[spec.com2])
      diff2 = rep(100, length(b.tot))
      diff2[spec.com2] = abs(b.st1[spec.com1] - b.tot[spec.com2])
      ord = order(diff, decreasing = T)
      w = ord[a+plus]
      ord2 = order(diff2, decreasing = T)
      w2 = ord2[a+plus2]
      if(diff[w]>=tol)
      {
        common[w] = F}  else {
          common[w] = T}  
      if(diff2[w2]>=tol)
      {
        common2[w2] = F}  else {
          common2[w2] = T} 
      scale1 = sd(b.tot[common2])/sd(b.st1[common])
      shift1 = mean(b.tot[common2])-mean(b.st1[common])*scale1
      b.st1 = shift1 + b.st1 * scale1
  #     scale2 = sd(b.tot[common2])/sd(b.st1[common])
  #     shift2 = mean(b.tot[common2])-mean(b.st1[common])*scale1
  #     b.st1 = shift1 + b.st1 * scale1
      if(sum(oldcommon == common)==length(b.st1) | sum(!common)>maxuniq) break else
      {a = a + 1}
    }
} 
  # Final iteration
if(!iterative) {
  common[excl.prior1]=F
  common2[excl.prior2]=F
}
  scale = sd(b.tot[common2])/sd(b1[common])
  shift = mean(b.tot[common2])-mean(b1[common])*scale
  newthres = (thres-shift)/scale
  a = rr1$a
  se.a = rr1$se.a
  d = rr1$d
# 
  # Max rs
  kk = length(b1)+1
  if(length(d)>2){
    if(d[2]<0){
      a = a[c(1,2:(kk+1))]
      se.a = se.a[c(1,2:(kk+1))]
    } else{
      a = a[-kk]
      se.a = se.a[-kk]
    }    
  }
  prevs.rs = matrix(NA, length(a), length(newthres))
  wt=rr1$wt
  if(!is.null(wt.spec)) wt = wt.spec
  XX=rr1$XX
  rv=rowSums(XX)
  k=ncol(XX)
  n_j = sapply(0:(k), function(i) sum(wt[rv == i], na.rm = T))
  f_j = n_j/sum(n_j)
  f_j[1] = 0
  for (i in 1:length(newthres)) {
    prevs.rs[,i] =  (1-pnorm(newthres[i], a, se.a))*f_j
  }
  prevs = colSums(prevs.rs)
  x.var=st
  y.var=shift+rr1$b*scale
  probs.rs=prevs.rs
  for (i in 1:length(newthres)) {
    probs.rs[,i] =  (1-pnorm(newthres[i], a, se.a))
  }
  probs.rs[1,]=0
  row.names(probs.rs)=as.numeric(names(table(rv)))
  if(truetresh) colnames(probs.rs) = c("FI_mod+", "FI_sev")
  if(plot){
    pdf("Equating_plot.pdf")
    range=range(x.var[spec.com2],y.var[spec.com1])
    plot(x.var[spec.com2],y.var[spec.com1],pch=16,col=2,xlab="Standard",ylab=rr1$country,
         xlim=range,ylim=range)
    abline(c(0,1))
    points(x.var[common2],y.var[common],pch=16,col="blue")
    text(x.var[spec.com2],y.var[spec.com1],colnames(rr1$XX)[spec.com1],cex=.6,pos=1)
    legend("topleft", pch=16, col=c("blue","red"), legend=c("Common","Unique"),
           cex=.7, bty="n", x.intersp=.5)
    dev.off()
  }
  cor.comm.items=cor(x.var[common2],y.var[common])
  names(common)=colnames(rr1$XX)
  # rownames(prevs.rs)=0:ncol(rr1$XX)
# Name thresholds and prevs if thresholds are official
if(truetresh){
  names(prevs)=c("FI_mod+", "FI_sev")
  names(newthres)=c("FI_mod+", "FI_sev")
}
# Write output file
  if(write.file){
    country=rr1$country
    print.title = paste(country, "equating output from R processing")
    write.table(print.title, paste("Equating_Output", country,".csv", sep=""), append = F, 
                sep = ",", eol = "\n", na = "NA", dec = ".", col.names=F, row.names=F)
    # Common items
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("Common items", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    comm.mat=common
    comm.mat[!common]="Unique"
    comm.mat[common]="Common"
    write.table(comm.mat, paste("Equating_Output", country,".csv", sep=""), append = T, 
                sep = ",", eol = "\n", na = "NA", dec = ".", col.names = F)
    # Correlation between common items
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("Correlation between common items", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    write.table(cor.comm.items, paste("Equating_Output", country,".csv", sep=""), append = T, 
                sep = ",", eol = "\n", na = "NA", dec = ".", col.names = F, row.names = F)
    # Adjusted thresholds
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("Adjusted thresholds on country metric", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    adj.thres=newthres
    write.table(adj.thres, paste("Equating_Output", country,".csv", sep=""), append = T, 
                sep = ",", eol = "\n", na = "NA", dec = ".", col.names = F)
    # Prevalence rates at adjusted thresholds
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("Prevalence rates at adjusted thresholds (%)", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    write.table(round(prevs*100,2), paste("Equating_Output", country,".csv", sep=""), append = T, 
                sep = ",", eol = "\n", na = "NA", dec = ".", col.names = F)
    # Probability of being beyond thresholds at each raw score
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("Probability at adjusted thresholds for each raw score", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    cat("\n", file = paste("Equating_Output", country,".csv", sep=""), append = TRUE)
    write.table(probs.rs, paste("Equating_Output", country,".csv", sep=""), append = T, 
                sep = ",", eol = "\n", na = "NA", dec = ".", col.names = F)
  }
  common.print=common
  common.print[common.print==TRUE] = "Common"
  common.print[common.print==FALSE] = "Unique"
  return(list("scale" = scale, "shift" = shift, 
              "common"=common.print, prevs = prevs, probs.rs=probs.rs,
              cor.comm.items=cor.comm.items, adj.thres=newthres, standard=st))  
}
