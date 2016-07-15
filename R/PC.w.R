# PCM.w partial credit model weighted

# PCM.w estimates polyotmous Rash partial credit model up to 4 cats
# TBN: recoding of XX matrix should be done previosly as 0=No,1=Yes for dich items,
# and 0-1-2-3-4 for polytomous items
# In the FIES polytomous, XX should have 8 columns and the last 2 should be coded as 
# 0-1-2...

## HELP

## XX
# Data matrix coded as 0=No, 1=Yes for dichotomous items and
# 0=No, 1=Only once or twice, 2=In some months but not every month, 3=Almost every month
# For the extended FIES, XX will need to have 8 columns. The first 6 are 0/1, and the last 2
# are 0/1/2/3 coded (hungry and whlday)

## wt
# Sampling weights. sum(wt) = sample size. If left unspecified, all the sampling units
# will be weighted in the same way

## extr
# Assumption on the extreme raw scores.It has to be specified as a vector. The first element is 
# the assumption on the 0 raw score, the second element is the assumption on the maximum raw score.
# For the extended FIES, raw score varies between 0 and 10.

## maxiter and minconv
# Convergence criteria

## country
# Country or data name

## write.file
# If TRUE, an output csv file will be written in the working directory

## recode
# This is a mandatory argument. Can be 0, 1, 2 or 3. recode = 1 aggregates "Almost every month" and 
# "In some months but not every month". recode=2  aggregates "In some months but not every month" and
# "Only once or twice". recode=3 aggregates "Never" and "Only once or twice".
# recode = 0 does no aggregation.

## write.iteration
# If TRUE, an information csv file on the iteration process will be written in the working directory 

PC.w=function(XX, wt=NULL, extr=NULL, maxiter=100,minconv=.00001,country=NULL,write.file=F,
            recode = 0, write.iteration=F){
  fn_rdata_in=XX
  if(is.null(wt)) wt = rep(1, nrow(XX))

###sections
##1 recoding
##2 read input data file
##3 recode input data
##4 calc missing information and overview information
##5 subset to complete nonexreme with weight and recalculate wt
##6 code response into riv matrix with 1 indicating at or above
##7 create model response matrix
##8 iterate to convergence
##9 calculate item-by-threshold estimation errors
##10 calculate rasch-thurstone parameters and fit statistics
#11 calculate overall item fit stats
##12 calculate raw score parameters and errors
##13 calculate Rasch reliability weighted by raw score and flat
##14 calculate conditional correlations
##15 write output file and R output list


##1 read recode file
nitem=ncol(XX)
#names are v0 to v5 and represent input values
#  but index values are 1 to 6
#  cell values are recoded values, with -1 missing

#calculate max v for each item and raw score max
i_max = apply(XX,2,function(x) max(x,na.rm=T))
rawmax=sum(i_max)

##2 read input data file
rdf_indata=fn_rdata_in
ncaseall=nrow(XX)

#extract item responses
ri_in_all=rdf_indata
i_names=colnames(ri_in_all)

#calc case weights
r_wt_all=wt
totwtall=sum(wt) #for no wt and redefine for possible future use

##3 recode input data and calculate distr of raw if any valid and all valid
which.pol = which(i_max>1)
n.pol = length(which.pol)
var.pol = XX[,which.pol]
# Max raw score on un-recoded data
# i_max = apply(XX,2,function(x) max(x,na.rm=T))
# rawmax=sum(i_max)

# Polytomous categories:
# Almost every month (3), During some months but not every month (2), 
# Only one or two times (1)
# recode=0 --> No aggregation
# recode=1 --> Almost every month (3) + During some months but not every month (2)
# recode=2 --> During some months but not every month (2) + Only one or two times (1)
# recode=3 --> aggregates "Never" and "Only once or twice"

if(n.pol!=0){
  if(recode == 1){
    for(i in 1:n.pol){
      v = var.pol[,i]
      v[v==3]=2
      XX[,which.pol[i]] = v
    } 
  } else if(recode == 2){
    for(i in 1:n.pol){
      v = var.pol[,i]
      v[v==2]=1
      v[v==3]=2
      XX[,which.pol[i]] = v
    } 
  } else if(recode == 3){
    for(i in 1:n.pol){
      v = var.pol[,i]
      v[v==1]=0
      v[v==2]=1
      v[v==3]=2
      XX[,which.pol[i]] = v
    } 
  }
} else 
{cat("Data include only dichotomous variables. Use the RM.w() function instead.")
break}

# Re-calculate max v for each item and raw score max on recoded data
i_max = apply(XX,2,function(x) max(x,na.rm=T))
rawmax=sum(i_max)
ri_in_all=XX
rr_value=as.numeric(names(table(rowSums(XX))))
rawmax=max(rr_value)

# Pseudo-raw score for extreme parameter calculation
if(!is.null(extr)){
  extr_lo=extr[1]
  extr_hi=extr[2]
} else {
  extr_lo=.5
  extr_hi=rawmax-0.5
}

# Raw score
rs = r_raw_all = rowSums(XX)
# r_nvalid = rs
# r_nvalid[is.na(rs)]=-1
r_nvalid = rowSums(!is.na(XX))
# r_nvalid[is.na(rs)]=-1
# 
ri_inz_all=as.matrix(ri_in_all)
ri_inz_all[is.na(ri_inz_all)]=-1
# ri_inz_all[,]=0
# r_nvalid=rep(0,ncaseall)
# r_raw_all=rep(0,ncaseall)
# # 
# r_nvalid = ri_inz_all+1
# r_nvalid[r_nvalid==NA]=-1
# 
# for (r in 1:ncaseall) {for (i in 1:nitem){x=ri_in_all[r,i];
#  if (is.na(x)) ri_inz_all[r,i]=-1 else{
#  if (x>=0 & x<=5) {ri_inz_all[r,i]=iv_recode[i,(x+1)]; 
#  if (ri_inz_all[r,i] >= 0) r_nvalid[r]=r_nvalid[r]+1} else ri_inz_all[r,i]=-1}}}
# 

# # 
# for (r in 1:ncaseall) if (r_nvalid[r]>0) 
#   r_raw_all[r]=sum(ri_inz_all[r,]*(ri_inz_all[r,]>=0))
# #corrected line above in 2014-09-11 revision
# 
# nrrall=rawmax+1
# rr_value_all=seq(0,rawmax) #values 0 to rawmax w index = value+1
# rr_ncase_anyv=rep(0,nrrall) #number of cases by raw score if any valid
# rr_wt_anyv=rep(0,nrrall) #weighted cases by raw score if any valid
# rr_pctwt_anyv=rep(0,nrrall) #weighted percent cases by raw score if any valid
# rr_ncase_allv=rep(0,nrrall) #number of cases by raw score if all valid
# rr_wt_allv=rep(0,nrrall) #weighted cases by raw score if all valid
# rr_pctwt_allv=rep(0,nrrall) #weighted percent cases by raw score if all valid
# 
# for (r in 1:ncaseall) {
#   if (r_nvalid[r]>0) {rr_ncase_anyv[r_raw_all[r]+1]=rr_ncase_anyv[r_raw_all[r]+1]+1;
#     rr_wt_anyv[r_raw_all[r]+1]=rr_wt_anyv[(r_raw_all[r]+1)]+r_wt_all[r]};
#   if (r_nvalid[r]==nitem){rr_ncase_allv[r_raw_all[r]+1]=rr_ncase_allv[r_raw_all[r]+1]+1;
#     rr_wt_allv[r_raw_all[r]+1]=rr_wt_allv[r_raw_all[r]+1]+r_wt_all[r]}}
# 
# totwt_anyv=sum(rr_wt_anyv)
# rr_pctwt_anyv=rr_wt_anyv*100/totwt_anyv
# 
# totwt_allv=sum(rr_wt_allv)
# rr_pctwt_allv=rr_wt_allv*100/totwt_allv

#calculate raw score and distribution by raw score for any valid and all valid
# Note: Raw score should vary from 0 to 10
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

##4 calc missing information and overview information
#    wt by raw score for all cases w any valid
#  NOTE: index values are 1 less than raw score for rrall variables
#    and 1 less than number valid for nvall variables
#    use value variables for identifying these

nvall=nitem+1
vall_value=rep(0,nvall) #
vall_value = 0:nitem
# for (i in 1:nitem) vall_value[i+1]=i

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

# i_nmiss=colSums(ri_inz_all==-1)
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

##7 create model response matrix
npat=1
for (i in 1:nitem) npat=npat*(i_max[i]+1)
npat=npat-2

#assign place values for columns in response matrix
i_ival=rep(0,nitem) #place value for indexing response matrix
for (i in 1:(nitem-1)){i_ival[i]=1; for (ii in (i+1):nitem) 
  i_ival[i]=i_ival[i]*(i_max[ii]+1)}
i_ival[nitem]=1

#assign values to items in response matrix and calc raw score
pi_x=mat.or.vec(npat,nitem)
for (p in 1:npat) {remval=p; for 
  (i in 1:nitem) {pi_x[p,i]=trunc(remval/i_ival[i]); 
  remval=remval-pi_x[p,i]*i_ival[i]}}

p_raw=rowSums(pi_x)


##8 iterate to convergence
#   note: parameters are reverse-signed until after iteration complete
it_param=mat.or.vec(nitem,3) #initial parameter values
it_paramrev=mat.or.vec(nitem,3) #revised parameter values

#to check convergence
nit_param=rep(0,nitem*3*100) #to see convergence
dim(nit_param)=c(100,nitem,3)

ic_lprob=mat.or.vec(nitem,4) #unconditional log-prob for category given theta 0
#  model log-probability matrix, VALUE IS 1 LESS THAN INDEX to accomodate zero

p_prob=rep(0,npat) #unconditional probability
rr_totprob=rep(0,rawmax-1) #sum of pattern probbilities in the raw score
p_probinraw=rep(0,npat) #probability of pattern given raw score
rriv_prob=rep(0,((rawmax-1)*nitem*3)) #by raw score, prob i>=v
dim(rriv_prob)=c(rawmax-1,nitem,3)
iv_totprob=mat.or.vec(nitem,3) #total accumulated probability weighted
iv_totpq=mat.or.vec(nitem,3) #total accumulated pq ex where either is zero
i_meant=rep(0,nitem) #mean threshold values for the item
i_meantrev=rep(0,nitem) #mean revised threshold values for the item

niterate=0
maxadj=1

#main iteration not indented
while (niterate<maxiter & maxadj>minconv){
niterate=niterate+1
it_param=it_paramrev

#calculate log-probability for item-category at theta zero
for(i in 1:nitem) {k=exp(it_param[i,1]); l=exp(it_param[i,2]); 
  m=exp(it_param[i,3]); if (i_max[i]==1)
  {px1=1/(1+1/k); ic_lprob[i,2]=log(px1); ic_lprob[i,1]=log(1-px1)}
  else if (i_max[i]==2) {px0=1/(1+k+k*l); px1=px0*k; px2=px1*l;
  ic_lprob[i,1]=log(px0); ic_lprob[i,2]=log(px1); ic_lprob[i,3]=log(px2)}
  else if(i_max[i]==3) {px0=1/(1+k+k*l+k*l*m); px1=px0*k; px2=px1*l; px3=px2*m
  ic_lprob[i,1]=log(px0); ic_lprob[i,2]=log(px1); ic_lprob[i,3]=log(px2)
  ic_lprob[i,4]=log(px3)}}

#calculate probability of each response pattern at theta zero
for (p in 1:npat) {logresp=0 #to accumulate log prob
  for (i in 1:nitem) logresp=logresp+ic_lprob[i,(pi_x[p,i]+1)]
  p_prob[p]=exp(logresp)}
  
#total probability of patterns in each raw score
for (rr in 1:rawmax-1) rr_totprob[rr]=sum((p_raw==rr)*p_prob)

#calculate probability given raw score
for (p in 1:npat) p_probinraw[p]=p_prob[p]/rr_totprob[p_raw[p]]

#calculate probability of each item >= v (category gt 0) by raw score
rriv_prob[]=0
for (p in 1:npat) for (i in 1:nitem) for (v in 1:i_max[i]) if (pi_x[p,i]>=v)
  rriv_prob[p_raw[p],i,v]=rriv_prob[p_raw[p],i,v]+p_probinraw[p];

#accumulate total item x category (gt 0) weighted by raw score
iv_totprob[]=0
for (rr in 1:(rawmax-1)) for (i in 1:nitem) for (v in 1:i_max[i])
  iv_totprob[i,v]=iv_totprob[i,v]+rriv_prob[rr,i,v]*rr_wt[rr]
 
#accumulate total pq where p ne 1 or 0
iv_totpq[]=0
for (rr in 1:(rawmax-1)) for (i in 1:nitem) for (v in 1:i_max[i])
  if (rriv_prob[rr,i,v]>0 & rriv_prob[rr,i,v]<1)
  iv_totpq[i,v]=iv_totpq[i,v]+rriv_prob[rr,i,v]*(1-rriv_prob[rr,i,v])*rr_wt[rr]

#calculate adjustment
it_adj=mat.or.vec(nitem,3)
it_paramrev=mat.or.vec(nitem,3) #revised parameters

for (i in 1:nitem) for (t in 1:i_max[i]) {v=t;
  it_adj[i,t]=(iv_wt[i,v]-iv_totprob[i,v])/iv_totpq[i,v];
  it_adj[is.na(it_adj)]=0
  if(it_adj[i,t]<(-.6)) it_adj[i,t]=-.6; #limits initial adjustment to avoid crash
  if(it_adj[i,t]>.6) it_adj[i,t]=.6;
  it_paramrev[i,t]=it_param[i,t]+it_adj[i,t]*.8}

#adjust to mean item zero
for (i in 1:nitem) {i_meant[i]=sum(it_param[i,])/i_max[i];
  i_meantrev[i]=sum(it_paramrev[i,])/i_max[i]}

paramrevmean=mean(i_meantrev)

for (i in 1:nitem) for (t in 1:i_max[i]) 
  it_paramrev[i,t]=it_paramrev[i,t]-paramrevmean

if (niterate<=100) nit_param[niterate,,]=it_paramrev

#calculate max adjustment
maxadj=0
for (i in 1:nitem) for (t in 1:i_max[i]) 
  if (abs(it_paramrev[i,t]-it_param[i,t])>maxadj) 
  maxadj=abs(it_paramrev[i,t]-it_param[i,t])
} #end of main iteration

#post-iteration calculations
it_param=-1*it_paramrev
p_logprob=log(p_prob)
minlogpprob=(min(p_logprob))
npatnotcalc=sum(p_logprob<(-34.5)) #beyond computer precision
propnotcalc=npatnotcalc/npat
rm(p_logprob)

#write iteration history to file in mrasch/mrwork
# write.table(nit_param, file="c:/mrasch/mrwork/mrpcm_ckiter.csv", append = F, 
if(write.iteration){
  write.table(nit_param, file="Iteration history.csv", append = F, 
              sep = ",", eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = TRUE)
}

##9 calculate item-by-threshold estimation errors
it_error=mat.or.vec(nitem,3) #estimation error item x threshold parameter

for (i in 1:nitem) for (t in 1:i_max[i]) it_error[i,t]=1/sqrt(iv_totpq[i,t]) 


##10 calculate rasch-thurstone parameters and fit statistics
it_rtparam=mat.or.vec(nitem,3) #rasch-thurstone 
it_rtinfit=mat.or.vec(nitem,3) #rasch-thurstone infit
it_rtoutfit=mat.or.vec(nitem,3) #rasch-thurstone outfit
it_rtniterate=mat.or.vec(nitem,3) #number of iterations for r-t threshold
it_rtadj=mat.or.vec(nitem,3) #final adjustment for r-t threshold

#calculate rasch-thurstone parameters
for (i in 1:nitem) for (t in 1:i_max[i]){
  rtadj=1; rtniterate=0; rtparamrev=it_param[i,t] #start here
  while (rtniterate<maxiter & rtadj>minconv){rtniterate=rtniterate+1;
  rtparam=rtparamrev
  #calculate category probabilities
    if (i_max[i]==1) pxtest=1/(1+1/exp(rtparam-it_param[i,1]))
    else if (i_max[i]==2) {k=exp(rtparam-it_param[i,1]);
    l=exp(rtparam-it_param[i,2]); px0=1/(1+k+k*l); px1=k*px0; px2=l*px1;
    if (t==1) pxtest=px1+px2 else if (t==2) pxtest=px2}
    else if (i_max[i]==3) {k=exp(rtparam-it_param[i,1]);
    l=exp(rtparam-it_param[i,2]); m=exp(rtparam-it_param[i,3]);
    px0=1/(1+k+k*l+k*l*m); px1=k*px0; px2=l*px1; px3=m*px2;
    if (t==1) pxtest=1-px0 else if (t==2) pxtest=px2+px3 else if (t==3) pxtest=px3}
    pqxtest=pxtest*(1-pxtest)
    rtadj=(pxtest-.5)/pqxtest
    if (rtadj>1) rtadj=1
    if (rtadj<(-1)) rtadj=-1
    rtparamrev=rtparam-.5*rtadj
    rtadj=abs(rtadj)
    }
  it_rtparam[i,t]=rtparamrev
  it_rtniterate[i,t]=rtniterate  
  it_rtadj[i,t]=rtadj
  }

#calculate rasch-thurstone-based fit stats
for (i in 1:nitem) for (t in 1:i_max[i]){errsq=0; eesq=0; outfit=0; fitwt=0;
#fitwt needs to be calculated because some categories are 0 or 1 given raw
  for (rr in 1:(rawmax-1)) 
  if (rr_wt[rr]>0 & rriv_prob[rr,i,t]>0 & rriv_prob[rr,i,t]<1){ 
    fitwt=fitwt+rr_wt[rr]; px=rriv_wtprop[rr,i,t]; epx=rriv_prob[rr,i,t];
    errsqx=px*(1-epx)^2+(1-px)*epx^2; eesqx=epx*(1-epx); 
    errsq=errsq+errsqx*rr_wt[rr]; eesq=eesq+eesqx*rr_wt[rr];
    outfit=outfit+errsqx/eesqx*rr_wt[rr]}
  it_rtinfit[i,t]=(errsq/eesq)
  it_rtoutfit[i,t]=(outfit/fitwt)}

#11 calculate overall item fit stats
i_infit=rep(0,nitem)
i_outfit=rep(0,nitem)
i_errsq=rep(0,nitem) #error squared x wt
i_eesq=rep(0,nitem) #expected error squared x wt
rriv_obspropinv=rep(0,((rawmax-1)*nitem*3)) #observed wt prop in v
#  note: zero category not needed
dim(rriv_obspropinv)=c((rawmax-1),nitem,3)
rriv_predpropinv=rep(0,((rawmax-1)*nitem*3)) #predicted prop in v
#  note: zero category not needed
dim(rriv_predpropinv)=c((rawmax-1),nitem,3)

#calculate obs wt prop and pred prop in each pos value by raw score and item
for (rr in 1:(rawmax-1)) if (rr_wt[rr]>0) for (i in 1:nitem) for (v in 1:i_max[i]){
  rriv_obspropinv[rr,i,v]=rriv_wtprop[rr,i,v]; 
  rriv_predpropinv[rr,i,v]=rriv_prob[rr,i,v];
  if (v<i_max[i]){
    rriv_obspropinv[rr,i,v]=rriv_obspropinv[rr,i,v]-rriv_wtprop[rr,i,v+1];
    rriv_predpropinv[rr,i,v]=rriv_predpropinv[rr,i,v]-rriv_prob[rr,i,v+1]}}
 
#accumulate fit components
for (rr in 1:(rawmax-1)) if (rr_wt[rr]>0) for (i in 1:nitem) {
  obsvz=1-rriv_wtprop[rr,i,1];
  predvz=1-rriv_prob[rr,i,1]; mnexpval=0; errsq=0; eesq=0;
  for (v in 1:i_max[i]) mnexpval=mnexpval+rriv_predpropinv[rr,i,v]*v;
  errsq=errsq+obsvz*mnexpval^2; #errsq if value 0
  eesq=eesq+predvz*mnexpval^2; #expected errsq if value 0
  for (v in 1:i_max[i]) {errsq=errsq+rriv_obspropinv[rr,i,v]*(v-mnexpval)^2;
    eesq=eesq+rriv_predpropinv[rr,i,v]*(v-mnexpval)^2}
  i_errsq[i]=i_errsq[i]+rr_wt[rr]*errsq;
  i_eesq[i]=i_eesq[i]+rr_wt[rr]*eesq;
  i_outfit[i]=i_outfit[i]+rr_wt[rr]*errsq/eesq}

#calculate infit and outfit
i_infit=i_errsq/i_eesq
i_outfit=i_outfit/totwtcnext

##12 calculate raw score parameters and errors
rr_value=as.numeric(names(table(rowSums(XX))))
# rr_param=rep(0,(rawmax+1)) #raw score parameter
# rr_error=rep(0,(rawmax+1)) #raw score error
# rr_niterate=rep(0,(rawmax+1)) #number of iterations to calc raw score param
rr_param=rr_error=rr_niterate=rep(0,length(rr_value))
rr_finaladj=rep(1,length(rr_value)) #final adjustment to raw score param
#rr_value=seq(0,rawmax)
rr_value[1]=extr_lo
if(extr_hi==99) rr_value[length(rr_value)]=rawmax-.5 else rr_value[length(rr_value)]=extr_hi

for (rr in 1:length(rr_value)) {paramrev=0;                    
  while (rr_niterate[rr]<maxiter & rr_finaladj[rr]>minconv){
    rr_finaladj[is.infinite(rr_finaladj)]=1e5
    rr_niterate[rr]=rr_niterate[rr]+1; param=paramrev; predraw=0; eesq=0; 
    for (i in 1:nitem) {
      for (t in 1:i_max[i]){
        if (t==1) k=exp(param-it_param[i,t])
        else if (t==2) l=exp(param-it_param[i,t])
        else if (t==3) m=exp(param-it_param[i,t])}
      if (i_max[i]==1) {
        px1=1/(1+1/k); px0=1-px1; predraw=predraw+px1; eesq=eesq+px0*px1}
      else if (i_max[i]==2){
        px0=1/(1+k+k*l); px1=k*px0; px2=l*px1
        mnx=px1+px2*2; predraw=predraw+mnx
        eesq=eesq+px0*mnx^2+px1*(1-mnx)^2+px2*(2-mnx)^2}
      else if (i_max[i]==3){
        px0=1/(1+k+k*l+k*l*m); px1=k*px0; px2=l*px1; px3=m*px2
        mnx=px1+px2*2+px3*3; predraw=predraw+mnx
        eesq=eesq+px0*mnx^2+px1*(1-mnx)^2+px2*(2-mnx)^2+px3*(3-mnx)^2}}
    rradj=(predraw-rr_value[rr])/eesq;
    rr_finaladj[rr]=abs(rradj);
    paramrev=param-rradj*.8}
  rr_param[rr]=paramrev
  rr_error[rr]=1/sqrt(eesq)}

##13 calculate Rasch reliability weighted by raw score and flat

#calculate mean respondent parameter for cnext
rr_param_cnext=rr_param[2:rawmax]
rr_error_cnext=rr_error[2:rawmax]

mnrrparam=sum(rr_wt*rr_param_cnext)/totwtcnext
moderrsq=sum(rr_wt*(rr_param_cnext-mnrrparam)^2)
errsq=sum(rr_wt*rr_error_cnext^2)
raschreliab=round(moderrsq/(moderrsq+errsq),4)
mnrrparamf=mean(rr_param_cnext)
moderrsqf=sum((rr_param_cnext-mnrrparamf)^2)
errsqf=sum(rr_error_cnext^2)
raschreliabf=round(moderrsqf/(moderrsqf+errsqf),4)


##14 calculate conditional correlations
#item by item matrices
ii_obscorr=mat.or.vec(nitem,nitem) #obs correlation
ii_obsyes=mat.or.vec(nitem,nitem) #obs wt resp in 1,1 this case
ii_obsyesprop=mat.or.vec(nitem,nitem) #obs wt prop resp in 1,1
ii_obsyestot=mat.or.vec(nitem,nitem) #obs total wt resp in 1,1
ii_predcorr=mat.or.vec(nitem,nitem) #predicted correlation
ii_predyesprop=mat.or.vec(nitem,nitem) #predicted wt prop resp in 1,1
ii_predyes=mat.or.vec(nitem,nitem) #predicted wt resp in 1,1 this case
ii_predyestot=mat.or.vec(nitem,nitem) #predicted total wt resp in 1,1
ii_residcorr=mat.or.vec(nitem,nitem) #residual correlation

#calculate observed correlations

for (r in 1:nccnext){for (i in 1:nitem){for (ii in 1:nitem){
 ii_obsyestot[i,ii]=
  ii_obsyestot[i,ii]+(ri_in[r,i]>=1)*(ri_in[r,ii]>=1)*r_wt[r]}}}

ii_obsyesprop=ii_obsyestot/totwtcnext

iv_pwt=iv_wt/totwtcnext
iv_pwt=iv_wt/totwtcnext
for (i in 1:nitem){for (ii in 1:nitem){
 ii_obscorr[i,ii]=(ii_obsyesprop[i,ii]-iv_pwt[i,1]*iv_pwt[ii,1])/
 sqrt(iv_pwt[i,1]*(1-iv_pwt[i,1])*iv_pwt[ii,1]*(1-iv_pwt[ii,1]));
 if (i==ii) ii_obscorr[i,ii]=1}}

#calculate predicted correlations

for (p in 1:npat){for (i in 1:nitem){for (ii in 1:nitem){
 ii_predyes[i,ii]=
  (pi_x[p,i]>=1)*(pi_x[p,ii]>=1)*p_probinraw[p]*rr_wt[p_raw[p]]}};
 ii_predyestot=ii_predyestot+ii_predyes}

ii_predyesprop=ii_predyestot/totwtcnext

for (i in 1:nitem){for (ii in 1:nitem){
 ii_predcorr[i,ii]=(ii_predyesprop[i,ii]-iv_pwt[i]*iv_pwt[ii])/
 sqrt(iv_pwt[i,1]*(1-iv_pwt[i,1])*iv_pwt[ii,1]*(1-iv_pwt[ii,1]));
 if (i==ii) ii_predcorr[i,ii]=1}}

for (i in 1:nitem){for (ii in 1:nitem){
 if (i==ii) ii_residcorr[i,ii]=1 else ii_residcorr[i,ii]=
 (ii_obscorr[i,ii]-ii_predcorr[i,ii])/(1-abs(ii_predcorr[i,ii]))}}


##15 write output csv file if write.file=T
#write file header
runtitle=country
ttext1=runtitle
ttext2=paste("\nPartial-credit Rasch model analysis")
# if(fn_csvdata_in!="x") ttext3=paste("\nInput data =,,",fn_csvdata_in) else
#   ttext3=paste("\nInput data directly from R datafile")
ttext=paste(ttext1,ttext2)
fn_output=paste(country,"_Polytomous_Output.csv",sep="")
if(write.file){
write(ttext, file=fn_output, append=F)
}
#write item parameters 
for (i in 1:nitem) for (t in 1:3) if (t>i_max[i]){
  it_rtparam[i,t]=NA; it_error[i,t]=NA; it_param[i,t]=NA;
  it_rtinfit[i,t]=NA; it_rtoutfit[i,t]=NA}
i_df_itemparam=data.frame(i_names,it_rtparam,it_error,it_param)

ttext1=paste("\nItem-threshold parameters")
ttext2=paste("\n N non extreme=,")
ttext3=paste(nccnext)
ttext22=paste("\n N total=,")
ttext33=paste(nrow(XX))
ttext2.1=paste("\n Tot any missing weighted=,")
ttext2.1.1 = paste(nrow(XX)-totwt_anyv)
ttext3.1=paste("\n Perc any missing weighted=,")
ttext3.1.1=paste(100-perc_wt_anyv)
ttext4=paste("\n,Rasch-Thurstone threshold,,,Rasch-Thurstone error,,,Rasch threshold")
ttext5=paste("\nItem,1,2,3,1,2,3,1,2,3")
i_dfheader_itemparam=paste(ttext1,ttext2, ttext3, ttext22, ttext33, ttext2.1, ttext2.1.1, 
                           ttext3.1, ttext3.1.1, ttext4, ttext5)
if(write.file){
write(i_dfheader_itemparam, file=fn_output, append=TRUE)
suppressWarnings(write.table(i_df_itemparam, file=fn_output, append = TRUE, 
  sep = ",", eol = "\n", na = "", dec = ".", row.names = FALSE,
  col.names = FALSE))}

#write item fit stats
i_df_fitstats=data.frame(i_names,it_rtinfit,it_rtoutfit,i_infit,i_outfit)
ttext1=paste("\nItem fit statistics")
ttext2=paste("\n N=,")
ttext3=paste(nccnext)
ttext4=paste("\n,Rasch-Thurstone infit,,,Rasch-Thurstone outfit,,,Overall")
ttext5=paste("\nItem,1,2,3,1,2,3,Infit,Outfit")
i_dfheader_fitstats=paste(ttext1,ttext2,ttext3,ttext4, ttext5)
if(write.file){
write(i_dfheader_fitstats, file=fn_output, append=TRUE)
suppressWarnings(write.table(i_df_fitstats, file=fn_output, append = TRUE, 
  sep = ",", eol = "\n", na = "", dec = ".", row.names = FALSE,
  col.names = FALSE))
}
#write raw score parameters, errors, and iteration information
rr_df_stats=data.frame(
 rr_value,as.numeric(rr_ncase_allv),rr_wt_allv,rr_param,rr_error,rr_niterate,rr_finaladj)
ttext1=paste("\nRespondent statistics by raw score (N and wt no missing)")
ttext2=paste("\nRaw score,Num cases,Wt cases,Parameter,Error,Iterations,Final adj")
rr_dfheader_stats=paste(ttext1,ttext2)
if(write.file){
write(rr_dfheader_stats, file=fn_output, append=TRUE)
suppressWarnings(write.table(rr_df_stats, file=fn_output, append = T, 
  sep = ",", eol = "\n", na = "", dec = ".", row.names = F,
  col.names = F))}
ttext1=paste("\n\nRasch reliability (weighted) = ",raschreliab)
ttext2=paste("\nRasch reliability w equal proportions in all raw scores = ",raschreliabf)
ttext=paste(ttext1,ttext2)
if(write.file){
write(ttext, file=fn_output, append=TRUE)
}
#write resid correlations to csv
suppressWarnings(rm(ii_outcorr))
ii_outcorr=data.frame(ii_residcorr,row.names=NULL)
colnames(ii_outcorr)=i_names
ii_outcorr=data.frame(i_names,ii_outcorr,row.names=NULL)

ttext=paste("\n\nInter-item residual correlations")
if(write.file){
write(ttext, file=fn_output, append=T)
suppressWarnings(write.table(ii_outcorr, file=fn_output, append = TRUE, quote = FALSE, 
  sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F,
  col.names = T))}

#write distribution by raw score if any valid and all valid
# rrall_df_rawdist=data.frame(rr_value,rr_ncase_anyv,rr_wt_anyv,rr_pctwt_anyv,
#   rr_ncase_allv,rr_wt_allv,rr_pctwt_allv)
rrall_df_rawdist=data.frame(rr_value,as.numeric(rr_ncase_anyv),rr_wt_anyv,rr_pctwt_anyv,
                            as.numeric(rr_ncase_allv),rr_wt_allv,rr_pctwt_allv)
ttext1=paste("\nDistribution of respondents by raw score - any valid and all valid")
ttext2=paste("\n,Any valid item response,,,All items valid response")
ttext3=paste("\nRaw score,Num case,Wt case,Pct wt,Num case, Wt case,Pct wt")
rrall_dfheader_rawdist=paste(ttext1,ttext2,ttext3)
if(write.file){
write(rrall_dfheader_rawdist, file=fn_output, append=T)
suppressWarnings(write.table(rrall_df_rawdist, file=fn_output, append = T, 
  sep = ",", eol = "\n", na = "", dec = ".", row.names = F,
  col.names = F))}

#write observed prop by category, cnext
tobs1=ri_in[1:rawmax-1,] #to get item names at column heads
for (rr in 1:rawmax-1) for (i in 1:nitem) tobs1[rr,i]=rriv_obspropinv[rr,i,1]
tobs2=ri_in[1:rawmax-1,] #to get item names at column heads
for (rr in 1:rawmax-1) for (i in 1:nitem) {
  if (i_max[i]<2) tobs2[rr,i]=NA 
  else tobs2[rr,i]=rriv_obspropinv[rr,i,2]}
tobs3=ri_in[1:rawmax-1,] #to get item names at column heads
for (rr in 1:rawmax-1) for (i in 1:nitem) {
  if (i_max[i]<3) tobs3[rr,i]=NA 
  else tobs3[rr,i]=rriv_obspropinv[rr,i,3]}

rr_value_cnext=seq(1:(rawmax-1))
Raw_score=rr_value_cnext #labels will be used for columns
rr_df_obsprop=data.frame(Raw_score,tobs1,tobs2,tobs3)

rr_dfheader_obsprop=paste("\nObserved item-by-category proportion of wt responses")
if(write.file){
write(rr_dfheader_obsprop, file=fn_output, append=TRUE)
suppressWarnings(write.table(rr_df_obsprop, file=fn_output, append = TRUE, 
  sep = ",", eol = "\n", na = "", dec = ".", row.names = F,
  col.names = T))
}
#write predicted prop by category, cnext

tpred1=ri_in[1:rawmax-1,] #to get item names at column heads
for (rr in 1:rawmax-1) for (i in 1:nitem) tpred1[rr,i]=rriv_predpropinv[rr,i,1]
tpred2=ri_in[1:rawmax-1,] #to get item names at column heads
for (rr in 1:rawmax-1) for (i in 1:nitem) {
  if (i_max[i]<2) tpred2[rr,i]=NA 
  else tpred2[rr,i]=rriv_predpropinv[rr,i,2]}
tpred3=ri_in[1:rawmax-1,] #to get item names at column heads
for (rr in 1:rawmax-1) for (i in 1:nitem) {
  if (i_max[i]<3) tpred3[rr,i]=NA 
  else tpred3[rr,i]=rriv_predpropinv[rr,i,3]}

Raw_score=rr_value_cnext #labels will be used for columns
rr_df_predprop=data.frame(Raw_score,tpred1,tpred2,tpred3)

rr_dfheader_predprop=paste("\nPredicted item-by-category proportion of wt responses")
if(write.file){
write(rr_dfheader_predprop, file=fn_output, append=TRUE)
suppressWarnings(write.table(rr_df_predprop, file=fn_output, append = TRUE, 
  sep = ",", eol = "\n", na = "", dec = ".", row.names = F,
  col.names = T))}

#write missing cases analysis
vall_value=0:nitem
vall_df_missing=data.frame(vall_value,vall_ncaseinv,vall_pctncaseinv,vall_wtcaseinv,
 vall_pctwtcaseinv,vall_pctncaseinv_anyv,vall_pctwtcaseinv_anyv)
ttext1=paste("\nDistribution of valid responses")
ttext2=paste("\nNum valid,Num cases,Pct cases,Wt cases,Wt pct,Pct if any valid,Wt pct if any valid")
vall_dfheader_missing=paste(ttext1,ttext2)
if(write.file){
write(vall_dfheader_missing,file=fn_output, append=T)
suppressWarnings(write.table(vall_df_missing, file=fn_output, append = TRUE, 
  sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
  col.names = FALSE))
}
#write missing by item if any valid
i_df_missing=data.frame(i_names,i_nmiss,i_pctnmiss,i_wtmiss,i_pctwtmiss)
i_dfheader_missing=paste("\nMissing by item if any valid\nItem,Num missing,Pct missing, Wt missing, Wt pct missing")
if(write.file){
write(i_dfheader_missing,file=fn_output,append=T)
suppressWarnings(write.table(i_df_missing, file=fn_output, append = TRUE, 
  sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,
  col.names = FALSE))
}
#write response by category, all cases
for (i in 1:nitem) for (c in 1:4) if (c>i_max[i]+1){
  ic_numcase_all[i,c]=NA;ic_wtcase_all[i,c]=NA; ic_pctwt_all[i,c]=NA}

i_df_resp=data.frame(i_names,i_nvalid,ic_numcase_all,ic_wtcase_all,ic_pctwt_all)
ttext1=paste("\nResponses in each category - cases with valid response to the item")
ttext2=paste("\n,,Number of cases,,,,Wt cases,,,,Pct of wt cases")
ttext3=paste("\nItem,Num valid,0,1,2,3,0,1,2,3,0,1,2,3")
i_dfheader_resp=paste(ttext1,ttext2,ttext3)
if(write.file){
write(i_dfheader_resp, file=fn_output, append=TRUE)
suppressWarnings(write.table(i_df_resp, file=fn_output, append = TRUE, 
  sep = ",", eol = "\n", na = "", dec = ".", row.names = FALSE,
  col.names = FALSE))}

#write response by category, complete nonextreme cases
for (i in 1:nitem) for (c in 1:4) if (c>i_max[i]+1) ic_numcase_cnext[i,c]=NA
 
i_df_resp_cnext=data.frame(i_names,ic_numcase_cnext)
ttext1=paste("\nResponses in each category - complete nonextreme responses")
ttext2=paste("\n N=,")
ttext3=paste(nccnext)
ttext4=paste("\nItem,0,1,2,3")
i_dfheader_resp_cnext=paste(ttext1,ttext2,ttext3,ttext4)

if(write.file){
write(i_dfheader_resp_cnext, file=fn_output, append=TRUE)
suppressWarnings(write.table(i_df_resp_cnext, file=fn_output, append = TRUE, 
  sep = ",", eol = "\n", na = "", dec = ".", row.names = FALSE,
  col.names = FALSE))}

#write specifications and analysis information
t1=paste("\nSpecifications and Information")
# if(fn_csvdata_in!="x") t2=paste("\n\nInput data =,,",fn_csvdata_in) else
#   t2=paste("\n\nInput data directly from R datafile")
t3=paste("\nnitem,,Number of items = ",nitem)
# t4=paste("\nloc_item1,,Location of first item = ",loc_item1)
wtanalysis=ifelse(sum(wt==1)==length(wt), FALSE, TRUE)
t5=paste("\nwtanalysis,,Use weights for analysis? = ",wtanalysis)
# if (wtanalysis==TRUE) 
#   t6=paste("\nloc_wt,,Location of weight variable = ",loc_wt) else t6=paste(" ")

fn_item_recode=if(recode==0){"No recoding"}else
  if(recode==1){"Almost every month + In some months but not every month"}else
    if(recode==2){"In some months but not every month + Only once or twice"}
t7=paste("\nfn_item_recode,,Item recoding used = ",fn_item_recode)
t8=paste("\nfn_output,,Output datafile = ",fn_output)
t9=paste("\nmaxiter,,Maximum iterations for item parameter estimation = ",maxiter)
t10=paste("\nminconv,,Convergence criterion = ",minconv)
t11=paste("\nMain estimation iterations = ",niterate)
t12=paste("\nLargest parameter adjustement at last iteration = ",maxadj)
t13=paste("\nSmallest log-probability of a modeled pattern = ",minlogpprob)
t14=paste("\nProportion of patterns with log-probability < -34.5 = ",propnotcalc)

# specinfo=paste(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14)
specinfo=paste(t1,t3,t5,t7,t8,t9,t10,t11,t12,t13,t14)
if(write.file){
write(specinfo, file=fn_output, append=TRUE)
}
#write recode matrix
# i_df_recode=data.frame(i_names,iv_recode)
# ttext1=paste("\n\nItem recode matrix (tabled value=rcode with -1=missing)")
# ttext2=paste("\n\n,Original code in input data")
# ttext3=paste("\nItem,0,1,2,3,4,5")

# i_dfheader_recode=paste(ttext1,ttext2,ttext3)
# write(i_dfheader_recode, file=fn_output, append=TRUE)
# suppressWarnings(write.table(i_df_recode, file=fn_output, append = T, 
#   sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F,
#   col.names = F))
# 
# tversion=paste("\n\n",mrpcm_version)
# write(tversion, file=fn_output, append=TRUE)

##15 create list for summary function output
t1="Rasch-Thurstone parameters"
i_df_rtparam=data.frame(i_names,it_rtparam)
t1.5=paste("N complete nonextreme = ",nccnext)
t1.6=paste("Number of iterations = ",niterate)
t1.7=paste("Largest adjustment in final iteration = ",maxadj)
t2="Parameter errors"
i_df_error=data.frame(i_names,it_error)
t3="Rasch-Thurstone infit"
i_df_rtinfit=data.frame(i_names,it_rtinfit)
t4="Rasch-Thurstone outfit"
i_df_rtoutfit=data.frame(i_names,it_rtoutfit)
t5="Overall infit and outfit"
i_df_overallfit=data.frame(i_names,i_infit,i_outfit)
raschreliabr=round(raschreliab,3)
raschreliabfr=round(raschreliabf,3)
t5.5=paste("Rasch reliability = ",raschreliabr)
t6=paste("Reliability flat = ",raschreliabfr)
t7=paste("Weighted distribution of cases w any valid response")
#rrall_df_rawdist
t8=paste("Inter-item residual correlations")
#ii_outcorr

mrpcm_out=list("country"=runtitle,"b"=i_df_rtparam,"se.b"=i_df_error,
#                "nit"=niterate,
#                "lastadj"=maxadj,
               "infit"=i_df_rtinfit,
               "outfit"=i_df_rtoutfit,"reliab"=raschreliab,"reliab.fl"=raschreliabfr,
               "RS_distr"=rrall_df_rawdist,"res.cor"=ii_outcorr,"a"= rr_param,
               "se.a"=rr_error,XX=XX,wt=wt, d=c(extr_lo,extr_hi),
               "N_valid_w"=totwt_anyv, "Perc_valid_w"=perc_wt_anyv,
               "N_tot"=nrow(XX), "N_compl"=nccnext
               )
return(mrpcm_out)

} #end of function
#end

# Try function on trial data
# rr.pol=PC.w(XX,wt,recode=0,country="Prova2",write.file=F)
