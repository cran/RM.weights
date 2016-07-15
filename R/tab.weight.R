tab.weight <-
function(variab = NULL, wt, XX = NULL){
  if(!is.null(XX)){
    rv = rowSums(XX)
    k = ncol(XX)
    tab1 = sapply(1:(k+1), function(i) sum(wt[rv==i-1],na.rm=T))
    tab2 = sapply(1:(k+1), function(i) sum(wt[rv==i-1],na.rm=T)/sum(wt,na.rm=T))
    tab3 = table(rv)
    tab4 = tab3/sum(tab3)
    tab5 = apply(XX, 2, sum, na.rm = T)/nrow(XX)
    tab6 = sapply(1:k, function(i) sum(wt[XX[,i] == 1], na.rm = T)/
                    sum(wt, na.rm = T))
    tab1 = round(tab1, 2)
    tab2 = round(tab2, 2)
    tab3 = round(tab3, 2)
    tab4 = round(tab4, 2)
    tab5 = round(tab5, 2)
    tab6 = round(tab6, 2)
    names(tab5) = names(tab6) = colnames(XX)
  } else {
    tab1 = tab2 = tab3 =  tab4 = tab5 = tab6 = rv = NULL
  }
  if(!is.null(variab)){
    if(is.list(variab)){
      var1 = variab[[1]]
      var2 = variab[[2]]
      lev1 = levels(var1)
      lev2 = levels(var2)
      tab.ext.w = sapply(1:length(lev1), function(i)
        sapply(1:length(lev2), function(j)
          sum(wt[(var1)==lev1[i] & (lev2)==lev2[j]], na.rm=T)
        ) )
      rownames(tab.ext.w) = lev2
      colnames(tab.ext.w) = lev1
      tab.ext.w = round(tab.ext.w, 2)
    } else {
      var1 = variab
      lev1 = levels(var1)
      tab.ext.w = sapply(1:length(lev1), function(i)
        sum(wt[(var1)==lev1[i]], na.rm=T) )
      names(tab.ext.w) = lev1
      tab.ext.w = round(tab.ext.w, 2)
    }
  } else tab.ext.w = NULL  
  return(list("RS.abs.w" = tab1, "RS.rel.w" = tab2, "RS.abs" = tab3, "RS.rel" = tab4, 
              "Perc.Yes" = tab5, "Perc.Yes.w" = tab6, rv = rv, tab.ext.w=tab.ext.w
  ))
}
