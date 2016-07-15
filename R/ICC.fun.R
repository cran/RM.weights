ICC.fun <-
function(b, plot = F){
  P.i <- function(x, .beta = 1){
    teil2a <- x - .beta 
    teil2 <- (1+exp(teil2a))
    teil3 <- exp(teil2a)
    LL <- teil3/teil2  
    return(LL)
  }
  x <- seq(-5, 5, 0.01)
  icc = sapply(1:length(b), function(i)  P.i(x, .beta = b[i]))
  if(plot) {
    plot(x, icc[,1], type = "l", col = 1, xlab = "Latent Trait",
         ylab = "Probability of affirmative response", 
         main = "Item Characteristic Curves")
    for(i in 2:ncol(icc)){
      lines(x, icc[,i], type = "l", col = i)
    }
  }
  return(list(icc = icc))
}
