EWaldtest <-
function(b1, b2, se1, se2) {
    z = (b1 - b2)/(sqrt(abs(se1^2 + se2^2)))
    p.value = 2 * pnorm(-abs(z))
    tab = cbind(z, p.value)
    rownames(tab) = 
  return(list(z = z, p = p.value, tab = round(tab, 3 )))
}
