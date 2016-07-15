\name{equating.fun}
\alias{equating.fun}
\title{Equating scales between different application contexts.}
\description{
This function calibrates the measure derived by a scale applied in one context (for example, country) to the metric of a reference scale, or standard (for example, the scale applied in another country, or in the same country but through a different survey, or to any other standard).
The main output is the prevalence rate in the country of interest calculated at specified thresholds along the latent trait. 
}
\usage{
equating.fun(rr1, st=NULL, tol = .35, spec.com1 = 1:8, spec.com2=1:8, 
thres = c(-0.25,  1.83), maxuniq=3, plot=F,iterative=T,excl.prior1,
excl.prior2)
}
\arguments{
  \item{rr1}{The Rasch model fitted with the \code{RM.w} function for the country of interest.}
  \item{st}{Item severity parameter estimates for the reference country or the standard. The number of items for the standard can differ from the number of items of the country of interest.}
  \item{tol}{Tolerance expressed in the metric of the standard, above which items are considered unique (i.e. not contributing to the common metric between the two countries). Default is \code{0.35}.}
  \item{spec.com1}{A priori numbered set of comparable items for the country of iterest to the standard. Items must be specified in the same order as they are used for the comparison.}
  \item{spec.com2}{A priori set of comparable items for the standard. The length of \code{spec.com1} and \code{spec.com2} should be the same.}
  \item{thres}{Thresholds (along the latent trait) corresponding to which the comparable prevalence rate of the phenomenon of interest is calculated.}
  \item{maxuniq}{A priori maximum number of unique items allowed.}
    \item{plot}{Logical argument. If \code{TRUE}, a pdf file with a plot of the equated item severities is produced in the working directory. Default is \code{FALSE}.}
      \item{iterative}{Logical argument. If \code{TRUE}, an iterative process to identify which item is unique is performed. Otherwise, unique items are specified in the arguments \code{excl.prior1} and \code{excl.prior2} and are taken as fixed.  Default is \code{TRUE}. }
      \item{excl.prior1}{Numbered set of fixed unique items for the country of interest. It needs to be specified only if \code{iterative = FALSE}.}
      \item{excl.prior2}{Numbered set of fixed unique items for the standard. It needs to be specified only if \code{iterative = FALSE}.}
}
\value{
\item{scale}{Scale parameter to be applied to the metric of the country of interest to be considered equivalent to the standard. }
\item{shift}{Shift parameter to be applied to the metric of the country of interest to be considered equivalent to the reference country. The standardized item parameters for the country of interest will be:
b.country.st=shift+scale*b.country
where b.country are estimated item severities for the country of interest.
}
\item{common}{Vector of logical, \code{TRUE} when the item in the country of interest is common.}
\item{prevs}{Prevalence rates for the country of interest, calculated on the metric of the standard.}
\item{prevs.rs}{Prevalence rates for the country of interest, calculated on the metric of the standard, calculated at each raw score threshold.}
\item{cor.comm.items}{Correlation between common items.}
}

\seealso{prob.assign}

\author{Sara Viviani \email{sara.viviani@fao.org}}
\examples{
\dontrun{
data(data.FAO_country1)
XX.country1 = data.FAO_country1[,1:8]
wt.country1 = data.FAO_country1$wt
rr1=RM.w(XX.country1, wt.country1)
data(data.FAO_country2)
XX.country2 = data.FAO_country2[,1:8]
wt.country2 = data.FAO_country2$wt
rr2=RM.w(XX.country2, wt.country2)
# Iterative equating
ee1=equating.fun(rr1, st=rr2$b, thres=c(-0.28, 1.85), spec.com1=c(1:8),spec.com2=c(1:8))
# Not iterative equating: fixing WORRIED unique
ee2=equating.fun(rr1, st=rr2$b, thres=c(-0.28, 1.85),spec.com1=c(1:8),
spec.com2=c(1:8), tol = .4,iterative=F,excl.prior1=1,excl.prior2=1)
# Prevalence rates
ee1$prevs
ee2$prevs
}
}


