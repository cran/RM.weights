\name{ICC.fun}
\alias{ICC.fun}
\title{Item caracteristic curves computation and plotting under CML.}
\description{
Computation and plotting of item caracteristic curves (ICCs) under conditional maximum likelihood (CML).
}
\usage{
ICC.fun(b, plot = F)
}
\arguments{
  \item{b}{Item severity vector estimated under CML.}
  \item{plot}{Logical. If \code{TRUE}, the ICCs are plotted.}
}
\value{
	\tabular{ll}{
	\code{icc} \tab Matrix of item characteristic curves. Coloumns are items
  and rows are severity levels. \cr
	}
}
\author{Sara Viviani \email{sara.viviani@fao.org}}
\seealso{
  \code{\link{RM.w}}
}
\examples{
\dontrun{
data(data.FAO_country2)
# Questionnaire data and weights
XX.country2 = data.FAO_country2[,1:8]
wt.country2 = data.FAO_country2$wt
# Fit weighted Rasch
rr.country2 = RM.w(XX.country2, wt.country2)
# ICC values
icc = ICC.fun(b)$icc
# ICC plot
icc = ICC.fun(b, plot = T)
}
}


