\name{RM.w}
\alias{RM.w}
\title{Estimation of the weighted Rasch model using CML}
\description{
This function computes item and raw score parameter estimates of a Rasch model for binary item responses by using weighted CML estimation. Input data should be a 0/1 matrix (1 = Yes). Residual correlation, fit statistics and corresponding standard errors, Rasch reliability and individual fit statistics are also reported.
}
\usage{
RM.w(.data, .w = NULL, .d=NULL, country=NULL, se.control = T, 
quantile.seq = NULL, write.file = F)
}
\arguments{
  \item{.data}{Input 0/1 data matrix or data frame; 
  affirmative responses must be coded as 1s. Rows represent respondents, columns represent items. Missing values must be coded as \code{NA}.}
  \item{.w}{Vector of sampling weights. 
  The length must be the same as the number of rows of \code{.data}. 
  If left unspecified, 
  all the individuals will be equally weighted (\code{.w = rep(1, nrow(.data)))}. }
  \item{.d}{Optional vector for the assumption on the extreme raw score parameters. 
  Default is \code{0.5} and \code{(k-0.5)},  \code{k} being the maximum number of items 
  (columns of \code{.data}).} 
    \item{country}{Optional (character) name of the dataset.}
    \item{se.control}{Are the measurement errors for respondent extreme parameters to be calculated based on raw scores \code{0.5} and \code{(k-0.5)}? If \code{FALSE}, measurement errors for those parameters are based on the pseudo-raw-scores specified in the \code{.d} vector. Default is \code{TRUE}.}
    \item{quantile.seq}{Quantiles corresponding to the observed and the expected 
    individual fit statistic distributions. Default is \code{seq(0,1,0.01)}.}
      \item{write.file}{If \code{TRUE}, a CSV file with the main results will be saved 
      in the working directory.}
      }
\value{
# TBD #
  A list with the following elements: 
  \tabular{ll}{
	\code{b1.rt} \tab R-T threshold for the first category. \cr
	\code{b2.rt} \tab R-T threshold for the second category. \cr 
	\code{b12.rt} \tab Vector of R-T thresholds for the first and second category. \cr 
	}
}
      
\author{Sara Viviani \email{sara.viviani@fao.org}, Mark Nord \email{mark.nord@fao.org}}
\details{
The weighted CML method is used to estimate the item parameters. Respondent parameters  are estimated post-hoc as maximum likelihood given the item parameters. Cases with missing responses to some items can be included in the input data, but will not be used to estimate the Rasch model. Missing responses should be coded as \code{NA} in the input data file.

As the parameters for the extreme raw scores (\code{0} and \code{k}), are undefined under the CML, some assumptions are needed unless the proportions of respondents with those raw scores are so small that they can be considered to be measured without error. 
Vector \code{.d}  gives the option to include up to four alternative assumptions on the pseudo-raw-scores on which calculation of parameters for extreme raw scores will be based. More in detail, \code{.d} can be a two, three or four element vector:
\itemize{
\item If \code{length(.d) = 4}, then the first two elements must refer to the assumptions for raw score 0, and the second two elements for raw score \code{k}. For instance \code{.d = c(0.1, 0.7, 7.1, 7.6)}, if the maximum raw score is 8.
\item If \code{length(.d) = 3}, then the first two elements can either refer to the assumptions upon raw score \code{0} or raw score \code{k}, and the last one is defined accordingly. For instance \code{.d = c(0.1, 7.1, 7.6)} or \code{.d = c(0.1, 0.7, 7.6)}, if the maximum raw score is 8.
\item If \code{length(.d) = 2}, then the first element must refer to the assumption for raw score \code{0}, and the second element for raw score \code{k}. For instance \code{.d = c(0.1, 7.6)}, if the maximum raw score is 8.
}
}
\examples{
\dontrun{
data(data.FAO_country1)
# Questionnaire data and weights
XX.country1 = data.FAO_country1[,1:8]
wt.country1 = data.FAO_country1$wt
# Fit weighted Rasch
rr.country1 = RM.w(XX.country1, wt.country1)
# Fit unweighted Rasch
rr.country1.nw = RM.w(XX.country1)

# Item severity
rr.country1$b

# Item standard error
rr.country1$se.b

# Respondent severity for each raw score
rr.country1$a

# Respondent measurement error for each raw score
rr.country1$se.a

# Item infit  
rr.country1$infit

# Item outfit 
rr.country1$outfit

# Rasch reliability based on observed distribution of cases across raw scores
rr.country1$reliab

# Rasch reliability  based on equal proportion of cases in each 
# non-extreme raw score (more comparable across datasets)
rr.country1$reliab.fl

# Respondent infit distribution: observed and expected
quantile.seq = c(0,.01,.02,.05,.10,.25,.50,.75,.90,.95,
   .98,.99,1)
q.infit = rr.country1$q.infit
q.infit.theor  = rr.country1$q.infit.theor
plot(quantile.seq, q.infit, type = "b", xlab = "Quantiles", 
ylab = "Observed infit", ylim = c(0, 6))
lines(quantile.seq, q.infit.theor, type = "b", col = 2)

# Checking conditional independence: residual correlation matrix 
rr.country1$res.cor

# Save outputs to csv file with data name
rr.country1 = RM.w(XX.country1, wt.country1, country = "country1", write.file = T)
}}


