\name{RM.w}
\alias{RM.w}
\title{Estimation of the weighted Rasch model using CML}
\description{
This function computes item and raw score parameter estimates of a Rasch model for binary item responses by using weighted CML estimation. Input data should be a 0/1 matrix (1 = Yes). Residual correlation, fit statistics and corresponding standard errors, Rasch reliability and individual fit statistics are also reported.
}
\usage{
RM.w(.data, .w = NULL, .d=NULL, country=NULL, se.control = TRUE, 
quantile.seq = NULL, write.file = FALSE, max.it=100)
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
      \item{max.it}{Maximum number of iterations item parameters' estimation.}
      }
\value{
  A list with the following elements: 
  \tabular{ll}{
	\code{country} \tab Name of the country or the application. \cr
	\code{b} \tab Item severity parameters. \cr 
	\code{a} \tab Raw score severity parameters. \cr 
\code{se.b} \tab Item severities' standard error. \cr 
\code{se.a} \tab Raw score parameters' standard error. \cr
\code{infit} \tab Items' infit statistics. \cr 
\code{outfit} \tab Items' outfit statistics. \cr 
\code{reliab} \tab Rasch reliability. \cr
\code{reliab.fl} \tab Rasch reliability (equally weighted across raw scores). \cr
\code{infit.person} \tab Respondents' infits. \cr
\code{outfit.person} \tab Respondents' outfits. \cr
\code{q.infit.theor} \tab Quintiles of theoretical person infits. \cr
\code{q.infit} \tab Quintiles of person infits. \cr
\code{q.outfit.theor} \tab Quintiles of theoretical person outfits. \cr
\code{q.outfit} \tab Quintiles of person outfits. \cr
\code{res.corr} \tab Residual correlation matrix. \cr
\code{se.infit} \tab Standard errors of infits. \cr
\code{mat.res} \tab Matrix of individual residuals for each item. \cr
\code{d} \tab Pseudo raw scores for the calculation of extreme raw score parameters. \cr
\code{XX} \tab Data matrix. \cr
\code{wt} \tab Vector of post-stratification sampling weights. \cr
\code{n.compl} \tab Number of complete, non-extreme cases. \cr
\code{wt.rs} \tab Weighted (absolute) distribution of respondents across raw scores. \cr
\code{wt.rs} \tab Weighted (relative) distribution of respondents across raw scores. \cr
\code{converged} \tab Character indicating if the item estimation algorithm has reached convergence
or the maximum number of iterations has been reached without convergence. \cr
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

# Checking conditional independence: residual correlation matrix 
rr.country1$res.cor

}


