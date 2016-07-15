\name{PC.w}
\alias{PC.w}
\title{Estimation of the weighted partial credit model using CML}
\description{
This function estimates item and raw score parameters of a partial credit Rasch model which may include both dichotomous and polytomous item responses by using weighted CML estimation. The maximum number of response categories allowed is 4. Input data should be a 0/1 matrix (1 = Yes) for the dichotomous items and 0/1/2 or 0/1/2/3 for the polytomous items. Residual correlation, fit statistics and Rasch reliability are also reported.
}
\usage{
PC.w(XX, wt=NULL, extr=NULL, maxiter=100,minconv=.00001,country=NULL,
write.file=F, recode = 0, write.iteration=F)
}
\arguments{
  \item{XX}{Input data matrix or data frame. Rows represent respondents, columns represent items. Dichotomous items are coded 0/1, polytomous items are coded 0/1/2 or 0/1/2/3; negative responses must be coded as 0s. Missing values must be coded as \code{NA}.}
  \item{wt}{Vector of sampling weights. 
  The length must be the same as the number of rows of \code{XX}. 
  If left unspecified, 
  all the individuals will be equally weighted (\code{wt = rep(1, nrow(XX)))}. }
  \item{extr}{Optional vector for the assumption on the extreme raw score parameters. Default is \code{0.5} and \code{(k-0.5)},  \code{k} being the maximum number of response categories.} 
    \item{maxiter}{Maximum number of iterations for convergence when parameter estimation is performed.}
    \item{minconv}{Minimum difference between parameter estimates considering two consecutive iterations for convergence.}
    \item{country}{Optional (character) name of the dataset.}
      \item{write.file}{If \code{TRUE}, a CSV file with the main model results will be saved in the working directory.}
      \item{recode}{It can be 0, 1, 2 or 3. \code{recode = 0} does not imply any category aggregation. \code{recode = 1} aggregates the two more severe categories (for example, Sometimes and Always). \code{recode=2}  aggregates the two middle categories (for example, Rarely and Sometimes). \code{recode=3}  aggregates the two less severe categories (for example, Never and Rarely). Default is \code{0}.}
\item{write.iteration}{If set \code{TRUE} a CSV file with information on the iteration process is produced in the working directory.}
      }
\value{
  A list with the following elements: 
  \tabular{ll}{
	\code{b} \tab Matrix of Rasch-Thurstone (R-T) thresholds (corresponding to item severity parameters in a Rasch model for dichotomous data). \cr
	\code{se.b} \tab Matrix of standard errors for R-T thresholds. \cr 
  \code{infit} \tab Infit statistics corresponding to each R-T threshold. \cr 
  \code{outfit} \tab Outfit statistics corresponding to each R-T threshold. \cr
  \code{reliab} \tab Rasch reliability. \cr
  \code{reliab.fl} \tab Rasch reliability (equally weighted across raw scores). To be used for comparisons. \cr
  \code{RS_distr} \tab Distribution of cases across raw scores (weighted and unweighted). \cr
  \code{res.cor} \tab Matrix of residual correlation. \cr
  \code{a} \tab Vector of raw score parameters. \cr
  \code{se.a} \tab Vector of raw score parameter errors (or measurement errors). \cr
  \code{XX} \tab Data matrix (after recoding). \cr
  \code{wt} \tab Vector of sampling weights. \cr
  \code{d} \tab Pseudo-raw-scores used for the estimation of extreme parameters. \cr
  \code{N_valid_w} \tab Number of cases with any valid response.  \cr 
  \code{Perc_valid_w} \tab Percentage of cases with any valid response. \cr 
  \code{N_tot} \tab Total sample size. \cr 
  \code{N_compl} \tab Number of complete non-extreme cases. \cr 
	}
}      
      
\author{Sara Viviani \email{sara.viviani@fao.org}, Mark Nord \email{mark.nord@fao.org}}
\details{
The weighted CML method is used to estimate the item and Rasch-Thurstone parameters. Respondent parameters  are estimated post-hoc as maximum likelihood given the item parameters. Cases with missing responses to some items can be included, but will not be used to estimate the Rasch model. Missing responses should be coded as \code{NA} in the input data file. If the total possible raw score is higher than 15, computer precision limits may distort results.

As the parameters for the extreme raw scores (\code{0} and \code{k}), are undefined under the CML, some assumptions are needed unless the proportions of respondents with those raw scores are so small that they can be considered to be measured without error. 
Vector \code{extr}  gives the option to include one assumption on the pseudo-raw-score on which raw score zero will be based and one assumption on the pseudo-raw-score on which raw score \code{k} will be based. The assumption for raw score zero can be a number between 0 and 1 (extremes excluded). The assumption for raw score \code{k} can be a number between \code{k-1} and \code{k} (extremes excluded).
}
\examples{
\dontrun{
data(data.FAO_country1)
# Questionnaire data and weights
XX.country1 = data.FAO_country1[,1:8]
XX.country1[,c(7,8)] = cbind(rbinom(nrow(XX.country1), 3, .5), rbinom(nrow(XX.country1), 3, .4)) 
wt.country1 = data.FAO_country1$wt
# Fit weighted partial credit model with aggregation of the two more severe categories
rr.country1.pol2 = PC.w(XX.country1, wt.country1, recode = 1)
# Write output file
rr.country1.pol2 = PC.w(XX.country1, wt.country1, recode = 1, write.file=T)
}
}

