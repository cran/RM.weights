\name{prob.assign}
\alias{prob.assign}
\title{Probabilistic assignment of cases to classes of severity along the latent trait.}
\description{
This function assigns cases probabilistically to classes of severity along the latent trait. This procedure is useful when different (cultural-geographical-linguistic) contexts are compared in terms of the prevalence of some phenomenon. 
}
\usage{
prob.assign(rr=NULL,  rwthres = NULL, sthres = NULL, 
eps.a = NULL, flex = list(a=NULL, se.a=NULL, d=NULL, XX=NULL, wt=NULL))
}
\arguments{
  \item{rr}{An object of \code{RM.w} type.}
  \item{rwthres}{Raw score thresholds. If this argument is specified, corresponding thresholds on the latent trait (\code{thres}) are provided as output.}
  \item{sthres}{Alternative argument to \code{rwthres}. Thresholds on the latent trait. If this argument is specified, probability of being beyond the specified thresholds (\code{sprob}) is provided as output.}
  \item{eps.a}{Tolerance for the algorithm that estimates thresholds on the latent trait corresponding to the specified thresholds in terms of raw score (\code{rwthres}). If left unspecified, the tolerance to convergence is \code{0.001}.}
  \item{flex}{Alternative argument to \code{rr}. If \code{rr=NULL}, than \code{flex} needs
  to be specified as a list with the following arguments: \code{a} (person parameters), \code{se.a} (person parameter standard errors),
  \code{d} (pseudo-raw-scores for estimation of extreme raw score parameters; if left unspecified leads to estimate the \code{k}-th raw score parameter at location \code{(k-1)+0.5}), \code{XX} (input data matrix or data frame), \code{wt} (sampling weights; if left unspecified leads to sampling weights equal to one for all subjects).}
}
\details{
The probabilistic assignment procedure is particularly useful when comparing results between different populations that might have interpreted some items of a scale differently. The distribution of the raw score is assumed to be a mixture of Gaussian densities, each centred on the raw score parameters (estimated using the Rasch model) and scaled with corresponding measurement errors. The resulting (complementary) cumulative prevalence, weighted by the proportion of individuals in each raw score, is as follows:
\deqn{
P(x)= \sum_{i=0}^N (1-\int_{-Inf}^x f_i(t)dt) p_i = 
\sum_{i=0}^N (1-F_i(x)) p_i,
}
where \eqn{f_i} is the probability density function (and \eqn{F_i} is the cumulative distribution function) of severity levels centred on raw score parameter \eqn{i} and scaled with the corresponding measurement error, \eqn{p_i} is the proportion observed in raw score \eqn{i}, and \eqn{N} is the maximum number of items. The function \eqn{P(x)} can be used to trace a continuous profile of prevalence on the latent trait for the phenomenon of interest. It allows to equate different metrics and it facilitates cross-cultural comparisons.
The \code{prob.assign} function can be used in two ways: providing the \code{rwthres} argument, i.e. specifying thresholds in terms of raw score, it provides the estimation of the thresholds on the latent trait that correspond to the raw score thresholds;  providing the \code{sthres} argument, i.e. specifying thresholds on the latent trait, it provides the probability of being beyond the specified thresholds for the population of interest.
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{sprob} \tab Estimated weighted probability of being beyond thresholds provided in \code{sthres} (\eqn{P(x)}). If \code{sthres} argument is not specified, the \code{sprob} value is \code{NULL}. \cr
	\code{thres} \tab Thresholds on the latent trait calculated corresponding to thresholds in terms of raw-score specified in \code{rwthres}. \cr 
	\code{f} \tab Probability of being beyond the \code{thres} thresholds distributed across raw scores. If more than one assumption on the exreme raw scores is made in \code{rr}, it is a list of a number of elements equal to the number of assumptions made. \cr 
	\code{p} \tab Empirical (weighted) distribution beyond the raw scores specified in \code{rwthres}. If only one assumption on the extreme raw scores is made, then 
  \code{colSums(f)} is approximately equal to \code{p}, where the order of approximation depends on 
  \code{eps.a}.\cr 
  \code{f_j} \tab Empirical (weighted) distribution across the raw scores. \cr 
	}
}

\author{Sara Viviani \email{sara.viviani@fao.org}, Mark Nord \email{mark.nord@fao.org}}
\seealso{
  \code{\link{RM.w}}
}
\examples{
\dontrun{
data(data.FAO_country1)
# Questionnaire data and weights
XX.country1 = data.FAO_country1[,1:8]
wt.country1 = data.FAO_country1$wt
# Fit weighted Rasch
rr.country1 = RM.w(XX.country1, wt.country1)

# Thresholds on the latent trait corresponding to a given raw score
pp.country1 = prob.assign(rr.country1, rwthres = c(3, 7))

# Table with prevalences and thresholds
tab = cbind("Raw score" = c(3, 7), "Latent trait" =  pp.country1$thres, 
"Prevalence" = colSums(pp.country1$f))
rownames(tab) = c("Thres 1","Thres 2")
tab

# Pre-defined thresholds on the latent trait
sthresh = c(-0.25, 1.81)
pp.country1.2 = prob.assign(rr.country1, sthres = sthresh)$sprob
# Probability of being beyond -0.25 on the latent trait in country 1 
pp.country1.2[1]*100
# Probability of being beyond 1.81 on the latent trait in country 1 
pp.country1.2[2]*100

# More than 2 extremes
# Fit the model
rr.country1.d = RM.w(XX.country1, wt.country1, .d = c(0.5, 7.5, 7.7))
# Probabilistic assignment
pp.country1.d = prob.assign(rr.country1.d, sthres = sthresh)$sprob
# Probability of being beyond -0.25 on the latent trait in country 1 
# using upper assumption on the extreme raw score parameter 8
pp.country1.d[[1]]*100
# Probability of being beyond -0.25 on the latent trait in country 1 
# using lower assumption on the extreme raw score parameter 8
pp.country1.d[[2]]*100

}

}
