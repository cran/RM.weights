\name{EWaldtest}
\alias{EWaldtest}

\title{Wald test to check sampling independence under CML.}

\description{
  This function performs a Wald test of sampling independence on item severity parameters.
}

\usage{
EWaldtest(b1, b2, se1, se2)
}
\arguments{
  \item{b1}{Item severity vector estimated on the first sample subset.}
  \item{b2}{Item severity vector estimated on the second sample subset (complementary to the one used to estimate \code{b1}).}
    \item{se1}{Item standard errors (first sub-sample).}
  \item{se2}{Item standard errors (second sub-sample).}
}
\details{
   This function computes the Wald test to check if the sampling independence assumption holds. The test statistics is 
   \deqn{
   z = \frac{\hat{b}_1 - \hat{b}_2}{\sqrt{Var(\hat{b}_1) + Var(\hat{b}_2)}},
   }
   where \eqn{\hat{b}_1} and \eqn{\hat{b}_2} are the ML item parameter estimates computed on two, randomly extracted and complementary, sub-samples of the overall sample, and \eqn{Var(\hat{b}_1)} and \eqn{Var(\hat{b}_2)} the corresponding variances.
   
   The null hypothesis is of parameter estimates equality while the alternative hypothesis is bilateral (\eqn{b_1 \neq b_2}).
   
   The sampling independence assumption implies that item severity parameter estimates do not depend on the analysed samples, but only on the severity of the items themselves. One way to check for this assumption is to extract two (complementary) random samples from the original sample and run the Wald test, or, more accurately, to extract \code{B} random sub-samples and observe the p-values distribution for each item. 
}
\value{
  A list with the following elements: 
  \tabular{ll}{
	\code{z} \tab The z-statistics of the Wald test. \cr
	\code{p} \tab Computed p-value of the Wald test. \cr
	\code{tab} \tab Descriptive table with test results. \cr 
	}
}
\author{Sara Viviani \email{sara.viviani@fao.org}}

\seealso{
  \code{\link{RM.w}}
}
\examples{
\dontrun{
data(data.FAO_country3)
# Questionnaire data and weights
XX.country3 = data.FAO_country3[,1:8]
wt.country3 = data.FAO_country3$wt

# Split the sample in two random sub-samples
n = nrow(XX.country3)
samp1 = sample(1:n, n/2)
samp2 = setdiff(1:n, samp1)
# Fit two Rasch models on the two sub-samples
rr1 = RM.w(XX.country3[samp1, ], wt.country3[samp1])
rr2 = RM.w(XX.country3[samp2, ], wt.country3[samp2])
# Test sampling independence
EWaldtest(rr1$b, rr2$b, rr1$se.b, rr2$se.b)$tab
}
}

