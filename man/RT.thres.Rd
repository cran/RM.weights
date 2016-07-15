\name{RT.thres}
\alias{RT.thres}
\title{Rasch-Thurstone thresholds for trichotomous scales under CML.}
\description{
  Computes the Rasch-Thurstone thresholds for the CML estimation, starting from Rasch thresholds of a trichotomous scale. 
}
\usage{
RT.thres(R.thres)
}
\arguments{
  \item{R.thres}{Vector of rasch thresholds of a trichotomous scale, i.e. points on the latent trait corresponding to which the caracteristic curves of two consecutive categories cross.}
}
\details{
The Rasch-Thurstone (R-T) threshold for a given response category of an item is defined as the point on the latent trait at which the combined probability of response in that category and all more severe categories is .5. Conceptually, an R-T threshold for a given response category for a trichotomous item is comparable with the item severity parameter estimated for a dichotomous item coded 1 for that and all more severe categories and 0 for all less severe categories.
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{b1.rt} \tab R-T threshold for the first category. \cr
	\code{b2.rt} \tab R-T threshold for the second category. \cr 
	\code{b12.rt} \tab Vector of R-T thresholds for the first and second category. \cr 
	}
}

\author{Sara Viviani \email{sara.viviani@fao.org}, Mark Nord \email{mark.nord@fao.org}}



