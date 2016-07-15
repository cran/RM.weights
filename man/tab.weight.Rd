\name{tab.weight}
\alias{tab.weight}
\title{Rasch descriptives and weighted tables.}
\description{
This function computes the main descriptive tables (number of affirmative responses to each item, weighted and unweighted distribution across raw scores, etc.), weighted and unweighted, for Rasch scale items and respondents. It also computes (univariate or bivariate) tables for user-specified variables.
}
\usage{
tab.weight(variab, wt, XX = NULL)
}
\arguments{
  \item{variab}{It could be a single variable, or a list of two variables, 
  of \code{factor} type.}
  \item{wt}{Vector of sampling weights. \code{length(variab) = length(wt)}.}
  \item{XX}{0/1 data matrix or data frame of responses to the questionnaire.
  Affirmative responses must be coded as 1s. If left unspecified, basic descriptive 
  statistics for Rasch analysis (number of Yes per item, weighted and unweighted 
  distribution across raw scores, etc.) are not computed.}
}
\value{
  \item{RS.abs.w}{Absolute weighted distribution across raw scores. Only computed when \code{XX} is not null.}
  \item{RS.rel.w}{Relative weighted distribution across raw scores. Only computed when \code{XX} is not null.}
  \item{RS.abs}{Absolute unweighted distribution across raw scores. Only computed when \code{XX} is not null.}
  \item{RS.rel}{Relative unweighted distribution across raw scores. Only computed when \code{XX} is not null.}
    \item{Perc.Yes}{Proportion of Yes responses to each item. Only computed when \code{XX} is not null.}
      \item{Perc.Yes.w}{Weighted proportion of Yes responses to each item. Only computed when \code{XX} is not null. Only computed when \code{XX} is not null.}
    \item{rv}{Vector of raw scores.}
     \item{tab.ext.w}{Weighted table for the variable \code{variab}. If \code{variab} is a list, a bivariate frequency table is computed.}
}


\author{Sara Viviani \email{sara.viviani@fao.org}}
\examples{
\dontrun{
data(data.FAO_country1)
XX.country1 = data.FAO_country1[,1:8]
wt.country1 = data.FAO_country1$wt
gender = data.FAO_country1$gender
urbanrural = data.FAO_country1$urbanrural
# Univariate weighted table by gender
tab.weight(gender, wt.country1)$tab.ext.w

# Bivariate weighted table by gender and urban/rural
tab.weight(list(gender,urbanrural), wt.country1)$tab.ext.w

# Fit Rasch descriptives
fit.descr = tab.weight(wt = wt.country1, XX = XX.country1)
# Weighted distribution across raw-scores (absolute and relative)
cbind("Abs.RS distrib." = fit.descr$RS.abs.w, 
"Rel.RS distrib." = fit.descr$RS.rel.w)
# Weighted and unweighted percentage of Yes per item
cbind("Weighted perc. of Yes" = fit.descr$Perc.Yes.w, 
      "Unweighted perc. of Yes" = fit.descr$Perc.Yes)
}
}


