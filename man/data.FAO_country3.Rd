\name{data.FAO_country3}
\alias{data.FAO_country3}
\title{Food insecurity data for a GWP pilot country (Country3).}
\description{
The dataset includes the FIES data (Food Insecurity Experience based Scale), sampling weights, and some demographic variables for a sample country (Country1). Data have been collected by the Gallup World Poll. 
}

\arguments{
  \item{WORRIED}{FIES question n. 1. During the past 12 months, was there a time when you were worried you would run out of food because of 
  	lack of money or other resources?}
  \item{HEALTHY}{FIES question n. 2. During the past 12 months, was there a time when you were unable to eat healty and nutritious food because of 
  	lack of money or other resources?} 
  \item{FEWFOOD}{FIES question n. 3. During the past 12 months, was there a time when you ate only few kinds of food because of 
  	lack of money or other resources?}
  \item{SKIPPED}{FIES question n. 4. During the past 12 months, was there a time when you had to skip a meal because of 
  	lack of money or other resources?}
  \item{ATELESS}{FIES question n. 5. During the past 12 months, was there a time when you ate less than you thought you should because of 
  	lack of money or other resources?}
  \item{RUNOUT}{FIES question n. 6. During the past 12 months, was there a time when your household ran out of food because of 
  	lack of money or other resources?}
  \item{HUNGRY}{FIES question n. 7. During the past 12 months, was there a time when you were hungry but did not eat because of 
  	lack of money or other resources?}
  \item{WHLDAY}{FIES question n. 8. During the past 12 months, was there a time when you went without eating for a whole day because of 
  	lack of money or other resources?}
  \item{gender}{Respondents' gender.}
  \item{age}{Respondents' age.}
  \item{education}{Respondents' education.}
  \item{urbanrural}{Respondents' area (urban or rural).}
  \item{wt}{Sampling weights.}
}
\author{Sara Viviani \email{sara.viviani@fao.org}}
\examples{
\dontrun{
data(data.FAO_country3)
# Matrix of questionnaire data and weights (for the Rasch analysis)
XX.country3 = data.FAO_country3[,1:8]
wt.country3 = data.FAO_country3$wt
}
}