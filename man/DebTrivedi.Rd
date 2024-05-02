\name{DebTrivedi}
\alias{DebTrivedi}
\docType{data}
\title{
DebTrivedi Dataset
}
\description{
The DebTrivedi dataset contains data on health care utilization and insurance coverage.
}
\usage{data("DebTrivedi")}
\format{
  A data frame with 4406 observations on the following 19 variables.
  \describe{
    \item{\code{ofp}}{a numeric vector indicating physicians office visits}
    \item{\code{ofnp}}{a numeric vector}
    \item{\code{opp}}{a numeric vector}
    \item{\code{opnp}}{a numeric vector}
    \item{\code{emer}}{a numeric vector}
    \item{\code{hosp}}{a numeric vector indicating number of hospital stays}
    \item{\code{health}}{a factor with levels \code{poor} \code{average} \code{excellent} indicating self-perceived health status}
    \item{\code{numchron}}{a numeric vector indicating number of chronic conditions}
    \item{\code{adldiff}}{a factor with levels \code{no} \code{yes}}
    \item{\code{region}}{a factor with levels \code{midwest} \code{noreast} \code{other} \code{west}}
    \item{\code{age}}{a numeric vector}
    \item{\code{black}}{a factor with levels \code{no} \code{yes}}
    \item{\code{gender}}{a factor with levels \code{female} \code{male}}
    \item{\code{married}}{a factor with levels \code{no} \code{yes}}
    \item{\code{school}}{a numeric vector indicating number of years of education}
    \item{\code{faminc}}{a numeric vector}
    \item{\code{employed}}{a factor with levels \code{no} \code{yes}}
    \item{\code{privins}}{a factor with levels \code{no} \code{yes} indicating private insurance indicator}
    \item{\code{medicaid}}{a factor with levels \code{no} \code{yes}}
  }
}
\details{
Deb and Trivedi (1997) conducted an analysis on data comprising 4406 individuals aged 66 and above, enrolled in Medicare, a public insurance program. These data were initially sourced from the US National Medical Expenditure Survey (NMES) for the years 1987/88. They are accessible through the data archive of the Journal of Applied Econometrics at http://qed.econ.queensu.ca/jae/1997-v12.3/deb-trivedi/. Additionally, the dataset was compiled for an R package associated with Kleiber and Zeileis (2008) and is also retrievable as DebTrivedi.rda from the Journal of Statistical Software, alongside Zeileis (2006). The primary aim is to develop models for medical care demand, encompassing physician/non-physician office and hospital outpatient visits, based on available patient covariates.
}
\source{
http://www.jstatsoft.org/v27/i08/paper
}
\references{
Zeileis, A. and Kleiber, C. and Jackma, S. (2008). "Regression Models for Count Data in R". JSS 27, 8, 1–25.
}
\examples{
data(DebTrivedi)
plot(table(DebTrivedi$ofp), main = "physicians office visit (N=4406)", ylab="frequency")
#the response variable “physicians office visits” is highly right skewed
}
\keyword{datasets}
