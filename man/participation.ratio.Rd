% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{participation.ratio}
\alias{participation.ratio}
\title{Calculate the participation ratio given the calculate Schmidt modes.}
\usage{
participation.ratio(modes)
}
\arguments{
\item{modes}{Schmidt modes}
}
\value{
participation ratio
}
\description{
Calculate the participation ratio given the calculate Schmidt modes.
}
\examples{
singlet<- matrix(c(0, sqrt(0.7), sqrt(0.3), 0), byrow = TRUE, nrow = 2)
modes<- schmidt.decompose(singlet)
participation.ratio(modes)

}
