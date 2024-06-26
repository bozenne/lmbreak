% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.R
\name{anova.lmbreak}
\alias{anova.lmbreak}
\title{Residual Sum of Squares}
\usage{
\method{anova}{lmbreak}(object, indiv = FALSE, psi = NULL, ...)
}
\arguments{
\item{object}{[lmbreak] breakpoint model}

\item{indiv}{[logical] should the contribution of each observation to the RSS be output.
Otherwise outputs the sum of the contributions.}

\item{psi}{[numeric vector] breakpoint values.}

\item{...}{Not used. For compatibility with the generic method.}
}
\value{
A numeric value (\code{indiv=FALSE}) or a numeric vector (\code{indiv=TRUE}) where each element corresponds to an individual.
}
\description{
Evaluate the residual sum of squares (RSS) at given breakpoints.
}
\examples{
set.seed(10)
df <- simBreak(c(1,100), breakpoint = c(0,50,100), slope = c(1,-1))
anova(lm(Y ~ X, data = df))

anova(lmbreak(Y ~ bp(X, "11"), data = df))
}
\keyword{methods}
