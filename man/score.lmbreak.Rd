% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.R
\name{score.lmbreak}
\alias{score.lmbreak}
\title{Breakpoint Score Equation}
\usage{
\method{score}{lmbreak}(x, indiv = FALSE, psi = NULL, ...)
}
\arguments{
\item{x}{[lmbreak] breakpoint model}

\item{indiv}{[logical] should the contribution of each observation to the score be output.
Otherwise outputs the sum of the contributions.}

\item{psi}{[numeric vector] breakpoint values.}

\item{...}{Not used. For compatibility with the generic method.}
}
\value{
A numeric vector (\code{indiv=FALSE}) where each element corresponds to a breakpoint.
or a matrix (\code{indiv=TRUE}) where each row corresponds to an individual and each column to a breakpoint.
}
\description{
First derivative of the residual sum of squares with respect to the breakpoints.
}
\examples{
if(require(lava)){
set.seed(10)
df <- simBreak(c(1,100), breakpoint = c(0,50,100), slope = c(1,-1))
score(lmbreak(Y ~ bp(X, "11"), data = df))
}
}
\keyword{methods}
