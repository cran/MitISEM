\name{Mit}
\alias{isMit}
\title{The 'mit' object}
\description{Function to check if \code{mit} generalized mixture of t densities is well-defined. The \code{mit} object is designed to be used in the rest of the \code{MitISEM} package functions
}
\usage{
isMit(mit)
}
\arguments{
  \item{mit}{
  an object to be tested
}
}
\value{%
  logical, \code{TRUE} if \code{mit} definition is correct, \code{FALSE} otherwise
}

\details{
  Argument \code{mit} is a list describing the mixture of Student-t distributions with the following components:
  \describe{%
    \item{\code{p}}{vector (of length \eqn{H}) of mixture probabilities.}
    \item{\code{mu}}{matrix (of size \eqn{H \times d}{Hxd}) containing
      the vectors of modes (in row) of the mixture components.}
    \item{\code{Sigma}}{matrix (of size \eqn{H \times d^2}{Hx(dxd)})
      containing the scale matrices (in row) of the mixture components.}
    \item{\code{df}}{vector (of length \eqn{H}) degree of freedom parameters for each Student-t component (double \eqn{>0}).}
  }
}
\examples{
  # a correct Mit definition returns 'TRUE'
  H      <- 2
  p      <- runif(H)
  p      =  p / sum(p) 
  mu     <- matrix(seq(1:H),H,1)
  Sigma  <- matrix(runif(H^2),H,H)
  df     <- seq(1:H)
  isMit(mit=list(p=p,mu=mu,Sigma=Sigma,df=df))

  # an incorrect Mit definition returns 'FALSE'
  mu   = t(mu)
  isMit(mit=list(p=p,mu=mu,Sigma=Sigma,df=df))
}
