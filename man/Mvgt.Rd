\name{Mvgt}
\alias{dmvgt}
\alias{rmvgt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{General student t distribution
%%  ~~function to do ... ~~
}
\description{
Density and random generation for the general studen t distribution
}
\usage{
dmvgt(theta, mit = list(), log=TRUE)
rmvgt(N,mit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{theta}{
  Vector of lenght \eqn{N} or \eqn{N\times k}{Nxk} matrix of quantiles. If \code{theta} is a matrix, each row is taken to be a quantile.
}
  \item{N}{
   number of observations
}
  \item{mit}{
   list defining the mixture components. See \code{\link{isMit}} for how it should be defined.
}
  \item{log}{
   logical; if TRUE (default), probabilities p are given as \eqn{\ln(p)}{ln(p)}.
}
}
\value{
   \code{dmvgt} returns vector of size \code{N} with density values for each row of \code{theta}
   
   \code{rmvgt} returns an \eqn{N\times k}{Nxk} matrix of draws from the \eqn{k}-variate mixture of student t densities
}
\seealso{\code{\link{isMit}}}
\examples{
  H      <- 2
  p      <- runif(H)
  p      =  p / sum(p) 
  mu     <- matrix(seq(1:H),H,1)
  Sigma  <- matrix(runif(H^2),H,H)
  df     <- seq(1:H)  
  Ndraws <- rmvgt(N=10,mit=list(p=p,mu=mu,Sigma=Sigma,df=df))
  pdraws <- dmvgt(theta=Ndraws,mit=list(p=p,mu=mu,Sigma=Sigma,df=df))
}
