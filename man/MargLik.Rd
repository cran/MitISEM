\name{MargLik}
\alias{MargLik}
\title{Marginal Likelihood calculation using Importance Sampling
and mixture of Student-\eqn{t}{t} densities as candidate
}
\description{
Calculation of marginal likelihoods using Importance Sampling, with a 
Mixture of Student-\eqn{t}{t} candidate density.
Calculated marginal likelihoods from two data samples can be used to get predictive likelihoods using 
\code{\link{PredLik}}.
}
\usage{
MargLik(N=1e4,mit,KERNEL,...)
}
\arguments{
  \item{N}{integer \eqn{> 100} number of draws for Importance Sampling
}
  \item{mit}{Mixture of Student-\eqn{t}{t} density for the full sample, list describing the mixture of Student-t. See \code{\link{isMit}}. The mixture density can be obtained from \code{\link{MitISEM}} or \code{\link{SeqMitISEM}}
}
\item{KERNEL}{Posterior kernel to be approximated. See *Details*. 
A \code{log} argument must exist. 
The function must return log-density if \code{log=TRUE}.
}
   \item{\dots}{other arguments to be passed to \code{KERNEL}
}
}
\value{%
  list containing:
  \item{ML.mean}{Marginal likelihood (posterior mean) \eqn{\times 10^{scale}}{x 10^{scale}}
  }
  \item{ML.NSE}{Numerical Standard Error for mean Marginal likelihood \eqn{\times 10^{scale}}{x 10^{scale}}
  }
  \item{scale}{integer \eqn{> 0}{> 0} providing the scaling for predictive likelihood.   
  (scaling may be necessary for numerical accuracy)
  }
}
\details{If \code{MargLik} is used to get the Marginal Likelihood of a single model, \code{KERNEL} must
be the exact posterior density (including the scaling constant) of parameters. 

If \code{MargLik} is used as an intermediate step, for instance for calculating predictive likelihoods,
\code{KERNEL} can be a posterior kernel or the exact posterior density of parameters.  See \code{\link{PredLik}}.
}
\examples{
mit  <- list(p=1,mu=matrix(1),Sigma=matrix(0.1),df=5)
data <- rnorm(100,1)
KERNEL <- function(theta,data,log=TRUE){
  if(is.vector(theta))
    theta = matrix(theta,nrow=1)
  r <- apply(theta,1,function(x,data)(sum(dnorm(data,x,log=log))),data=data)
}
MargLik(N=1000,mit=mit,KERNEL=KERNEL,data=data)
}
\seealso{\code{\link{isMit},\link{PredLik},\link{MitISEM},\link{SeqMitISEM}}}
