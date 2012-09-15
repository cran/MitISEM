\name{PredLik}
\alias{PredLik}
\title{Predictive Likelihood calculation using Importance Sampling
and mixture of Student-\eqn{t}{t} densities as candidate
}
\description{
Calculation of predictive likelihoods using Importance Sampling, 
given subsample and  full data sample and Mixture of Student-\eqn{t}{t} 
candidate density. Predictive likelihood is calculated using the 
marginal likelihood from full sample and subsample. 
See \code{\link{MargLik}}.
}
\usage{
PredLik(N=1e4,mit.fs,mit.ss,KERNEL,data.fs,data.ss,...)
}
\arguments{
  \item{N}{integer \eqn{> 100} number of draws for Importance Sampling
}
  \item{mit.fs}{Mixture of Student-\eqn{t}{t} density for the full sample, list describing the mixture of Student-t. See \code{\link{isMit}}. The mixture density can be obtained from \code{\link{MitISEM}} or \code{\link{SeqMitISEM}}
}
  \item{mit.ss}{Mixture of Student-\eqn{t}{t} density for subsample. Must be defined as \code{mit.fs}.
}
\item{KERNEL}{Function/posterior to be approximated. 
\code{data} and \code{log} arguments must exist. The function must return log-density if \code{log=TRUE}. All data should be loaded in argument \code{data}
}
\item{data.fs}{Full data, vector (length \eqn{T1}{T1}) or matrix (size \eqn{T1\times m}{T1xm}) with data values, \eqn{T1} observations and \eqn{m} data series.
}
\item{data.ss}{Sample of data, vector (length \eqn{T2}{T2}) or matrix (size \eqn{T2\times m}{T2xm}) with data values, \eqn{T2} observations and \eqn{m} data series. \eqn{T2 <T1}{T2 <T1}.
}
   \item{\dots}{other arguments to be passed to \code{KERNEL}
}
}
\value{%
  list containing:
  \item{PL}{Predictive likelihood \eqn{\times 10^{scale}}{x 10^{scale}}
  }
  \item{scale}{integer \eqn{> 0}{> 0} providing the scaling for predictive likelihood.   
  (scaling may be necessary for numerical accuracy)
  }
}
\details{
  Argument \code{KERNEL}
}
\references{
Eklund, J. and Karlsson, S. (2007).
Forecast combination and model averaging using predictive measures. \emph{Econometric Reviews}, 26, 329-363.

Min, C. and Zellner, A. (1993).
Bayesian and non-Bayesian methods for combining models and forecasts with applications to forecasting international growth rates.
\emph{Journal of Econometrics}, 56, 89-118.
}
\seealso{\code{\link{isMit},\link{MargLik},\link{MitISEM},\link{SeqMitISEM}}}
