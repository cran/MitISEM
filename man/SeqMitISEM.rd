\name{SeqMitISEM}
\alias{SeqMitISEM}
\title{Sequential approximation using Mixture of 
Student-\eqn{t}{t} distributions using Importance Sampling 
weighted Expectation Maximization steps
}
\description{
Approximates a \eqn{k}{k} dimensional function/kernel 
using mixture of student-\eqn{t}{t} distributions for the 
initial data sample and updated data samples sequentially
}
\usage{
SeqMitISEM(data,KERNEL,mu0,Sigma0=NULL,df0=1,control.MitISEM=list(),control.seq=list(),...)
}
\arguments{
  \item{data}{matrix (size \eqn{T\times m}{Txm}) with data values, \eqn{T} observations and \eqn{m} data series
  }
  \item{KERNEL}{Function/posterior to be approximated. \code{data} and \code{log} arguments should exis. The function must return log-density if \code{log=TRUE}. All data should be loaded in argument \code{data}
  }
  \item{mu0}{vector of length \eqn{k} starting points. They should be defined as in \code{\link{MitISEM}}
}
  \item{Sigma0,df0}{(optional) initial scale and degrees of freedom for the student t density. They should be defined as in \code{\link{MitISEM}} 
}
  \item{control.MitISEM}{(optional)
   control parameters passed to \code{MitISEM}. See \code{\link{MitISEM}}.
}
   \item{control.seq}{
   control parameters for sequential updating of the \code{MitISEM} approximation. See *Details*.
}
   \item{\dots}{other arguments to be passed to \code{KERNEL}.
}
}
\seealso{\code{\link{MitISEM}}}
\details{
  The optional argument \code{control.seq} can provide several optimization parameters:
\describe{%
    \item{\code{T0}}{integer (\eqn{<T}) number of observations. Default: \code{round(T/2)}.}
    \item{\code{tau}}{vector of length \eqn{t} with iterative number of observations to add to the sample. Its elements should be positive integers, and \eqn{T0+max(tau)<=T} should hold. Default: \code{tau=1}, one single observation is added to the sample for Sequential MitISEM.
    }
    \item{\code{tol.seq}}{ double in \eqn{(0,1)} convergence criteria for sequential Coefficient of Variation convergence
	Default: \code{tol.seq=0.2}.
    }
    \item{\code{method}}{\eqn{{0,1,2}} method to select initial data sample. 
       if \code{method=0} initial sample is randomly selected.
       if \code{method=1} first \code{T0} observations are taken as initial sample (default).
       if \code{method=2} last \code{T0} observations are taken as initial sample.
   }
    \item{\code{trace}}{logical to print partial output. Default: \eqn{trace=FALSE}, no tracing information.
   }
   }
}
\references{
Hoogerheide, L. and Opschoor, A. and Van Dijk, H. K. (2012).
A Class of Adaptive Importance Sampling Weighted {EM} Algorithms for Efficient and Robust Posterior and Predictive Simulation. \emph{Journal of Econometrics}, in press.
\url{http://www.sciencedirect.com/science/article/pii/S0304407612001583}
}
