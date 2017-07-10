\name{MitISEM}
\alias{MitISEM}
\title{Mixture of Student-\eqn{t}{t} distributions using Importance Sampling 
weighted Expectation Maximization steps
}
\description{
Approximates a \eqn{k}{k} dimensional function/kernel by a mixture of 
student-\eqn{t}{t} distributions using Importance Sampling weighted Expectation 
Maximization steps.
}
\usage{
MitISEM(KERNEL,mu0,Sigma0=NULL,df0=1,mit0=NULL,control=list(),...)
}
\arguments{
  \item{KERNEL}{Function to be approximated. First argument should be the parameter matrix. 
  A \code{log} argument should exist such that the function returns log-density if \code{log=TRUE}
  }
  \item{mu0}{vector of length \eqn{k}, starting points for approximation
}
  \item{Sigma0}{(optional) \eqn{k \times k}{kxk} dimensional positive definite symmetric initial scale matrix. Default: matrix is obtained by the \code{BFGS} algorithm.
}
  \item{df0}{(optional) double \eqn{>0} initial degree of freedom of the Student-t density. Default: \eqn{df0=1}.
}
  \item{mit0}{(optional) initial mixture density defined. See *Details*.
}
  \item{control}{(optional)
   list of control parameters for IS and EM optimization and stopping rule of the \code{H} component mixture of t densities. See *Details*
}
   \item{\dots}{other arguments to be passed to \code{KERNEL}
}
}
\value{
list containing:
    \item{mit}{(list) optimal mixture density  with \eqn{H}{H} mixture components. See *Details*.}
    \item{CV}{vector of length \eqn{H}{H} with coefficient of variation obtained from each addition of mixture components.}
     \item{time}{(double) processed time.}
     \item{summary}{(data.frame) summary information on construction of components, processed time and CV.}
}
\references{%
Basturk, N., Grassi, S., Hoogerheide, L., Opschoor, A. and Van Dijk, H. K. (2017) 
The R Package MitISEM: Efficient and Robust Simulation Procedures for Bayesian Inference. 
\emph{Journal of Statistical Software}, 79(1): 1-39. \doi{10.18637/jss.v079.i01}.

Hoogerheide L., Opschoor, A. and Van Dijk, H. K. (2012) 
A Class of Adaptive Importance Sampling Weighted {EM} Algorithms for Efficient and Robust Posterior and Predictive Simulation. 
\emph{Journal of Econometrics}, 171(2): 101-120.
\url{http://www.sciencedirect.com/science/article/pii/S0304407612001583}.
}

\seealso{\code{\link{isMit}}}
\details{
Providing \code{mit0} argument makes arguments \code{mu0}, \code{Sigma0} and \code{df0} obsolete. Argument \code{mit0} (if provided) should include the following components (see \code{\link{isMit}}):
 \describe{%
    \item{\code{p}}{vector (of length \eqn{H}) of mixture probabilities.}
    \item{\code{mu}}{matrix (of size \eqn{H \times k}{Hxk}) containing
      the vectors of modes (in row) of the mixture components.}
    \item{\code{Sigma}}{matrix (of size \eqn{H \times k^2}{Hx(kxk)})
      containing the scale matrices (in row) of the mixture components.}
    \item{\code{df}}{vector (of length \eqn{H}) degree of freedom parameters of the Student-t components. Each element should be above \eqn{0}}
  }

Value \code{mit} has the same structure as \code{mit0}, where \eqn{H}{H} and parameters of the mixture density are optimized. 

  Optional argument \code{control} can provide several optimization parameters:
\describe{%
    \item{\code{N}}{integer (>100) number of draws used in the simulation. Default: \code{N=1e5}.}
    \item{\code{robust.N}}{logical indicating if robust draws are used if \code{robust.N=TRUE} (default), simulations are repeated to get \code{N} draws with finite \code{KERNEL} values.}
    \item{\code{Hmax}}{ integer\eqn{>0} maximum number of components. Default: \code{H=10}.}
    \item{\code{StopMethod}}{string, \code{CV} (default) or \code{AR} defining type of stopping criteria for MitISEM approximation. \code{CV} method stops the algorithm if the coefficient of variation in IS weights converges. \code{AR} method stops the algorithm if the (expected) acceptance rate given the current MitISEM approximation and function \code{KERNEL} converges.}
      \item{\code{CVtol}}{double in \eqn{(0,1)} convergence criteria for \code{CV} method. Higher values lead to faster convergence but worse approximation. Default: \eqn{CVtol=0.1}, algorithm stops if \eqn{StopMethod=CV} and the change in coefficient of variation is below 10\eqn{\%}.}
      \item{\code{ARtol}}{double in \eqn{(0,1)} convergence criteria similar to \code{CVtol}, used if \eqn{StopMethod=AR}. Default: \eqn{ARtol=0.1}.}
      \item{\code{trace}}{logical to print partial output. Default: \eqn{trace=FALSE}, no tracing information.}
      \item{\code{trace.init}}{logical to print output of the first student-t optimization. Default: \eqn{trace=FALSE}, no tracing information.}
      \item{\code{maxit.init}}{double, maximum number of iterations in the first student-t optimization. Default: \eqn{maxit.init=1e4}.}
	\item{\code{reltol.init}}{double, relative tolerance in the first student-t optimization. Default: \eqn{reltol.init=1e-8}.}
	\item{\code{maxit.EM}}{integer\eqn{>0}, maximum number of iterations for the EM algorithm. Default: \eqn{maxit.EM=1000}.}
	\item{\code{tol.EM}}{double\eqn{>0}, tolerance for EM steps' convergence, Default: \eqn{tol.EM=0.001}.}
	\item{\code{trace.EM}}{logical to print partial output during the IS-EM algorithm. Default: \eqn{trace.EM=FALSE}, no tracing information.}
	\item{\code{optim.df}}{logical to optimize degrees of freedom of the Student-t components. Default: \eqn{optim.df=TRUE} df are optimized. Note: Keeping degrees of freedom in low values may be desired if the approximation is used in a rejection sampling. If \eqn{optim.df=FALSE}, degree of freedom of all student t components are fixed at \code{df0}.}
       \item{\code{inter.df}}{increasing vector of length \eqn{2} range of search values for df optimization, active if \eqn{optim.df=TRUE}. Default: \eqn{inter.df= (0.01,30)}.}
       \item{\code{tol.df}}{double \eqn{>0}, tolerance for degree of freedom optimization, active if \eqn{optim.df=TRUE}. Default: \eqn{tol.df=0.0001}}
	\item{\code{maxit.df}}{integer \eqn{>0} maximum number of iterations for degree of freedom optimization, active if \eqn{optim.df=TRUE}. Default: \eqn{maxit.df=1e3}.}
	\item{\code{trace.df}}{logical to print partial output during degree of freedom optimization, active if \eqn{optim.df=TRUE}. Default: \eqn{trace.df=FALSE}.} 
	\item{\code{tol.pr}}{double in [0,1), minimum probability required to keep mixture components. Default: \eqn{tol.pr=0} all mixture components are kept.}
	\item{\code{ISpc}}{double in (0,1), fraction of draws to construct new component. Default: \eqn{ISpc=0.1}.}
	\item{\code{Pnc}}{double in (0,1), initial probability of the new component. Default: \eqn{Pnc=0.1}.}
  }
}
\examples{
require(graphics)
set.seed(1234);
# define Gelman Meng Kernel
GelmanMeng <- function(x, A = 1, B = 0, C1 = 3, C2 = 3, log = TRUE){
	    if (is.vector(x))
	    x <- matrix(x, nrow = 1)
	    r <- -.5 * (A * x[,1]^2 * x[,2]^2 + x[,1]^2 + x[,2]^2
                 - 2 * B * x[,1] * x[,2] - 2 * C1 * x[,1] - 2 * C2 * x[,2])
	    if (!log)
	    r <- exp(r)
	    as.vector(r)
}
# get MitISEM approximation
mu0 <-c(3,4)
app.MitISEM <- MitISEM(KERNEL=GelmanMeng,mu0=mu0,control=list(N=2000,trace=TRUE))
mit=app.MitISEM$mit

# plot approximation (components and full approximation)
MitISEM.plot.comps <- function(mit,x1,x2,log=FALSE){
  Mitcontour <- function(x1,x2,mit,log=FALSE){
    dmvgt(cbind(x1,x2),mit=mit,log=log)
  } 
  H <- length(mit$p)
  K <- ncol(mit$mu)
  cols <- 1:H
  for (h in 1:H){
    mit.h  <-list(p=1,mu=matrix(mit$mu[h,],1,K),
                Sigma=matrix(mit$Sigma[h,],1,(K^2)),df=mit$df[h])
    z      <- outer(x1,x2,FUN=Mitcontour,mit=mit.h)
    contour(x1,x2,z,col=h,lty=h,labels="",add=(h!=1),
            xlab="x1",ylab="x2",main='MitISEM approximation components')
  }	
  legend("topright",paste("component ",1:H),lty=cols,col=cols,bty='n')
  z <- outer(x1,x2, FUN=Mitcontour,mit=mit)
  image(x1,x2,z,las=1,col=gray((20:0)/20),main='MitISEM approximation')
}
x1 <- seq(-2,6,0.05)
x2 <- seq(-2,7,0.05)
MitISEM.plot.comps(mit,x1,x2)

\dontrun{
  # Bayesian inference of the GARCH model using MitISEM and Importance Sampling
  library(AdMit) # required for Importance Sampling
  library(tseries) # required for loading the data
  # load data : downloaded on 2013/01/18
  prices <- as.vector(get.hist.quote("^GSPC",quote="AdjClose",start="1998-01-02",end="2002-12-26"))
  data  <- 100 * (prices[-1] - prices[-length(prices)]) / (prices[-length(prices)]) 
  prior.GARCH<-function(omega,beta,alpha,
                        mu,log=TRUE){
    c1 <- (omega>0 & omega <1 & beta>=0 & alpha>=0)
    c2 <- (beta + alpha< 1)
    c3 <- (mu>-1 & mu<1)
    r1 <- c1 & c2 & c3
    r2 <- rep.int(-Inf,length(omega))
    r2[r1==TRUE] <- 0
    if (!log)
      r2 <- exp(r2)
    cbind(r1,r2)
  }
  post.GARCH <- function(theta,data,h1,log=TRUE){
    if (is.vector(theta))
      theta <- matrix(theta, nrow = 1)
    omega <- theta[,1]
    beta <- theta[,2]
    alpha <- theta[,3]
    mu <- theta[,4]
    N <- nrow(theta)
    pos <- 2:length(data)
    prior <- prior.GARCH(omega=omega,beta=beta,alpha=alpha,mu=mu)
    d <- rep.int(-Inf,N)
    for (i in 1:N){
      if (prior[i,1] == TRUE){
        h <- c(h1, omega[i] + alpha[i] * (data[pos-1]-mu[i])^2)
        for (j in pos){
          h[j] <- h[j] + beta[i] * h[j-1]
        }
        tmp <- dnorm(data[pos],mu[i],sqrt(h[pos]),log=TRUE)
        d[i] <- sum(tmp) + prior[i,2]
      }
    }
    if (!log) d <- exp(d)
    as.numeric(d)
  }
  theta <- c(.08, .86, .02, .03) # initial parameters for MitISEM
  names(theta)<-c("omega","beta","alpha","mu")
  h1 <- var(data) # initial data variance
  # MitISEM GARCH approximation
  cat("MitISEM GARCH results",fill=TRUE)
  cat('--------------------------',fill=TRUE)
  set.seed(1111)
  app.GARCH <- MitISEM(KERNEL=post.GARCH,
                       mu0=theta, control=list(trace=TRUE),h1=h1,
                       data=data)
  print(app.GARCH$summary)
  # Importance Sampling using MitISEM candidate
  cat('Importance Sampling result from MitISEM candidate',fill=TRUE)
  cat('---------------------------------------------------',fill=TRUE)
  set.seed(1111)
  IS.MitISEM.GARCH <- AdMitIS(N = 10e4,data=data,h1=h1,
                              KERNEL=post.GARCH,mit=app.GARCH$mit)
  print(IS.MitISEM.GARCH)
}

}
