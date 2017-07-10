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
SeqMitISEM(data,KERNEL,mu0,Sigma0=NULL,df0=1,control.MitISEM=list(),control.seq=list(),
...)
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
Basturk, N., Grassi, S., Hoogerheide, L., Opschoor, A. and Van Dijk, H. K. (2017) 
The R Package MitISEM: Efficient and Robust Simulation Procedures for Bayesian Inference. 
\emph{Journal of Statistical Software}, 79(1): 1-39. \doi{10.18637/jss.v079.i01}.

Hoogerheide L., Opschoor, A. and Van Dijk, H. K. (2012) 
A Class of Adaptive Importance Sampling Weighted {EM} Algorithms for Efficient and Robust Posterior and Predictive Simulation. 
\emph{Journal of Econometrics}, 171(2): 101-120.
\url{http://www.sciencedirect.com/science/article/pii/S0304407612001583}.
}
\examples{
\dontrun{
  # Sequential MitISEM application for SP500 data
  # Calculates 50 predictive likelihoods for the mGARCH(1,1) model, SP500 data
  # For details see: 'The R package MitISEM: Efficient and Robust Simulation Procedures 
  # for Bayesian Inference' by N. Basturk, S. Grassi, L. Hoogerheide, A. Opschoor, 
  # H.K. van Dijk.
  library(tseries)
  source("PostmGARCH.R") # posterior of the model under flat priors
  
  # load data
  prices <- as.vector(get.hist.quote("^GSPC",quote="AdjClose",start="1998-01-02",
    end="2002-12-26"))
  y <- 100 * (prices[-1] - prices[-length(prices)]) /  (prices[-length(prices)])
  # Prior and posterior densities for the mixture of GARCH(1,1) model with 
  # 2 mixture components
  prior.mGARCH<-function(omega, lambda, beta, alpha, p, mu, log=TRUE){
    c1 <- (omega>0 & omega<1 & beta>=0 & alpha>=0)
    c2 <- (beta + alpha< 1)
    c3 <- (lambda>=0 & lambda<=1)
    c4 <- (p>0.5 & p<1)
    c5 <- (mu>-1 & mu<1)
    r1 <- c1 & c2 & c3 & c4 & c5  
    r2 <- rep.int(-Inf,length(omega))
    tmp <- log(2) # ln(1 / ( p(beta,alpha) * p(p) * p(mu))
    r2[r1==TRUE] <- tmp
    if (!log)
      r2 <- exp(r2)
    cbind(r1,r2)
  }
  post.mGARCH <- function(theta, data, h1, log = TRUE){
    if (is.vector(theta))
      theta <- matrix(theta, nrow = 1)
    omega <- theta[,1]
    lambda <- theta[,2]
    beta <- theta[,3]
    alpha <- theta[,4]
    p <- theta[,5]
    mu <- theta[,6]
    N <- nrow(theta)
    pos <- 2:length(data) # # observation index (removing 1st)
    prior <- prior.mGARCH(omega=omega,lambda=lambda,beta=beta,alpha=alpha,
      p=p,mu=mu)
    d <- rep.int(-500000,N)#fixme
    for (i in 1:N){
      if (prior[i,1] == TRUE){
        h <- c(h1, omega[i] + alpha[i] * (data[pos-1]-mu[i])^2)
        for (j in pos){
          h[j] <- h[j] + beta[i] * h[j-1]
        }
        sigma <- 1 / (p[i] + ((1-p[i]) / lambda[i]))
        tmp1 <- dnorm(data[pos],mu[i],sqrt(h[pos]*sigma),log=T)
        tmp2 <- dnorm(data[pos],mu[i],sqrt(h[pos]*sigma/lambda[i]),log=T)     
        tmp <- log(p[i] * exp(tmp1) + (1-p[i]) * exp(tmp2))
        d[i] <- sum(tmp) + prior[i,2]
      }
    }
    if (!log) 
      d <- exp(d)
    as.numeric(d)
  }
  # define data subsample
  y.ss <- y[1:626] 
  # initial data variance
  h1   <- var(y) # initial variance
  N <- 1e3 # number of draws for predictive likelihood
  mu0 <- c(0.08, 0.37, 0.86, 0.03, 0.82, 0.03)  # initial parameters for MitISEM
  names(mu0) <- c("omega","lambda","beta","alpha","p","mu")
  set.seed(1234)
  cat("starting training subsample estimation", fill=TRUE) 
  mit.ss <- MitISEM(KERNEL = post.mGARCH, mu0 = mu0, data = y.ss, h1 = h1, 
    control=list(trace=TRUE))$mit
  cat("starting full sample estimation", fill=TRUE) 
  mit.fs <- MitISEM(KERNEL = post.mGARCH, mu0 = mu0, data = y,    h1 = h1, 
    control=list(trace=TRUE))$mit
  cat("starting predictive likelihood calculation", fill=TRUE) 
  N <- 1000  # number of simulations for IS 
  rep <- 50  # times to replicate application
  set.seed(1111)
  Mcompare.MitISEM <- PredLik(N,mit.fs,mit.ss,post.mGARCH,y,y.ss,h1=h1)
  # REPLICATE PRED LIKELIHOOD CALCULATION SEVERAL TIMES
  for(i in 2:rep){
    tmp <- PredLik(N,mit.fs,mit.ss,post.mGARCH,y,y.ss,h1=h1)
    Mcompare.MitISEM=mapply(rbind,Mcompare.MitISEM,tmp,SIMPLIFY=FALSE)
    if(i%%4 == 0)
      cat("rep MitISEM",i,fill=TRUE);
  }
  # REPORT MEAN AND STANDARD DEVIATION
  Means.MitISEM <- mapply(colMeans,Mcompare.MitISEM,SIMPLIFY=FALSE)
  scales <- rep(0,2)
  tmp <- Means.MitISEM[[1]]
  while(floor(tmp)==0){
    scales[i] = scales[i]+1
    tmp = tmp * 10
  }
  # average predictive likelihood and NSE from 50 repetitions
  Adj.Mcompare.MitISEM = Mcompare.MitISEM
  NSE.MitISEM <- sqrt(apply(Adj.Mcompare.MitISEM[[1]],2,var)/rep)
  table1 <- c(colMeans(Adj.Mcompare.MitISEM$PL),NSE.MitISEM)
  table1 =  rbind(rep(Adj.Mcompare.MitISEM$scale[1],2),table1)
  rownames(table1) = c("scale (10^scale)","value")
  colnames(table1) = c("Pred Lik","NSE")
  cat("Pred. Likelihood and NSE values are multiplied by 10^(scale)", fill = TRUE)
  print(round(table1,4))
  cat("number of student t components for full sample and training sample estimation",
    fill = TRUE)
  table2 <- cbind(length(mit.ss$p), length(mit.fs$p))
  colnames(table2) <- c("full sample", "training sample")
  print(round(table2,0))
}
}
