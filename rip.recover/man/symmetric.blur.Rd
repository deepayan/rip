\name{symmetric.blur}
\alias{symmetric.blur}
\title{Estimate symmetric blur kernel}
\description{
  Estimate symmetric blur kernel assuming Gaussian image gradient prior
}
\usage{
symmetric.blur(y, kdim = round(dim(y)/3), resize = 1,
               g.method = c("autoreg", "nonpar", "none"), g = NULL,
               rho.along = 0.3, rho.across = 0.6, eta.sq = 0.005^2,
               corr.grad = TRUE, np.span = 2/3, trim = TRUE,
               zap.digits = 1.5, edgetaper = FALSE)
}
\arguments{

  \item{y}{Input blurred image, a \code{"rip"} object or something that
    can be coerced to one. Color images are handled by first converting
    to grayscale using \code{\link[rip.opencv]{rip.desaturate}}.}

  \item{kdim}{Extent of the kernel. The actual size of the estimated
    kernel is \code{2 * kdim + 1}, before being trimmed if requested.}

  \item{resize}{Factor by which the estimate is to be resized before
    being returned.}

  \item{g.method}{Method used to compute variance of Fourier
    coefficients. Ignored if \code{g} is explicitly specified. The
    \code{"none"} method assumes constant variance, corresponding to a
    prior that assumed uncorrelated gradients. The \code{"autoreg"}
    method uses a variance model corresponding to a prior with gradients
    having a lag-1 auto-regressive structure. The \code{"nonpar"} method
    estimates the variance by smoothing the absolute DFT coefficients;
    this is of dubious usefulness but included to allow
    experimentation.}

  \item{g}{Optional list with components named \code{h} and \code{v},
    giving the variances of DFT coefficients for horizontal and vertical
    gradients respectively. Should be of same size as \code{y}.}

  \item{rho.along}{Correlation coefficient defining auto-regressive
    gradient prior model along the direction of the gradient.}

  \item{rho.across}{Correlation coefficient defining auto-regressive
    gradient prior model across the direction of the gradient.}

  \item{eta.sq}{Variance of the error terms (after taking gradients)}

  \item{corr.grad}{Whether the error terms should be assumed to be
    correlated after taking gradients. If \code{FALSE}, the error
    gradients are assumed to be independent; if \code{TRUE}, the errors
    in \code{y} are assumed to be independent, and correlation in
    gradient errors computed accordingly. In practice, there is not much
    difference for reasonable values of \code{eta.sq}.}

  \item{np.span}{Passed to \code{g.nonpar} when \code{g.method="nonpar"}.}

  \item{trim}{Logical flag; whether the estimated kernel should be
    "trimmed" by setting small elements to zero and removing all-zero
    boundary rows.}
  
  \item{zap.digits}{The \code{digits} argument to
    \code{\link{zapsmallp}} to use when \code{trim=TRUE}.}

  \item{edgetaper}{Logical, whether \code{\link{rip.edgetaper}} should
    be applied on the gradients before further processing.}

}

\value{
  A \code{"rip"} object giving the estimated blur kernel.
}

\details{ Estimates the blur kernel from a blurred image assuming that
  the kernel is symmetric and that image gradients follow a Gaussian
  prior.
}
