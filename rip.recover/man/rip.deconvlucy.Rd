\name{rip.deconvlucy}
\alias{rip.deconvlucy}
\title{Richardson-Lucy Algorithm}
\description{
Image deconvoltion using the Richardson-Lucy algorithm
}
\usage{
rip.deconvlucy(y, k, niter = 10, x.init = NULL, keep.intermediate = FALSE)
}
\arguments{
  
  \item{y}{The blurred image, a \code{"rip"} object or something that
    can be coerced into one. Multi-channel (color) images are
    supported.}

  \item{k}{The known blur kernel, a single-channel \code{"rip"} object
    or matrix.}

  \item{niter}{The number of Richardson-Lucy (EM) iterations to
    perform.}

  \item{x.init}{Initial estimate of the latent image, in the same form
    as \code{y}. There is usually no benefit to providing this (by
    default the initial estimate is flat image filled with the average
    of \code{y}), and the option is only provided for experimentation.}

  \item{keep.intermediate}{A logical flag to indicate whether
    intermediate estimates should be saved. This is useful to track
    performance; increasing the number of iterations does not always
    improve quality.}

}

\value{
  A \code{"rip"} object representing the estimate of the latent
  image. If \code{keep.intermediate = TRUE}, the return value
  additionally has an attribute named \code{"history"}, which is a list
  containing the result of each iteration.
}

\details{
  The Richardson-Lucy algorithm is used to recover the underlying latent
  image given a blurred image and the blur kernel as input, using a
  Poisson noise model. The algorithm here uses an efficient
  implementation of correlation filters from the OpenCV library. To
  minimize boundary artifacts, the REPLICATE boundary type is used for
  filtering.
}
