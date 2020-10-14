\name{make.kernel}
\alias{make.kernel}
\title{Construct a two-dimensional blur kernel}
\description{
  Constructs a two-dimensional blur kernel from a set of
  predetermined parametric family.
}
\usage{
make.kernel(h = 1,
            kern = c("epanechnikov", "rectangular", "triangular", "biweight",
                     "gaussian", "cosine", "optcosine"),
            dim = 2, normalize = TRUE)
}
\arguments{

  \item{h}{Numeric scalar giving bandwidth (in pixels).}

  \item{kern}{Character string giving parametric family.}

  \item{dim}{Integer scalar giving kernel dimension; must be 1 or 2.}

  \item{normalize}{Logical flag indicating whether the kernel should be
    normalized to have total sum 1.}

}
\value{

  A \code{"rip"} object giving the kernel, or, if \code{dim = 1}, a
  vector.

}

\details{

  The parametric families supported are the same ones as in
  \code{\link{density}}, but the interpretation of the bandwidth is
  different; it gives the distance from the origin outside which the
  kernel is 0, except for the Gaussian kernel where it is twice the
  variance parameter.

  The center pixel of the kernel spans the [-0.5, 0.5]^2 square, and
  other pixels are shifted accordingly. To compute the kernel value
  on a given pixel, the kernel function is evaluated on a roughly
  10x10 grid within the pixel and then averaged. Bandwidth values
  less than 0.5 essentially produce degenerate kernels.

}

\examples{

image(make.kernel(h = 10, kern = "gaussian"))
image(make.kernel(h = 10, kern = "triangular"))
  
}

