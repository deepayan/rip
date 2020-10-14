\name{rip.deconv}
\alias{rip.deconv}
\alias{rip.denoise}
\alias{rip.upscale}
\title{Non-blind deconvolution and super-resolution}
\description{
  Single-image non-blind deconvolution and super-resolution (upscaling)
  using Gaussian and hyper-Lapacian image gradient priors.
}
\usage{
rip.deconv(y, k, method = c("iterative", "direct"),
           patch = switch(method, direct = 100, iterative = 300),
           overlap = 20, ..., super.factor = 1)
}
\arguments{
  
  \item{y}{The blurred and / or downscaled image to be reconstructed, a
    \code{"rip"} object or something that can be coerced into one. Color
    images are supported (as three-channel) objects, and are split into
    individual channels before processing. Missing values (NA) are
    allowed, but only for \code{method = "direct"}.}

  \item{k}{The known blur kernel, a single-channel \code{"rip"} object
    or matrix. Mandatory for \code{rip.deconv}, but optional for
    \code{rip.denoise} and \code{rip.upscale}. See details.}

  \item{method}{Character string giving the method to be used for
    solving the linear equations derived from (weighted) least squares
    problems. The \code{"iterative"} method uses a memory-efficient
    conjugate-gradient approach through \code{\link[stats]{optim}} to
    obtain an approximate solution, using functions in the OpenCV
    library for efficient convolution and filtering operations. The
    \code{"direct"} method uses the \CRANpkg{Matrix} package to
    explicitly construct the coefficient matrix as a sparse matrix and
    performs its sparse Cholesky decomposition to obtain exact solutions
    of the linear equations. Depending on the sparsity of \code{k}, this
    may require large amounts of memory. Missing values in \code{y} are
    only supported for the \code{"direct"} method.}

  \item{patch}{Integer vector, replicated to be of length 2, giving the
    size of patches into which the input image is split to avoid
    processing large images.}

  \item{overlap}{Integer vector, replicated to be of length 2, giving
    the amount of overlap between contiguous images.}

  \item{...}{Further arguments, to be passed on to the underlying
    (unexported) workhorse functions. Supported arguments are:}

  \item{super.factor}{Integer, factor by which input image is to be
    upscaled.}
  
}

\value{
  A \code{"rip"} object representing the estimate of the latent image.
}

\details{

  These are the main functions of this package, supporting single-image
  non-blind (known blur kernel) deconvolution and super-resolution using
  a Bayesian regularization approach. As blurring and downsampling are
  linear operations, the blurred image y can be expressed in terms of
  the latent image \code{x} as \code{y = T %*% x + error}, where
  \code{T} is a sparse matrix determined by \code{k}. However, the
  problem is generally under-constrained, so requires regularization to
  solve. These funtions assume the following prior on \code{x}:
  horizontal and vertical lag-1 gradients are independent, marginally
  either mean-zero Gaussian or hyper-Laplacian (density proportional to
  \code{exp(abs(x)^alpha}}), and possibly correlated according to a
  simple 2-D lag-1 auto-regressive model.

  The estimate of \code{x} is the posterior mode, equivalent to a
  generalized form of Ridge-type regularization in the Gaussian
  case. The general case is solved using IRLS. The error term is
  generally assumed to be Gaussian, but robust loss functions can also
  be used via IRLS.

  Estimation involves solving a sparse linear problem. The size of the
  problem is large, and is solved (usually after dividing large images
  into smaller patches to keep individual problems manageable) either
  using an approximate conjugate-gradient method or a direct sparse
  Cholesky decomposition. The latter approach is potentially memory
  intensive (depending on the sparsity of \code{k}) and should be used
  with caution when \code{k} is being estimated.

  Missing values in \code{y} are allowed, but only with the
  \code{"direct"} method. This allows imputation or inpainting, with
  missing pixels restored by their posterior modes. However,

  Although \code{k} is assumed to be known, this is rarely true in
  practice. Estimating complicated blur kernels is a difficult problem
  not yet handled by this package, but external methods are available. A
  simple Fourier-domain method, assuming a symmetric kernel and Gaussian
  gradient prior, is implemented in \code{\link{symmetric.blur}}, and
  used when \code{kbw < 1}. Kernels in parametric form can be generated
  using \code{\link{make.kernel}}.

}
