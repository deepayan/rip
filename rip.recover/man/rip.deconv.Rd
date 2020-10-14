\name{rip.deconv}
\alias{rip.deconv}
\alias{rip.denoise}
\alias{rip.upscale}
\title{Non-blind deconvolution and super-resolution}
\description{
  Single-image non-blind deconvolution and super-resolution (upscaling)
  using Gaussian and hyper-Lapacian image gradient priors. Supports
  missing values (using the direct method only), as well as robust loss
  functions instead of squared error loss for image fidelity.
}
\usage{
rip.deconv(y, k, method = c("iterative", "direct"),
           patch = switch(method, direct = 100, iterative = 300),
           overlap = 20, ..., super.factor = 1)

rip.upscale(y, factor = 2,
            ksigma = 1, kbw = ksigma, k = NULL,
            method = c("iterative", "direct"),
            patch = round(switch(method, direct = 100, iterative = 300) / factor),
            overlap = round(20 / factor), ...)

rip.denoise(y, ksigma = 1, kbw = ksigma, k = NULL,
            method = c("iterative", "direct"),
            patch = switch(method, direct = 100, iterative = 300),
            overlap = 20, ...)

}
\arguments{
  
  \item{y}{The blurred and / or downscaled image to be reconstructed, a
    \code{"rip"} object or something that can be coerced into one. Color
    images are supported (as three-channel) objects, and are split into
    individual channels before processing. Missing values (NA) are
    allowed, but only for \code{method = "direct"}.}

  \item{k}{The known blur kernel applied to the latent image, a
    single-channel \code{"rip"} object or matrix. Mandatory for
    \code{rip.deconv}, but optional for \code{rip.denoise} and
    \code{rip.upscale}. See details.}

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
    (unexported) workhorse functions. Useful arguments are:
    
    \describe{

      \item{\code{cache.dir}:}{ The name of a directory to store cached
	sparse matrices. Used to store downscaling operators if
	specified; by default they are computed as necessary (the
	penalty is relatively small). }

      \item{\code{yerror}:}{ Character string giving loss function used
	to measure image fidelity. Possible values are \code{"normal"}
	for squared error loss, \code{"huber"} for Huber loss,
	\code{"bisquare"} for Tukey bisquare loss, and \code{"poisson"}
	for squared error loss with variance proportional to mean. The
	last option is not really useful, as there is an arbitrary scale
	factor has no natural estimate. }

      \item{\code{huber.k = 1.345}:}{ Parameter for Huber loss. The
	scale parameter is estimated using the MAD in each iteration. }

      \item{\code{bisquare.c = 4.685}:}{ Parameter for bisquare
	loss. The scale parameter is estimated using the MAD in each
	iteration.  }

      \item{\code{wt.thres = 0.01}:}{ IRLS weights are restricted to be
	no smaller than this, both for image error and hyper-Lapacian
	prior weights.}

      \item{\code{niter.irls = 5}:}{ Number of IRLS iterations.  }

      \item{\code{verbose = FALSE}:}{ Logical flag; whether progress
	should be indicated. Useful for long-running computations. }

      \item{\code{cg.update}:}{ Character string giving type of
	conjugate-gradient update (passed on to \code{optim}). Possible
	values are \code{"FR"} for Fletcher-Reeves (the default),
	\code{"PR"} for Polak-Ribiere, and \code{"BS"} for
	Beale-Sorenson. Polak-Ribiere is slightly slower than
	Fletcher-Reeves but sometimes gives better results. }

      \item{\code{maxit}:}{ The maximum number of conjugate-gradient
	iterations in \code{optim}. In fact, other components of the
	\code{control} argument of \code{optim} can also be specified. }

    }
  }

  \item{super.factor, factor}{Integer, factor by which input image is to
    be upscaled.}

  \item{kbw, ksigma}{For \code{rip.denoise} and \code{rip.upscale}, a
    numeric value giving the bandwidth (in units of pixels in the input
    image) for an Epanechnikov kernel as constructed using
    \code{\link{symmetric.blur}}. For \code{rip.upscale}, the kernel
    that applies on the latent image has bandwidth \code{kbw * factor}.
    Two special values are \code{0} to indicate no blurring and a
    negative value to estimate the kernel from the image itself assuming
    a correlated gradient prior. For the last case, explicitly using
    \code{\link{symmetric.blur}} allows more control, especially on the
    size of the kernel.}

}

\value{
  A \code{"rip"} object representing the estimate of the latent image.
}

\details{

  These are the main functions of this package, supporting single-image
  non-blind (known blur kernel) deconvolution and super-resolution using
  a Bayesian regularization approach. As blurring and downsampling are
  linear operations, the blurred image \code{y} can be expressed in
  terms of the latent image \code{x} as \code{y = T x + error}, where
  \code{T} is a sparse matrix determined by \code{k}. However, the
  problem is generally under-constrained, so requires regularization to
  solve. These funtions assume the following prior on \code{x}:
  horizontal and vertical lag-1 gradients are independent, marginally
  either mean-zero Gaussian or hyper-Laplacian (density proportional to
  \code{exp(abs(x)^alpha)}), and possibly correlated according to a
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
