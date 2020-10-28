

## Richardson-Lucy: Add borderType option for filter()?: REPLICATE
## (default) is probably best. Check if WRAP works

##' Image deconvoltion using the Richardson-Lucy algorithm
##'
##' The Richardson-Lucy algorithm is used to recover the underlying
##' latent image given a blurred image and the blur kernel as input,
##' using a Poisson noise model. The algorithm here uses an efficient
##' implementation of correlation filters from the OpenCV library. To
##' minimize boundary artifacts, the REPLICATE boundary type is used
##' for filtering.
##'
##' @title Richardson-Lucy Algorithm
##' @param y The blurred image, a \code{"rip"} object or something
##'     that can be coerced into one. Multi-channel (color) images are
##'     supported.
##' @param k The known blur kernel, a single-channel \code{"rip"}
##'     object or matrix.
##' @param niter The number of Richardson-Lucy (EM) iterations to
##'     perform.
##' @param x.init Initial estimate of the latent image, in the same
##'     form as \code{y}. There is usually no benefit to providing
##'     this (by default the initial estimate is flat image filled
##'     with the average of \code{y}), and the option is only provided
##'     for experimentation.
##' @param keep.intermediate A logical flag to indicate whether
##'     intermediate estimates should be saved. This is useful to
##'     track performance; increasing the number of iterations does
##'     not always improve quality.
##' @return A \code{"rip"} object representing the estimate of the
##'     latent image. If \code{keep.intermediate = TRUE}, the return
##'     value additionally has an attribute named \code{"history"},
##'     which is a list containing the result of each iteration.
##' 
rip.deconvlucy <-
    function(y, k, niter = 10,
             x.init = NULL, keep.intermediate = FALSE)
        ## FIXME: add option to split and merge, but be careful about multi-channel input
{
    y <- as.rip(y)
    k <- as.rip(k)
    if (nchannel(k) != 1) stop("'k' must have exactly one channel")
    if (keep.intermediate)
    {
        h <- vector(mode = "list", length = niter)
    }
    x <- y
    x[] <- if (is.null(x.init)) mean(y) else x.init
    .k <- rip.flip(k)
    for (i in seq_len(niter))
    {
        ## FIXME: check if valid / full convolution would be more correct
        r <- y / rip.filter(x, .k)
        r[!is.finite(r)] <- 0 # FIXME: is this the right thing to do?
        x[] <- x * rip.filter(r, k)
        if (keep.intermediate) h[[i]] <- x
    }
    if (keep.intermediate) attr(x, "history") <- h
    x
}

rip.deconv.rl <- rip.deconvlucy # older name for compatibility, remove eventually



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Wiener filter for (non-blind) image deconvolution
##' @param y blurred image
##' @param k blur kernel of psf
##' @param S Estimate of the power spectra of the original image
##' @param N Estimate of the power spectra of the noise
##' @return 
##' @author kaustav nandy
rip.deconvwnr <- function(y, k, S, N)
{
    dim.y <- dim(y)
    K <- rip.dft(y, pad = dim(y))
    Y <- rip.dft(y)
    G <- K / (Mod(K)^2 + S / N)
    rip.dft(Y * G, inverse = TRUE)
}

## deconvolution / super-resolution using hyper-Laplacian and possibly
## correlated gradient prior. Workhorse functions doing the actual
## implementation are in sparse-matrix.R and conjugate-gradient.R

## Separate interfaces for denoising, deconvolution and
## super-resolution.

## Overlap should be at least kernel size for deconvolution, but 20
## should be OK 


nonblind.hyperlap.direct <-
    function(y, k, alpha, lambda, rho,
             patch, overlap,
             cache.dir = getOption("rip.cache.dir"),
             ...,
             super.factor = 1)
{
    force.full.rank <- FALSE
    rank.offset <- 1 - force.full.rank
    rank.convtype <- if (force.full.rank) "same" else "valid"
    ## message("Using cache: ", cache.dir)
    patch <- rep(patch, length.out = 2) # needed as arguments to linear2sparse below
    overlap <- rep(overlap, length.out = 2)
    dlatent <- super.factor * patch + dim(k) - 1 # size of latent image
    ## Note: if cache.dir is non-NULL, a md5-based identifier is used in conv2sparse
    Tk <- conv2sparse(k, dlatent[1], dlatent[2], cache.dir = cache.dir)
    if (super.factor < 1 || round(super.factor) != super.factor)
        stop("'super.factor' must be a positive integer")
    if (super.factor > 1)
    {
        ## str(list(Tk = dim(Tk), patch = patch, super.factor =
        ## super.factor)) kernel must be modified to incorporate
        ## subsampling. Where to start is arbitrary, we choose to set
        ## the first element to TRUE, because that's what we do in the
        ## iterative case to resolve ambiguity in defining Tk'.
        subsamp <- rep(FALSE, super.factor)
        ## subsamp[ceiling(0.5 * super.factor)] <- TRUE
        subsamp[1] <- TRUE
        Tsubsamp <- linear2sparse(FUN = function(x) x[subsamp, subsamp],
                                  super.factor * patch[1], super.factor * patch[2],
                                  cache = sprintf("subsample-%g", super.factor),
                                  cache.dir = cache.dir)
        Tk <- Tsubsamp %*% Tk
    }
    ## IID gradients (this reduces individual rank, but sum seems to
    ## be full rank. Any guarantee? Holds for AR adjustment too?
    ## Otherwise use 'same' colvolution --- shouldn't matter much)
    Td.h <- conv2sparse(rip.flip(rip.grad$x), dlatent[1], dlatent[2],
                        cache = "gradh-iid", cache.dir = cache.dir,
                        conv.type = rank.convtype)
    Td.v <- conv2sparse(rip.flip(rip.grad$y), dlatent[1], dlatent[2],
                        cache = "gradv-iid", cache.dir = cache.dir,
                        conv.type = rank.convtype)
    ## Update with AR model if applicable
    if (rho$across != 0 || rho$along != 0)
    {
        ## older DFT transform version
        ## Td.h <- ar2sparse(c(rho$across, rho$along), dlatent[1], dlatent[2]-1,
        ##                   match.sd = FALSE, ## rip.flip(rip.grad$x), # pre-filter
        ##                   cache = sprintf("gradh-%g-%g", rho$along, rho$across),
        ##                   cache.dir = cache.dir,
        ##                   zap.digits = 3) %*% Td.h
        ## Td.v <- ar2sparse(c(rho$along, rho$across), dlatent[1]-1, dlatent[2],
        ##                   match.sd = FALSE, ## rip.flip(rip.grad$y), # pre-filter
        ##                   cache = sprintf("gradv-%g-%g", rho$along, rho$across),
        ##                   cache.dir = cache.dir,
        ##                   zap.digits = 3) %*% Td.v
        ## newer convolution version
        hconv <- with(rho, as.rip(matrix(c(along * across, -along, -across, 1), 2, 2)))
        vconv <- with(rho, as.rip(matrix(c(along * across, -across, -along, 1), 2, 2)))
        hconv[] <- hconv / sqrt(sum(hconv^2)) # normalize variance
        vconv[] <- vconv / sqrt(sum(vconv^2)) # normalize variance
        Td.h <- conv2sparse(rip.flip(hconv), dlatent[1], dlatent[2] - rank.offset,
                            cache = sprintf("arh-%g-%g", rho$along, rho$across),
                            cache.dir = cache.dir,
                            conv.type = "same") %*% Td.h
        Td.v <- conv2sparse(rip.flip(vconv), dlatent[1] - rank.offset, dlatent[2],
                            cache = sprintf("arv-%g-%g", rho$along, rho$across),
                            cache.dir = cache.dir,
                            conv.type = "same") %*% Td.v
    }
    direct.irls(Tk, y, lambda = lambda, alpha = alpha, Td.h = Td.h, Td.v = Td.v,
                super.factor = super.factor, patch = patch, overlap = overlap,
                latent.dim = dlatent, ...)
}

nonblind.hyperlap.iterative <-
    function(y, k, alpha, lambda, rho,
             patch, overlap,
             cache.dir = getOption("rip.cache.dir"),
             decorrelate = c("frequency", "image"),  # frequency is slower but does not cause shifts
             cluster = NULL,
             ...,
             super.factor = 1)
{
    decorrelate <- match.arg(decorrelate)
    ## message("Using cache: ", cache.dir)
    patch <- rep(patch, length.out = 2) # needed as arguments to linear2sparse below
    overlap <- rep(overlap, length.out = 2)
    dlatent <- super.factor * patch + dim(k) - 1 # size of latent image
    ## Note: if cache.dir is non-NULL, a md5-based identifier is used in conv2sparse
    if (super.factor < 1 || round(super.factor) != super.factor)
        stop("'super.factor' must be a positive integer")
    ## Update with AR model if applicable
    if (rho$across != 0 || rho$along != 0)
    {
        if (decorrelate == "frequency")
        {
            g <- g.autoreg(d = dlatent, rho = rho)
            g.h <- g$h
            g.v <- g$v
            ## g.h <- as.rip(outer(var.AR1(rho$across, dlatent[1]),
            ##                     var.AR1(rho$along, dlatent[2]-1)))
            ## g.v <- as.rip(outer(var.AR1(rho$along, dlatent[1]-1),
            ##                     var.AR1(rho$across, dlatent[2])))
        }
        else # (decorrelate == "image")
        {
            g.h <- with(rho, as.rip(matrix(c(along * across, -along, -across, 1), 2, 2)))
            g.v <- with(rho, as.rip(matrix(c(along * across, -across, -along, 1), 2, 2)))
            ## flip for convolution instead of filtering, and normalize variance
            g.h[] <- rip.flip(g.h) / sqrt(sum(g.h^2))
            g.v[] <- rip.flip(g.v) / sqrt(sum(g.v^2))
        }
    }
    else g.h <- g.v <- 1
    if (is.null(cluster))
        iterative.irls(k, y, lambda = lambda, alpha = alpha, g.h = g.h, g.v = g.v,
                       super.factor = super.factor, patch = patch, overlap = overlap,
                       latent.dim = dlatent, ...)
    else
        iterative.irls.parallel(k, y, lambda = lambda, alpha = alpha, g.h = g.h, g.v = g.v,
                       super.factor = super.factor, patch = patch, overlap = overlap,
                       latent.dim = dlatent, ..., cluster = cluster)
}


ktrim0odd <- function(k)
{
    ## remove any all-zero boundary rows / columns, and then add some
    ## (back) to ensure that kernel size is odd
    while (rowSums(k)[1] == 0) { k <- k[-1, , drop = FALSE] }
    while (rowSums(k)[nrow(k)] == 0) { k <- k[-nrow(k), , drop = FALSE] }
    while (colSums(k)[1] == 0) { k <- k[, -1, drop = FALSE] }
    while (colSums(k)[ncol(k)] == 0) { k <- k[, -ncol(k), drop = FALSE] }
    if (nrow(k) %% 2 == 0) k <- rbind(k, 0)
    if (ncol(k) %% 2 == 0) k <- cbind(k, 0)
    as.rip(k)
}


rip.deconv.channel <- function(y, k, FUN, alpha, lambda, rho,
                               patch, overlap, ..., super.factor = 1)
{
    FUN(y, k, alpha = alpha,
        lambda = lambda, rho = rho,
        patch = patch, overlap = overlap,
        ..., super.factor = super.factor)

    
}

## NOTE: patch decreased to maximum of image size if larger


##' Single-image non-blind deconvolution and super-resolution
##' (upscaling) using Gaussian and hyper-Lapacian image gradient
##' priors. Supports missing values (using the direct method only), as
##' well as robust loss functions instead of squared error loss for
##' image fidelity.
##'
##' These are the main functions of this package, supporting
##' single-image non-blind (known blur kernel) deconvolution and
##' super-resolution using a Bayesian regularization
##' approach. \code{rip.upscale} and \code{rip.denoise} are simple
##' wrappers for \code{rip.deconv} with additional ways to specify the
##' blur kernel, and pass along additional arguments to it.
##'
##' As blurring and downsampling are linear operations, the blurred
##' image \code{y} can be expressed in terms of the latent image
##' \code{x} as \code{y = T x + error}, where \code{T} is a sparse
##' matrix determined by \code{k}. However, the problem is generally
##' under-constrained, so requires regularization to solve. These
##' funtions assume the following prior on \code{x}: horizontal and
##' vertical lag-1 gradients are independent, marginally either
##' mean-zero Gaussian or hyper-Laplacian (density proportional to
##' \code{exp(abs(x)^alpha)}), and possibly correlated according to a
##' simple 2-D lag-1 auto-regressive model.
##'
##' The estimate of \code{x} is the posterior mode, equivalent to a
##' generalized form of Ridge-type regularization in the Gaussian
##' case. The general case is solved using IRLS. The error term is
##' generally assumed to be Gaussian, but robust loss functions can
##' also be used via IRLS.
##'
##' Estimation involves solving a sparse linear problem. The size of
##' the problem is large, and is solved (usually after dividing large
##' images into smaller patches to keep individual problems
##' manageable) either using an approximate conjugate-gradient method
##' or a direct sparse Cholesky decomposition. The latter approach is
##' potentially memory intensive (depending on the sparsity of
##' \code{k}) and should be used with caution when \code{k} is being
##' estimated.
##'
##' Missing values in \code{y} are allowed, but only with the
##' \code{"direct"} method. This allows imputation or inpainting, with
##' missing pixels restored by their posterior modes. However, 
##'
##' Although \code{k} is assumed to be known, this is rarely true in
##' practice. Estimating complicated blur kernels is a difficult
##' problem not yet handled by this package, but external methods are
##' available.
##'
##' \code{rip.denoise} and \code{rip.upscale} are expected to work
##' with images that have only mild burring, usually due to slight
##' lack of focus or limitations of camera optics. A simple
##' Fourier-domain method, assuming a symmetric kernel and Gaussian
##' gradient prior, is implemented in \code{\link{symmetric.blur}},
##' and used when \code{kbw < 1} for \code{rip.denoise} and
##' \code{rip.upscale}. Kernels in parametric form can be generated
##' using \code{\link{make.kernel}}.
##' 
##' @title Non-blind deconvolution and super-resolution
##' @param y The blurred and / or downscaled image to be
##'     reconstructed, a \code{"rip"} object or something that can be
##'     coerced into one. Color images are supported (as
##'     three-channel) objects, and are split into individual channels
##'     before processing. Missing values (NA) are allowed, but only
##'     for \code{method = "direct"}.
##' @param k The known blur kernel applied to the latent image, a
##'     single-channel \code{"rip"} object or matrix. Mandatory for
##'     \code{rip.deconv}, but optional for \code{rip.denoise} and
##'     \code{rip.upscale}. See details.
##' @param method Character string giving the method to be used for
##'     solving the linear equations derived from (weighted) least
##'     squares problems. The \code{"iterative"} method uses a
##'     memory-efficient conjugate-gradient approach through
##'     \code{\link[stats]{optim}} to obtain an approximate solution,
##'     using functions in the OpenCV library for efficient
##'     convolution and filtering operations. The \code{"direct"}
##'     method uses the \CRANpkg{Matrix} package to explicitly
##'     construct the coefficient matrix as a sparse matrix and
##'     performs its sparse Cholesky decomposition to obtain exact
##'     solutions of the linear equations. Depending on the sparsity
##'     of \code{k}, this may require large amounts of memory. Missing
##'     values in \code{y} are only supported for the \code{"direct"}
##'     method.
##' @param alpha The hyper-Laplacian parameter. \code{alpha = 2}
##'     corresponds to Gaussian, and \code{alpha = 1} to double
##'     exponential or Laplacian. \code{alpha = 0.8} is commonly used
##'     to model natural images, and often referred to as a "sparse"
##'     prior because it puts relatively more weight to 0 gradients.
##' @param lambda The regularization or tuning parameter.
##' @param rho List with components \code{along} and \code{across}
##'     specifying local dependence of gradients as a 2-D
##'     auto-regressive process in terms of correlation along and
##'     across the gradient direction. \code{along = 0.3} and
##'     \code{across = 0.6} are good values to try.
##' @param patch Integer vector, replicated to be of length 2, giving
##'     the size of patches into which the input image is split to
##'     avoid processing large images.
##' @param overlap Integer vector, replicated to be of length 2,
##'     giving the amount of overlap between contiguous images.
##' @param ... Further arguments, to be passed on to the underlying
##'     (unexported) workhorse functions. Useful arguments are:
##'
##'     \describe{
##' 
##'     \item{\code{cache.dir}:}{ The name of a directory to store
##'       cached sparse matrices. Used to store downscaling operators
##'       if specified; by default they are computed as necessary (the
##'       penalty is relatively small). }
##'
##'     \item{\code{yerror}:}{ Character string giving loss function
##'       used to measure image fidelity. Possible values are
##'       \code{"normal"} for squared error loss, \code{"huber"} for
##'       Huber loss, \code{"bisquare"} for Tukey bisquare loss, and
##'       \code{"poisson"} for squared error loss with variance
##'       proportional to mean. The last option is not really useful,
##'       as there is an arbitrary scale factor has no natural
##'       estimate. }
##'
##'     \item{\code{huber.k = 1.345}:}{ Parameter for Huber loss. The
##'       scale parameter is estimated using the MAD in each
##'       iteration. }
##' 
##'     \item{\code{bisquare.c = 4.685}:}{ Parameter for bisquare
##'       loss. The scale parameter is estimated using the MAD in each
##'       iteration.  }
##' 
##'     \item{\code{wt.thres = 0.01}:}{ IRLS weights are restricted to
##'     be no smaller than this, both for image error and
##'     hyper-Lapacian prior weights.}
##' 
##'     \item{\code{niter.irls = 5}:}{ Number of IRLS iterations.  }
##' 
##'     \item{\code{verbose = FALSE}:}{ Logical flag; whether progress
##'     should be indicated. Useful for long-running computations. }
##'
##'     \item{\code{cg.update}:}{ Character string giving type of
##'       conjugate-gradient update (passed on to
##'       \code{optim}). Possible values are \code{"FR"} for
##'       Fletcher-Reeves (the default), \code{"PR"} for
##'       Polak-Ribiere, and \code{"BS"} for
##'       Beale-Sorenson. Polak-Ribiere is slightly slower than
##'       Fletcher-Reeves but sometimes gives better results. }
##' 
##'     \item{\code{maxit}:}{ The maximum number of conjugate-gradient
##'     iterations in \code{optim}. In fact, other components of the
##'     \code{control} argument of \code{optim} can also be specified. }
##' 
##' 
##' }
##'
##' @param super.factor Integer, factor by which input image is to be
##'     upscaled.
##' @param factor Integer, factor by which input image is to be
##'     upscaled.
##' @param kbw For \code{rip.denoise} and \code{rip.upscale}, a
##'     numeric value giving the bandwidth (in units of pixels in the
##'     input image) for an Epanechnikov kernel as constructed using
##'     \code{\link{symmetric.blur}}. For \code{rip.upscale}, the
##'     kernel that applies on the latent image has bandwidth
##'     \code{kbw * factor}.  Two special values are \code{0} to
##'     indicate no blurring and a negative value to estimate the
##'     kernel from the image itself assuming a correlated gradient
##'     prior. For the last case, explicitly using
##'     \code{\link{symmetric.blur}} allows more control, especially
##'     on the size of the kernel.
##' @return A \code{"rip"} object representing the estimate of the
##'     latent image.
##' 
rip.deconv <- function(y, k, method = c("iterative", "direct"),
                       alpha = 2, lambda = 0.01,
                       rho = list(along = 0, across = 0),
                       patch = switch(method, direct = 100, iterative = 300),
                       overlap = 20, ..., super.factor = 1)
{
    stopifnot(inherits(y, "rip"))
    method <- match.arg(method)
    k <- ktrim0odd(k)
    if (any(dim(k) %% 2 == 0)) stop("kernel must have odd dimensions")
    patch <- pmin(patch, dim(y)) # patch size should not exceed size of y
    ## FIXME: TODO adjust overlap so that it is as uniform as possible
    ## without increasing number of patches (instead of only the last
    ## overlap being large, as it currently is).
    ## 
    ## Note that latent image is larger than y if kernel is more than
    ## 1x1, although by convention we estimate xhat of the same size
    ## as x. The latent image (for each patch) has size latent.dim <-
    ## patch + dim(k) - 1
    FUN <- switch(method,
                  direct = nonblind.hyperlap.direct, 
                  iterative = nonblind.hyperlap.iterative)
    nc <- nchannel(y)
    if (nc == 1) ## grayscale
        rip.deconv.channel(y, k, FUN = FUN,
                           alpha = alpha, lambda = lambda, rho = rho,
                           ..., 
                           patch = patch, overlap = overlap,
                           super.factor = super.factor, label = "[GRAY]")
    else if (nc != 3) stop("Input image 'y' must have exactly 1 or 3 channels")
    else
    {
        rgb.layers <- vector(mode = "list", length = nc)
        rgb.array <- as.array(y)
        for (i in 1:nc)
            rgb.layers[[i]] <-
                rip.deconv.channel(as.rip(rgb.array[,,i]), k, FUN = FUN,
                                   alpha = alpha, lambda = lambda, rho = rho,
                                   ..., 
                                   patch = patch, overlap = overlap,
                                   super.factor = super.factor,
                                   label = c("[Rgb]", "[rGb]", "[rgB]")[i])
        ## assume without checking that all layers have same size
        ans <- array(0, dim = c(dim(rgb.layers[[1]]), nc))
        for (i in 1:nc) ans[,,i] <- rgb.layers[[i]]
        as.rip(ans)
    }
}



rip.upscale <-
    function(y, factor = 2,
             ksigma = 1, kbw = ksigma, # for back-compatibility, remove ksigma eventually
             k = NULL, # kernel for the upscales image
             method = c("iterative", "direct"),
             patch = round(switch(method, direct = 100, iterative = 300) / factor),
             overlap = round(20 / factor), ...)
        ## add option for kernel to be uniform on factor x factor square
{
    method <- match.arg(method)
    if (is.null(k)) 
        k <- if (kbw < 0) symmetric.blur(y, kdim = pmin(round(dim(y) / 3), overlap),
                                         resize = factor)
             else if (kbw == 0) as.rip(matrix(1,factor,factor))
             else make.kernel(h = factor * kbw, kern = "epanechnikov")
    k[] <- k / sum(k)
    ## str(k)
    rip.deconv(y, k, method = method, patch = patch, overlap = overlap,
               ..., super.factor = factor)
}


rip.denoise <-
    function(y, 
             ksigma = 1, kbw = ksigma, # for back-compatibility, remove ksigma eventually
             k = NULL,
             method = c("iterative", "direct"),
             patch = switch(method, direct = 100, iterative = 300),
             overlap = 20, ...)
{
    method <- match.arg(method)
    if (is.null(k)) 
        k <- if (kbw < 0) symmetric.blur(y, kdim = pmin(round(dim(y) / 3), overlap))
             else if (kbw == 0) as.rip(matrix(1,1,1))
             else make.kernel(h = kbw, kern = "epanechnikov")
    k[] <- k / sum(k)
    ## str(dim(k))
    rip.deconv(y, k, method = method, patch = patch, overlap = overlap, ...)
}


