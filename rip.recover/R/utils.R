

## FIXME: fix bug in R nextn() and use it

## internal function: smallest integer larger than m which has no
## prime factors larger than 5.

optDFTlength <- function(m)
    ## can be a vector
{
    hasBadFactor <- function(n)
    {
        while (n %% 2 == 0) n <- n %/% 2
        while (n %% 3 == 0) n <- n %/% 3
        while (n %% 5 == 0) n <- n %/% 5
        n > 1
    }
    getOptLen <- function(n) {
        while (hasBadFactor(n)) n <- n + 1L
        n
    }
    sapply(m, getOptLen)
}

rip.grad <- list(x = as.rip(rbind(c(-1, 1))), y = as.rip(cbind(c(-1, 1))))


## attempt to approximate MATLAB's edgetaper

rip.edgetaper <- function(x, sd = 5, boundary = 3 * sd, ...)
{
    ## gaussian kernel with 3 * 2 * sd
    gk <- rip.cv$filter$getGaussianKernel(6 * sd, sd, 5)
    gk <- as.rip(gk %*% t(gk))
    gk[] <- gk / sum(gk)
    b <- rip.filter(x, gk, ...)
    boundary <- round(boundary)
    if (2 * boundary > min(dim(x))) stop("boundary is too large for given image")
    line.begin <- seq(0, 1, length.out = round(boundary))
    line.end <- rev(line.begin)
    alpha.row <- c(line.begin, rep(1, nrow(x) - 2 * boundary), line.end)
    alpha.col <- c(line.begin, rep(1, ncol(x) - 2 * boundary), line.end)
    alpha <- outer(alpha.row, alpha.col)
    x * alpha + b * (1-alpha)
}


##' Constructs a two-dimensional blur kernel from a set of
##' predetermined parametric family.
##'
##' The parametric families supported are the same ones as in
##' \code{\link{density}}, but the interpretation of the bandwidth is
##' different; it gives the distance from the origin outside which the
##' kernel is 0, except for the Gaussian kernel where it is twice the
##' variance parameter.
##'
##' The center pixel of the kernel spans the [-0.5, 0.5]^2 square, and
##' other pixels are shifted accordingly. To compute the kernel value
##' on a given pixel, the kernel function is evaluated on a roughly
##' 10x10 grid within the pixel and then averaged. Bandwidth values
##' less than 0.5 essentially produce degenerate kernels.
##' 
##' @title Construct a two-dimensional blur kernel
##' @param h Numeric scalar giving bandwidth (in pixels).
##' @param kern Character string giving parametric family.
##' @param dim Integer scalar giving kernel dimension; must be 1 or 2.
##' @param normalize Logical flag indicating whether the kernel should
##'     be normalized to have total sum 1.
##' @return A \code{"rip"} object giving the kernel, or, if \code{dim
##'     = 1}, a vector.
make.kernel <-
    function(h = 1,
             kern = c("epanechnikov","rectangular", "triangular",
                      "biweight", "gaussian", "cosine", "optcosine"),
             dim = 2, normalize = TRUE)
{
    if (h <= 0) stop("Bandwidth 'h' must be positive")
    ## only limit to standard 'finite' kernels. For Gaussian, sd is taken to be h/2.
    ## consider adding cauchy (with parameter angle) to model out-of-focus blur.
    kern <- match.arg(kern)
    klim <- ceiling(h-0.5) + 0.5 # kernel goes from -[klim] to [klim]
    ## evaluate on finer grid, then 'integrate' between +/- 0.5
    x <- seq(-klim, klim, length.out = 10 * (2 * klim))
    ax <- abs(x)
    k <- switch(kern,
                gaussian = dnorm(x, sd = h / 2),
                rectangular = ifelse(ax < h, 0.5/h, 0),
                triangular = ifelse(ax < h, (1 - ax/h)/h, 0),
                epanechnikov = ifelse(ax < h, 3/4 * (1 - (ax/h)^2)/h, 0),
                biweight = ifelse(ax < h, 15/16 * (1 - (ax/h)^2)^2/h, 0),
                cosine = ifelse(ax < h, (1 + cos(pi * x/h))/(2 * h), 0),
                optcosine = ifelse(ax < h, pi/4 * cos(pi * x/(2 * h))/h, 0))
    if (!(dim %in% c(1, 2))) stop("Only 'dim=1' and 'dim=2' are allowed.")
    if (dim == 1)
    {
        f <- cut(x, seq(-klim, klim, by = 1), include.lowest = TRUE)
        kk <- return(tapply(k, f, mean))
        if (normalize) kk[] <- kk / sum(kk)
        kk
    }
    else if (dim == 2)
    {
        ## could just take product measure, but will do slightly better
        f <- cut(x, seq(-klim, klim, by = 1), include.lowest = TRUE)
        ## str(list(f = as.character(f), x = x, k = k))
        g <- expand.grid(f1 = f, f2 = f, KEEP.OUT.ATTRS=FALSE)
        g$k <- as.vector(outer(k, k))
        kk <- xtabs(k ~ f1 + f2, g) / xtabs( ~ f1 + f2, g) # mean
        if (normalize) kk[] <- kk / sum(kk)
        as.rip(kk)
    }
}

##' Set small values in a numeric vector to zero, similar to
##' \code{\link{zapsmall}}, except that non-zero elements are retained
##' as is instead of being rounded.
##'
##' @title Set small values to zero
##' @param x A numeric vector, matrix, or array.
##' @param digits Numeric giving "number" of decimal digits to retain
##' @param prop Proportion (fraction) of maximum below which elements
##'     are to be set to zero. Overrides \code{digits}.
##' @param threshold Threshold absolute value below which elements are
##'     to be set to zero. Overrides \code{prop}.
##' @return Modified version of \code{x}, retaining attributes, with
##'     small elements set to zero.
zapsmallp <- function(x, digits = 2, prop = 10^(-digits), threshold = max(abs(x)) * prop)
{
    x[abs(x) < threshold] <- 0
    x
}

