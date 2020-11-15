
## Assuming k is symmetric, want to estimate it crudely from an input
## image given some g (or rho) from its Mod(DFT)

## The nonpar method is not really OK with the input image itself; add
## option to provide different image

##' Estimate symmetric blur kernel assuming Gaussian image gradient prior
##'
##' Estimates the blur kernel from a blurred image assuming that the
##' kernel is symmetric and that image gradients follow a Gaussian
##' prior.
##' 
##' @title Estimate symmetric blur kernel
##' @param y Input blurred image, a \code{"rip"} object or something
##'     that can be coerced to one. Color images are handled by first
##'     converting to grayscale using
##'     \code{\link[rip.opencv]{rip.desaturate}}.
##' @param kdim Extent of the kernel. The actual size of the estimated
##'     kernel is \code{2 * kdim + 1}, before being trimmed if
##'     requested.
##' @param resize Factor by which the estimate is to be resized before
##'     being returned.
##' @param g.method Method used to compute variance of Fourier
##'     coefficients. Ignored if \code{g} is explicitly specified. The
##'     \code{"none"} method assumes constant variance, corresponding
##'     to a prior that assumed uncorrelated gradients. The
##'     \code{"autoreg"} method uses a variance model corresponding to
##'     a prior with gradients having a lag-1 auto-regressive
##'     structure. The \code{"nonpar"} method estimates the variance
##'     by smoothing the absolute DFT coefficients; this is of dubious
##'     usefulness but included to allow experimentation.
##' @param g Optional list with components named \code{h} and
##'     \code{v}, giving the variances of DFT coefficients for
##'     horizontal and vertical gradients respectively. Should be of
##'     same size as \code{y}.
##' @param rho.along Correlation coefficient defining auto-regressive
##'     gradient prior model along the direction of the gradient.
##' @param rho.across Correlation coefficient defining auto-regressive
##'     gradient prior model across the direction of the gradient.
##' @param eta.sq Variance of the error terms (after taking gradients)
##' @param corr.grad Whether the error terms should be assumed to be
##'     correlated after taking gradients. If \code{FALSE}, the error
##'     gradients are assumed to be independent; if \code{TRUE}, the
##'     errors in \code{y} are assumed to be independent, and
##'     correlation in gradient errors computed accordingly. In
##'     practice, there is not much difference for reasonable values
##'     of \code{eta.sq}.
##' @param np.span Passed to \code{g.nonpar} when
##'     \code{g.method="nonpar"}.
##' @param trim Logical flag; whether the estimated kernel should be
##'     "trimmed" by setting small elements to zero and removing
##'     all-zero boundary rows.
##' @param zap.digits The \code{digits} argument to
##'     \code{\link{zapsmallp}} to use when \code{trim=TRUE}.
##' @param edgetaper Logical, whether \code{\link{rip.edgetaper}}
##'     should be applied on the gradients before further processing.
##' @return A \code{"rip"} object giving the estimated blur kernel.
symmetric.blur <-function(y, kdim = round(dim(y) / 3), resize = 1,
                          g.method = c("autoreg", "nonpar", "none"),
                          g = NULL,
                          rho.along = 0.3, rho.across = 0.6,
                          eta.sq = 0.005^2,
                          corr.grad = TRUE,
                          np.span = 2/3,
                          trim = TRUE, zap.digits = 1.5,
                          edgetaper = FALSE)
{
    y <- rip.desaturate(y) # in case y is a color image
    ## even size leads to ambiguity in center pixel, so avoid
    if (nrow(y) %% 2 == 0) y <- as.rip(y[-1, ])
    if (ncol(y) %% 2 == 0) y <- as.rip(y[, -1])
    g.method <- match.arg(g.method)
    ## g can be provided as list(h = ..., v = ...)
    if (is.null(g))
        g <- switch(g.method,
                    autoreg = g.autoreg(y, rho = c(rho.along, rho.across), valid = FALSE),
                    nonpar = g.nonpar(y, valid = FALSE, span = np.span),
                    none = list(h = 1, v = 1))
    Wh <- g$h
    Wv <- g$v
    ## filter with same size (don't drop row / column)
    yh <- rip.filter(y, rip.grad$x)
    yv <- rip.filter(y, rip.grad$y)
    if (edgetaper) {
        yh <- rip.edgetaper(yh)
        yv <- rip.edgetaper(yv)
    }
    Yh <- Mod(rip.ndft(yh))
    Yv <- Mod(rip.ndft(yv))
    symk <- function(Y, W = 1, H = 1, kdim)
    {
        M <- sqrt(pmax(Y^2 - eta.sq * H, 0) / W)
        ##   sqrt(pmax((Y^2/W) - eta.sq * H, 0)) # old - wrong ~ eta=0
        k2 <- rip.shift(pmax(rip.ndft(M + 0i, inverse = TRUE), 0))
        cpos <- (dim(k2) + 1) / 2 # guaranteed to be integer, see above
        rlim <- seq(-kdim[1], kdim[1]) + cpos[1]
        clim <- seq(-kdim[2], kdim[2]) + cpos[2]
        as.rip(k2[rlim, clim])
    }
    if (corr.grad)
    {
        Hh <- h.theoretical(dim(Yh))$h
        Hv <- h.theoretical(dim(Yv))$v
    }
    else Hh <- Hv <- 1
    symk.hat <- 0.5 * (symk(Yh, W = Wh, H = Hh, kdim = kdim) +
                       symk(Yv, W = Wv, H = Hv, kdim = kdim))
    if (resize != 1) symk.hat <- rip.resize(symk.hat, fx = resize)
    if (trim) symk.hat <- ktrim0odd(zapsmallp(symk.hat, digits = zap.digits))
    symk.hat[] <- symk.hat / sum(symk.hat)
    symk.hat
}


