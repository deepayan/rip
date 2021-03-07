
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
                          trim = TRUE, zap.digits = 1.5, show.msc = FALSE,
                          resize.mode = c("image", "kernel"),
                          edgetaper = FALSE)
{
    resize.mode <- match.arg(resize.mode)
    y <- rip.desaturate(y) # in case y is a color image
    if (resize != 1 && resize.mode == "image") y <- rip.resize(y, fx = resize, method = "linear")
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
    Yh2 <- Yh^2 # actually only need Yh2 and Yv2, 
    Yv2 <- Yv^2 # so could simplify code below
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
    ## Can we trim in a more sophisticated way by using AIC / BIC? The
    ## idea is that we try different values of kdim and zap.digits (in
    ## terms of percentage of maximum below which we set to 0), and
    ## choose the best one. This is also ultimately ad hoc, because
    ## the 'cropped / trimmed' k is not the actual constrained MLE.

    ## To do this, we need to calculate log-likelihood for any given
    ## candidate k. We then need an estimate of sigma as well
    ## (ensuring sum(k) == 1), but as we haven't divided by
    ## sum(symk.hat) yet, we can take sigma=1.

    ## In practice, this seems to work well for choosing size (using
    ## BIC), but not so well for zeroing small values

    maxk <- max(symk.hat)

    kAICsize <- function(size) {
        k <- symk.hat
        cpos <- (dim(k) + 1) / 2 # guaranteed to be integer, see above
        rlim <- seq(-size, size) + cpos[1]
        clim <- seq(-size, size) + cpos[2]
        k <- k[rlim, clim]
        npar <- length(k)
        Kh <- Mod(rip.dft(k, pad = dim(Yh))) # NOTE: not rip.ndft() here; otherwise
        Kv <- Mod(rip.dft(k, pad = dim(Yv))) # NOTE: convolution theorem will not hold
        Lh <- eta.sq * Hh + Wh * Kh^2
        Lv <- eta.sq * Hv + Wv * Kv^2
        logL <-
            sum(dchisq(Yh2 / Lh, df = 1, log = TRUE)) +
            sum(dchisq(Yv2 / Lv, df = 1, log = TRUE)) -
            sum(log(Lh)) - sum(log(Lv))
        c(size = size,
          logL = logL,
          AIC = 2 * npar - 2 * logL,
          BIC = npar * log(length(Yh2) + length(Yv2)) - 2 * logL)
    }
    if (show.msc)
    {
        msm <- as.data.frame(t(sapply(seq(max(kdim)-1, 1), kAICsize)))
        print(msm)
    }
    ## plot(AIC ~ size, data = msm)
    ## plot(BIC ~ size, data = msm)
    ## plot(logL ~ size, data = msm)

    ## using model selection criteria for proportion below which to zap (unused)
    ## kAICprop <- function(prop) {
    ##     ## message(prop)
    ##     k <- symk.hat
    ##     k[k < prop * maxk] <- 0
    ##     npar <- sum(k > 0)
    ##     Kh <- Mod(rip.dft(k, pad = dim(Yh))) # NOTE: not rip.ndft() here; otherwise
    ##     Kv <- Mod(rip.dft(k, pad = dim(Yv))) # NOTE: convolution theorem will not hold
    ##     Lh <- eta.sq * Hh + Wh * Kh^2
    ##     Lv <- eta.sq * Hv + Wv * Kv^2
    ##     logL <-
    ##         sum(dchisq(Yh2 / Lh, df = 1, log = TRUE)) +
    ##         sum(dchisq(Yv2 / Lv, df = 1, log = TRUE)) -
    ##         sum(log(Lh)) - sum(log(Lv))
    ##     c(prop = prop,
    ##       logL = logL,
    ##       AIC = 2 * npar - 2 * logL,
    ##       BIC = npar * log(length(Yh2) + length(Yv2)) - 2 * logL)
    ## }
    if (trim)
    {
        ## if (FALSE)
        ## {
        ##     msm <- as.data.frame(t(sapply(seq(0, 0.1, by = 0.005), kAICprop)))
        ##     which.aic <- which.min(msm$AIC)
        ##     aic.prop <- msm$AIC[which.aic]
        ##     ## keep <- with(msm, unique(c(which.max(logL), which.min(AIC), which.min(BIC))))
        ##     ## print(msm[keep,])
        ##     ## print(msm)
        ##     symk.hat[symk.hat < aic.prop * maxk] <- 0
        ##     symk.hat <- ktrim0odd(symk.hat)
        ## }
        ## else 
        symk.hat <- ktrim0odd(zapsmallp(symk.hat, digits = zap.digits))
    }    
    if (resize != 1 && resize.mode == "kernel")
        symk.hat <- rip.resize(symk.hat, fx = resize)
    symk.hat[] <- symk.hat / sum(symk.hat)
    symk.hat
}





