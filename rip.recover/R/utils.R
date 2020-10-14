

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


## Note: interpretation of bandwidth different from density()

make.kernel <-
    function(h = 1,
             kern = c("epanechnikov","rectangular", "triangular",
                      "biweight", "gaussian", "cosine", "optcosine"),
             ..., dim = 2)
{
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
    ## return(list(x = x, y = k))
    if (dim == 1)
    {
        f <- cut(x, seq(-klim, klim, by = 1), include.lowest = TRUE)
        kk <- return(tapply(k, f, mean))
        kk[] <- kk / sum(kk)
    }
    else if (dim == 2)
    {
        ## could just take product measure, but will do slightly better
        f <- cut(x, seq(-klim, klim, by = 1), include.lowest = TRUE)
        ## str(list(f = as.character(f), x = x, k = k))
        g <- expand.grid(f1 = f, f2 = f, KEEP.OUT.ATTRS=FALSE)
        g$k <- as.vector(outer(k, k))
        kk <- xtabs(k ~ f1 + f2, g) / xtabs( ~ f1 + f2, g) # mean
        ## kk[] <- kk / sum(kk)
        as.rip(kk)
    }
    else stop("Only 'dim=1' and 'dim=2' are allowed.")
}

## sort of like zapsmall, but doesn't round

zapsmallp <- function(x, digits = 2, prop = 10^(-digits), threshold = max(abs(x)) * prop)
{
    x[abs(x) < threshold] <- 0
    x
}

