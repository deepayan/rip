
## Workhorse functions for solving underconstrained linear image
## estimation problems using L2 penalty / Gaussian prior on image
## gradients, using direct solution of linear equations using sparse
## matrix calculations (using the Matrix package)


## Solve linear equation of the form Ax = b. Essentially solve(A, b),
## but forcing Cholesky(perm = TRUE, super = TRUE) as this is
## considerably faster and more memory-efficient in our use-cases.

## Note that we do often solve multiple b with the same A, so re-using
## a one-time decomposition may be worthwhile, but this is only
## helpful for the initial L2 step which is usually fast anyway, and
## not for the IRLS steps, where A changes with current estimate of
## b. In any case, this optimization happens automatically if the
## calling function is careful not to rewrite A, because just calling
## Cholesky(A) automatically populates A@factors for subsequent
## re-use.

## On the other hand, the sparsity structure of A often remains
## unchanged for consecutive calls, especially during IRLS updates,
## suggesting that we may benefit from using update(C, parent =
## newA). However, experiments show that this optimization is only
## marginally beneficial at best, so it is has not been implemented.

solve.full <- function(A, b)
{
    C <- try(Cholesky(A, perm = TRUE, super = TRUE), silent = TRUE)
    if (inherits(C, "try-error"))
        C <- Cholesky(A, perm = TRUE, super = FALSE)
    solve(C, b)
}


## Divide the A matrix into 2x2 blocks and solve diagonals. This is an
## approximation to solve(A, b) that uses less memory, and is an
## alternative (if we generalize to having more blocks) to pre-dividing
## the image into patches. However, this is not really that useful,
## and is only retained for possible experimentation.

solve.by.parts <- function(A, b, parts = 2, overlap = 0.1,
                           randperm = FALSE, verbose = FALSE)
{
    if (parts != 2) warning("'parts != 2' not yet supported; ignored.")
    if (randperm)
    {
        rperm <- sample(nrow(A))
        invperm <- order(rperm)
        A <- A[rperm, rperm]
        b <- b[rperm]
    }
    ## Assume image patch is square (without checking). Otherwise need
    ## to refine calculations.
    N <- nrow(A)
    M <- round(sqrt(N)) # should really have original image dimensions
    H <- round((0.5 + overlap) * M)
    x <- b # placeholder
    solveBatch <- function(part)
    {
        part <- as.vector(part)
        xpart <- solve.full(A[part, part], b[part])
        x[part] <<- xpart
        NULL
    }
    if (verbose) cat("Solving by parts")
    ## block 1,1 = 1:H, 1:H
    part11 <- outer(0:(H-1), seq(1, by = M, length.out = H), "+")
    solveBatch(part11)
    ## block 2,1 = H:M, 1:H (modulo rounding)
    part <- part11 + M - H
    solveBatch(part)
    ## block 1,2 = 1:H, H:M
    part <- part11 + (M - H) * M
    solveBatch(part)
    ## block 2,1 = H:M, H:M
    part <- part11 + (M - H) * (M + 1)
    solveBatch(part)
    if (randperm) x[invperm] else x
}



direct.gaussian <-
    function(A, b, ...,
             latent.dim = NULL,
             verbose = FALSE,
             by.parts = FALSE,
             randperm = FALSE,
             restrict01 = TRUE)
{
    x <- if (by.parts) solve.by.parts(A, b, randperm = randperm)
         else solve.full(A, b)
    if (restrict01)
    {
        x[x < 0] <- 0
        x[x > 1] <- 1
    }
    ## this is the dim of the x to be estimated, product should match ncol(A)
    dim(x) <- latent.dim
    stopifnot(prod(dim(x)) == nrow(A))
    x
}

## Main workhorse:
## - splits input image,
## - sets up linear equations
## - performs IRLS iterations
## - merges results back

direct.irls <-
    function(Tk, y, lambda = 0.001, alpha = 2, Td.h, Td.v,
             super.factor = 1,
             patch = round(100 / super.factor),
             overlap = round(20 / super.factor),
             latent.dim = NULL,
             ...,
             x.start = NULL,
             yerror = c("normal", "huber", "bisquare", "poisson"),
             huber.k = 1.345, bisquare.c = 4.685, # must be > 1.548 but not checked
             wt.thres = 0.01, niter.irls = 5, verbose = FALSE, label = "")
{
    stopifnot(inherits(y, "rip"))
    yerror <- match.arg(yerror)
    if (nchannel(y) != 1)
        stop("Multi-channel images not supported; call each channel separately.")
    ysplit <- splitImage(y, patch = patch, overlap = overlap)
    xsplit <- ysplit # placeholder
    if (alpha == 2 && !is.null(x.start))
        stop("Unexpected input: x.start specified with alpha=2")
    ## First solve unweighted L2 problem, unless x.start is already specified
    if (is.null(x.start))
    {
        ## The A matrix is the same for all splits (unless there are
        ## NA values), so sompute it once and pass on
        A <- crossprod(Tk) + 2 * lambda * (crossprod(Td.h) + crossprod(Td.v))
        for (j in seq_along(ysplit))
        {
            if (verbose)
                cat(sprintf("\r%s[A=%g, L=%g, SR=%d][ Least square estimate: split %d / %d ]           ",
                            label, alpha, lambda, super.factor, j, length(ysplit)))
            yy <- as.vector(ysplit[[j]])
            xsplit[[j]] <-
                if (anyNA(yy))
                {
                    keep <- !is.na(yy)
                    Tkk <- Tk[keep, ]
                    AA <- crossprod(Tkk) + 2 * lambda * (crossprod(Td.h) + crossprod(Td.v))
                    direct.gaussian(AA, b = crossprod(Tkk, yy[keep]),
                                    ...,
                                    latent.dim = latent.dim,
                                    verbose = verbose)
                }
                else
                    direct.gaussian(A, b = crossprod(Tk, yy),
                                    ...,
                                    latent.dim = latent.dim,
                                    verbose = verbose)
        }
        if (verbose) cat("\r                                                                                         \r")
        if (alpha == 2 && yerror == "normal")
            return(unsplitImage(xsplit, enlarge.factor = super.factor))
    }
    else
    {
        ## Careful: the attributes of xsplit are already set
        ## (corresponding to y), but results are possibly larger by
        ## super.factor. If x.start is specified, they must be split
        ## into accordingly larger pieces. (FIXME: check that this is
        ## always the right thing to do)
        xsplit[] <- splitImage(x.start,
                               patch = super.factor * patch,
                               overlap = super.factor * overlap)
    }
    ## Next do IRLS if (a) alpha != 2, or (b) yerror != "normal". 
    ## Process one sub-image (patch) at a time 
    for (j in seq_along(ysplit))
    {
        xx <- as.vector(xsplit[[j]])
        yy <- as.vector(ysplit[[j]])
        for (i in seq_len(niter.irls))
        {
            if (verbose)
                cat(sprintf("\r%s[A=%g, L=%g, SR=%d][ IRLS iteration: % d / %d for split %d / %d ]           ",
                            label, alpha, lambda, super.factor, i, niter.irls, j,
                            length(ysplit)))
            ## Prior weights based on current estimate of x:
            u.h <- Td.h %*% xx
            u.v <- Td.v %*% xx
            sparsewts.h <- pmax(abs(u.h), wt.thres)^(alpha-2)
            sparsewts.v <- pmax(abs(u.v), wt.thres)^(alpha-2)
            Wh <- Diagonal(nrow(Td.h), sqrt(as.vector(sparsewts.h)))
            Wv <- Diagonal(nrow(Td.v), sqrt(as.vector(sparsewts.v)))
            ## y-error weights (yerror = poisson or yerror = huber)
            wt.huber <- function(u)
            {
                pmin(1, huber.k / abs(u))
            }
            wt.bisquare <- function(u)
            {
                (1 - pmin(1, abs(u/bisquare.c))^2)^2
            }
            Wy <- switch(yerror,
                         normal = Diagonal(nrow(Tk)),
                         poisson = { # FIXME: check
                             mu <- as.vector(Tk %*% xx)
                             mu[] <- mu / mean(mu, na.rm = TRUE) # make 1 on average
                             ww <- 1 / mu
                             ## str(range(ww, finite = TRUE))
                             Diagonal(x = sqrt(pmax(ww, wt.thres)))
                         },
                         huber = {
                             ee <- as.vector(yy - Tk %*% xx) # may contain NAs
                             ## s <- 0.005 OR
                             s <- median(abs(ee), na.rm = TRUE) / 0.6745 # MAD
                             ## str(s)
                             ww <- wt.huber(ee / s)
                             Diagonal(x = sqrt(pmax(ww, wt.thres)))
                             ## Diagonal(x = sqrt(ww))
                         },
                         bisquare = {
                             ee <- as.vector(yy - Tk %*% xx) # may contain NAs
                             ## s <- 0.005 # OR
                             s <- median(abs(ee), na.rm = TRUE) / 0.6745 # MAD
                             ## str(s)
                             ww <- wt.bisquare(ee / s)
                             ## str(ww)
                             Diagonal(x = sqrt(pmax(ww, wt.thres)))
                             ## Diagonal(x = sqrt(ww))
                         })
            xx[] <-
                if (anyNA(yy))
                {
                    keep <- !is.na(yy)
                    Tkk <- Wy[keep, keep] %*% Tk[keep, ]
                    A <- crossprod(Tkk) + alpha * lambda * (
                        crossprod(Wh %*% Td.h) + crossprod(Wv %*% Td.v)
                    )
                    direct.gaussian(A, b = crossprod(Tkk, Wy[keep, keep] %*% yy[keep]),
                                    ...,
                                    latent.dim = latent.dim,
                                    verbose = verbose)
                }
                else
                {
                    Tkk <- Wy %*% Tk
                    A <- crossprod(Tkk) + alpha * lambda * (
                        crossprod(Wh %*% Td.h) + crossprod(Wv %*% Td.v)
                    )
                    direct.gaussian(A, b = crossprod(Tkk, Wy %*% yy), 
                                    ...,
                                    latent.dim = latent.dim,
                                    verbose = verbose)
                }
        }
        xsplit[[j]][] <- xx
    }
    if (verbose) cat("\r                                                                                         \r")
    unsplitImage(xsplit, enlarge.factor = super.factor)
}

