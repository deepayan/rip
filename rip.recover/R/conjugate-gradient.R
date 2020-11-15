

## Workhorse functions for solving underconstrained linear image
## estimation problems using L2 penalty / Gaussian prior on image
## gradients, using iterative methods (specifically
## conjugate-gradient) using the optim() function.

## TODO: add robust loss (see sparse-matrix.R)

## To solve Ax=b iteratively, the conjugate gradient algorithm
## minimizes 1/2 x'Ax - b'x. Instead of a custom implementation of the
## algorithm, we want to do it using optim(method = "CG"). For this,
## we need functions
## 
## fn(x) : evaluate objective, which is = 0.5 * x'Ax - b'x
## gr(x) : evaluate gradient, which is  = Ax - b
##
## (Although technically C-G doesn't need fn(x), optim() does require it.) 
## This can be done without actually computing A.

## NOTE: What is Tk' etc when we look at 'valid' convolution? Turns
## out that this is equivalent to 'full' convolution. See
## ../tests/conv2matrix.R

Ax <- function(x, k, lambda, g.h = 1, g.v = 1, W.h = 1, W.v = 1, W.y = 1)
    ## Note: the hyper-Laplacian parameter is incorporated in lambda
{
    ## save(x, k, file = "/tmp/bar.rda")
    Ax <- conv2full(W.y * conv2valid(x, k), rip.flip(k))
    ## horizontal gradient
    z.h <- conv2valid(x, rip.flip(rip.grad$x))
    u.h <- stationary2iid(z.h, g.h)
    Ax.h <- conv2full(stationary2iid(W.h * u.h, g.h), rip.grad$x)
    ## vertical gradient
    z.v <- conv2valid(x, rip.flip(rip.grad$y))
    u.v <- stationary2iid(z.v, g.v)
    Ax.v <- conv2full(stationary2iid(W.v * u.v, g.v), rip.grad$y)
    Ax + lambda * (Ax.h + Ax.v)
}


superAx <- function(x, k, lambda, g.h = 1, g.v = 1, W.h = 1, W.v = 1, W.y = 1, factor)
    ## Note: the hyper-Laplacian parameter is incorporated in lambda
{
    Ax <- convup2full(W.y * convdown2valid(x, k, factor), rip.flip(k), factor)
    ## horizontal gradient
    z.h <- conv2valid(x, rip.flip(rip.grad$x))
    u.h <- stationary2iid(z.h, g.h)
    Ax.h <- conv2full(stationary2iid(W.h * u.h, g.h), rip.grad$x)
    ## vertical gradient
    z.v <- conv2valid(x, rip.flip(rip.grad$y))
    u.v <- stationary2iid(z.v, g.v)
    Ax.v <- conv2full(stationary2iid(W.v * u.v, g.v), rip.grad$y)
    Ax + lambda * (Ax.h + Ax.v)
}

setup_optim <- function(y, k, lambda, g.h, g.v, W.h, W.v, W.y,
                        latent.dim, xinit = NULL, super.factor = 1)
{
    x <- as.rip(matrix(0.5, latent.dim[1], latent.dim[2]))
    if (!is.null(xinit)) x[] <- xinit
    if (super.factor == 1)
    {
        ## str(list(y = y, k = k))
        b <- conv2full(W.y * y, rip.flip(k))
        Ax <- Ax(x, k, lambda = lambda,
                 g.h = g.h, g.v = g.v,
                 W.h = W.h, W.v = W.v, W.y = W.y)
        ## Re-Evaluate Ax only if x has been updated
        ## FIXME: check speed savings, if any (otherwise simpler to always evaluate)
        updateAx <- function(par)
        {
            ## calls <<- calls + 1
            if (!all(par == x))
            {
                x[] <<- par
                Ax[] <<- Ax(x, k, lambda = lambda,
                            g.h = g.h, g.v = g.v,
                            W.h = W.h, W.v = W.v, W.y = W.y)
            }
        }
    }
    else
    {
        b <- convup2full(W.y * y, rip.flip(k), super.factor)
        Ax <- superAx(x, k, lambda = lambda,
                      g.h = g.h, g.v = g.v,
                      W.h = W.h, W.v = W.v, W.y = W.y,
                      factor = super.factor)
        ## Re-Evaluate Ax only if x has been updated
        ## FIXME: check speed savings, if any (otherwise simpler to always evaluate)
        updateAx <- function(par)
        {
            ## calls <<- calls + 1
            if (!all(par == x))
            {
                x[] <<- par
                Ax[] <<- superAx(x, k, lambda = lambda,
                                 g.h = g.h, g.v = g.v,
                                 W.h = W.h, W.v = W.v, W.y = W.y,
                                 factor = super.factor)
            }
        }
    }
    ## calls <- 1
    fn <- function(par)
    {
        updateAx(par)
        ## str(list(fn.x = x, fn.b = b))
        0.5 * sum(x * Ax) - sum(x * b)
    }
    gr <- function(par)
    {
        updateAx(par)
        ## str(list(gr.x = x, gr.b = b))
        Ax - b
    }
    list(par = as.vector(x), fn = fn, gr = gr)
}


iterative.gaussian <-
    function(y, k, lambda, g.h, g.v, W.h, W.v, W.y = 1, cg.update, 
             latent.dim = NULL, super.factor = 1, xinit = NULL,
             ...,
             verbose = FALSE,
             restrict01 = TRUE)
{
    s <- setup_optim(as.rip(y), k, lambda = lambda, g.h = g.h, g.v = g.v,
                     W.h = W.h, W.v = W.v, W.y = W.y, super.factor = super.factor,
                     latent.dim = latent.dim, xinit = xinit)
    cg.type <- switch(cg.update, FR = 1, PR = 2, BS = 3)
    ocontrol <- list(type = cg.type, ...)
    o <- with(s, optim(par = par, fn = fn, gr = gr,
                       method = "CG", control = ocontrol))
    ## if (verbose) message("Convergence code: ", o$convergence)
    x <- o$par
    if (restrict01)
    {
        x[x < 0] <- 0
        x[x > 1] <- 1
    }
    dim(x) <- latent.dim
    x
}


## Main workhorse:
## - splits input image,
## - computes weights
## - performs IRLS iterations
## - merges results back

iterative.irls <-
    function(k, y, lambda = 0.001, alpha = 2, g.h = 1, g.v = 1,
             super.factor = 1,
             patch = round(100 / super.factor),
             overlap = round(20 / super.factor),
             latent.dim = NULL,
             full.latent = FALSE,
             ...,
             cg.update = c("PR", "FR", "BS"),
             x.start = NULL,
             yerror = c("normal", "huber", "bisquare", "poisson"),
             huber.k = 1.345, bisquare.c = 4.685, # must be > 1.548 but not checked
             wt.thres = 0.01, niter.irls = 5, verbose = FALSE, label = "")
{
    stopifnot(inherits(y, "rip"))
    yerror <- match.arg(yerror)
    if (nchannel(y) != 1)
        stop("Multi-channel images not supported; call each channel separately.")
    cg.update <- match.arg(cg.update)
    ysplit <- splitImage(y, patch = patch, overlap = overlap)
    xsplit <- ysplit # placeholder
    if (alpha == 2 && !is.null(x.start))
        stop("Unexpected input: x.start specified with alpha=2")
    ## First solve unweighted L2 problem, unless x.start is already specified
    if (is.null(x.start))
    {
        for (j in seq_along(ysplit))
        {
            if (verbose)
                cat(sprintf("\r%s[A=%g, L=%g, SR=%d (%s)][ Least square estimate: split %d / %d ]           ",
                            label, alpha, lambda, super.factor,
                            if (identical(g.h, 1) && identical(g.v, 1)) "IID" else "DEP",
                            j, length(ysplit)))
                xsplit[[j]] <-
                    iterative.gaussian(ysplit[[j]], k, lambda = 2 * lambda,
                                       g.h = g.h, g.v = g.v, W.h = 1, W.v = 1,
                                       ...,
                                       super.factor = super.factor,
                                       cg.update = cg.update,
                                       latent.dim = latent.dim,
                                       verbose = verbose)
        }
        if (verbose) cat("\r                                                                                          \r")
        if (alpha == 2 && yerror == "normal")
            return(unsplitImage(xsplit, enlarge.factor = super.factor, full = full.latent))
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
        xx <- xsplit[[j]]
        yy <- ysplit[[j]]
        ## x <- as.rip(xsplit[[j]])
        for (i in seq_len(niter.irls))
        {
            if (verbose)
                cat(sprintf("\r%s[A=%g, L=%g, SR=%d (%s)][ IRLS iteration: % d / %d for split %d / %d ]           ",
                            label, alpha, lambda, super.factor,
                            if (identical(g.h, 1) && identical(g.v, 1)) "IID" else "DEP",
                            i, niter.irls, j, length(ysplit)))
            ## Prior weights based on current estimate of x:
            u.h <- stationary2iid(conv2valid(xx, rip.flip(rip.grad$x)))
            u.v <- stationary2iid(conv2valid(xx, rip.flip(rip.grad$y)))
            sparsewts.h <- pmax(abs(u.h), wt.thres)^(alpha-2)
            sparsewts.v <- pmax(abs(u.v), wt.thres)^(alpha-2)
            ## y-error weights (yerror = poisson or yerror = huber)
            wt.huber <- function(u)
            {
                pmin(1, huber.k / abs(u))
            }
            wt.bisquare <- function(u)
            {
                (1 - pmin(1, abs(u/bisquare.c))^2)^2
            }
            Wy <- switch(yerror,        # NOTE: no need to take sqrt() here
                         normal = 1,
                         poisson = { # FIXME: check
                             mu <- conv2valid(xx, k)
                             mu[] <- mu / mean(mu, na.rm = TRUE) # make 1 on average
                             ww <- 1 / mu
                             pmax(ww, wt.thres)
                         },
                         huber = {
                             ee <- yy - conv2valid(xx, k)
                             ## s <- 0.005 OR
                             s <- median(abs(ee), na.rm = TRUE) / 0.6745 # MAD
                             ## str(s)
                             ww <- wt.huber(ee / s)
                             pmax(ww, wt.thres)
                         },
                         bisquare = {
                             ee <- yy - conv2valid(xx, k)
                             ## s <- 0.005 # OR
                             s <- median(abs(ee), na.rm = TRUE) / 0.6745 # MAD
                             ## str(s)
                             ww <- wt.bisquare(ee / s)
                             ## str(ww)
                             pmax(ww, wt.thres)
                         })
            xx[] <-
                iterative.gaussian(ysplit[[j]], k, lambda = alpha * lambda,
                                   g.h = g.h, g.v = g.v,
                                   W.h = sparsewts.h, W.v = sparsewts.v, W.y = Wy,
                                   ...,
                                   ## xinit = x, # FIXME: check if useful
                                   super.factor = super.factor,
                                   cg.update = cg.update,
                                   latent.dim = latent.dim,
                                   verbose = verbose)
        }
        xsplit[[j]] <- xx
    }
    if (verbose) cat("\r                                                                                          \r")
    unsplitImage(xsplit, enlarge.factor = super.factor, full = full.latent)
}


