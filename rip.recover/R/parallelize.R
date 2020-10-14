
## Parallelized interface using the parallel package. This needs more
## work to be useable.

## TODO: Change to give a public API for split() unsplit() and
## otherwise more general

## Must be able to work standalone. Use clusteEvalQ()

processPatch <-
    function(y, k, lambda, alpha, g.h = 1, g.v = 1,
             super.factor = 1,
             latent.dim = NULL,
             ...,
             cg.update,
             poisson.variance = FALSE,
             wt.thres = 0.01, niter.irls = 5, verbose = FALSE, label = "")
{
    require(rip.opencv)
    require(rip.recover)
    iterative.gaussian <- rip.recover:::iterative.gaussian
    stationary2iid <- rip.recover:::stationary2iid
    conv2valid <- rip.recover:::conv2valid
    ## if (verbose)
    ##     cat(sprintf("\r%s[A=%g, L=%g, SR=%d (%s)][ Least square estimate: ]           ",
    ##                 label, alpha, lambda, super.factor,
    ##                 if (identical(g.h, 1) && identical(g.v, 1)) "IID" else "DEP"))
    x <-
        iterative.gaussian(y, k, lambda = 2 * lambda,
                           g.h = g.h, g.v = g.v, W.h = 1, W.v = 1,
                           ...,
                           super.factor = super.factor,
                           cg.update = cg.update,
                           latent.dim = latent.dim,
                           verbose = verbose)
    ## if (verbose) cat("\r                                                                                          \r")
    if (alpha == 2 && !poisson.variance) return(x)
    for (i in seq_len(niter.irls))
    {
        ## if (verbose)
        ##     cat(sprintf("\r%s[A=%g, L=%g, SR=%d (%s)][ IRLS iteration: % d / %d ]           ",
        ##                 label, alpha, lambda, super.factor,
        ##                 if (identical(g.h, 1) && identical(g.v, 1)) "IID" else "DEP",
        ##                 i, niter.irls))
        ## Weights based on current estimate of x:
        u.h <- stationary2iid(conv2valid(x, rip.flip(rip.grad$x)))
        u.v <- stationary2iid(conv2valid(x, rip.flip(rip.grad$y)))
        sparsewts.h <- pmax(abs(u.h), wt.thres)^(alpha-2)
        sparsewts.v <- pmax(abs(u.v), wt.thres)^(alpha-2)
        ## FIXME: add poisson weights as well to include in IRLS
        ## Note: sqrt() not needed for weights in this case
        tmp <- 
            try(iterative.gaussian(y, k, lambda = alpha * lambda,
                               g.h = g.h, g.v = g.v,
                               W.h = sparsewts.h, W.v = sparsewts.v,
                               ...,
                               ## xinit = x, # FIXME: check if useful
                               super.factor = super.factor,
                               cg.update = cg.update,
                               latent.dim = latent.dim,
                               verbose = verbose))
        x[] <- if (inherits(tmp, "try-error")) 1 else tmp
    }
    x
}



iterative.irls.parallel <-
    function(k, y, lambda = 0.001, alpha = 2, g.h = 1, g.v = 1,
             super.factor = 1,
             patch = round(100 / super.factor),
             overlap = round(20 / super.factor),
             latent.dim = NULL,
             cluster = NULL,
             ...,
             cg.update = c("PR", "FR", "BS"),
             x.start = NULL,
             poisson.variance = FALSE,
             wt.thres = 0.01, niter.irls = 5, verbose = FALSE, label = "")
{
    stopifnot(inherits(y, "rip"))
    if (poisson.variance) stop("'poisson.variance = TRUE' not implemented yet.")
    if (nchannel(y) != 1) stop("Multi-channel images not supported yet.")
    cg.update <- match.arg(cg.update)
    ysplit <- splitImage(y, patch = patch, overlap = overlap)
    xsplit <- ysplit # placeholder
    if (!is.null(x.start))
        warning("'x.start' not supported with parallelization... ignoring")
    xsplit[] <-
        parLapplyLB(cl = cluster, ysplit, processPatch,
                    k = k, lambda = lambda, alpha = alpha, g.h = g.h, g.v = g.v,
                    super.factor = super.factor,
                    latent.dim = latent.dim,
                    cg.update = cg.update,
                    ...,
                    poisson.variance = poisson.variance,
                    wt.thres = wt.thres, niter.irls = niter.irls,
                    verbose = verbose, label = label)
    if (verbose) cat("\r                                                                                          \r")
    unsplitImage(xsplit, enlarge.factor = super.factor)
}


