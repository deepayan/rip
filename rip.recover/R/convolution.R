

## conv2valid() gives convolution retaining only pixels which are
## fully defined by pixels in x (i.e., does not need to extrapolate
## values of x outside the boundary). Multiple kernels are allowed to
## operate sequentially. Similar to MATLAB conv2(x, k, "valid")

conv2valid <- function(x, ...)
{
    for (k in list(...)) x <- rip.conv(x, k, "valid")
    return (x)
}

## similarly for full and same

conv2full <- function(x, ...)
{
    for (k in list(...)) x <- rip.conv(x, k, "full")
    return (x)
}

conv2same <- function(x, ...)
{
    for (k in list(...)) x <- rip.conv(x, k, "same")
    return (x)
}

## conv2full with flipped k is the 'transpose' of conv2valid (viewed
## as a linear operation).

## For super-resolution convolution is followed by sub-sampling rows /
## columns. The 'transpose' of this operation is also required for the
## conjugate gradient method

convdown2valid <- function(x, k, factor)
{
    ## str(list(x = x, k = k, factor = factor))
    ## sub-sampling must include the first element (after convolution)
    ## by convention, as otherwise the transpose operation becomes
    ## ambiguous.
    sub <- logical(factor); sub[1] <- TRUE
    as.rip(conv2valid(x, k)[sub, sub, drop = FALSE])
}

convup2full <- function(x, k, factor)
{
    ## We need dim() of the output z. We know that
    ##
    ## (zdim - dim(k) + 1)[sub, sub] == dim(y),
    ##
    ## so there is some ambiguity due to dependence on sub; e.g.
    ## (1:3)[c(T,F,T)] and (1:4)[c(T,F,T,F)] both have length 2. To
    ## resolve this, we assume everywhere that
    zdim <- factor * dim(x) + dim(k) - 1
    z <- matrix(0, zdim[1], zdim[2])
    kp <- rip.pad(k, dim(k) + factor - 1, offset = rep(0, 2))
    for (i in seq_len(factor))
        for (j in seq_len(factor))
        {
            isub <- logical(factor); isub[i] <- TRUE
            jsub <- logical(factor); jsub[j] <- TRUE
            ksub <- kp[isub, jsub, drop = FALSE]
            zup <- conv2full(x, as.rip(ksub))
            ## str(list(x = dim(x), z = dim(z), k = dim(k), kp = dim(kp),
            ##          zsub = dim(z[sub, sub]),
            ##          ksub = dim(ksub),
            ##          zup = dim(zup)))
            z[isub, jsub] <- zup
        }
    z
}


