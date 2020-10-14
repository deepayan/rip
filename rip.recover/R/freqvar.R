
## variances of Fourier coefficients as function of frequency

## We use the general model where a stationary process z (e.g., image
## gradient) is related to an IID process u through a coordinate-wise
## multiplication of their DFT's, i.e., DFT(z) = DFT(u) * g. This
## function transforms such a stationary process to the corresponding
## IID process.


stationary2iid <- function(z, g = 1, match.sd = FALSE)
{
    if (identical(g, 1)) return(z)
    if (identical(dim(z), dim(g))) # g is variance of DFT coefficients
    {
        u <- rip.dft(rip.dft(z) / sqrt(g), inverse = TRUE)
        if (match.sd) u[] <- (sqrt(mean(z^2)) / sqrt(mean(u^2))) * u;
    }
    else # g is a decorrelating filter
    {
        ## FIXME: this is really only used for the AR filter, where
        ## conv2same() causes ambiguous centering with even
        ## dimensions. Using conv2same() causes problems with the
        ## conjugate gradient method. The solution below is a
        ## temporary workaround that needs to be eventually fixed.
        
        ## u <- conv2same(z, g) ## causes shift
        ## u <- conv2valid(z, g) ## reduces dimension
        ## This seems to work OK, but that's a quirk of how conv2valid
        ## is implemented, so be careful
        u <- as.rip(rbind(0, cbind(0, conv2valid(z, g)))) 
    }
    u
}

## Assuming stationary AR(1) type process with kernel corr(X_t, X_{t+s}) = rho^s

var.AR1 <- function(rho, N, omega = seq_len(N) - 1)
{
    i <- 1i
    theta <- 2 * pi * i * omega / N
    a <- exp(-theta)
    ## when varAR2 uses, rho can be 1
    if (rho == 1)
    {
        if (FALSE)  {
            ## Kaustav's code
            if (a != 1) out <- ((1 - a^N) / (1 - a) * (1 - a^(-N)) / (1 - 1 / a))
            else out <- N^2
        }
        else {
            ## Is this what we actually want?
            out <-
                ifelse(a == 1, N^2,
                       ((1 - a^N) / (1 - a) * (1 - a^(-N)) / (1 - 1 / a)))
        }
    }
    else if (rho == 0) {
        out <- rep(N, N)
    }
    else
    {
        ar <- a * rho
        ard <- a / rho
        rad <- rho / a
        out <- (N * (1 - 1 / (1 - 1 / ar) + rad / (1 - rad)) +
                    1 / (1 - 1 / ar) * (1 - ar^N) / (1 - ar) -
                        (rad^N - 1) / ((1 - rad) * (1 - ard)))
    }
    Re(out / N)
}

## Functions to estimate the DFT variance g (separately for h and v)
## from an AR model, or from an input image by smoothing its DFT (by
## gaussian blurring)

g.autoreg <- function(x, d = dim(x), rho = c(0.3, 0.6), valid = TRUE)
{
    ## x: input image, used only for its dimension
    rho <- as.list(rep(rho, length.out = 2))
    if (is.null(names(rho))) names(rho) <- c("along", "across")
    list(h = as.rip(outer(var.AR1(rho$across, d[1]),
                          var.AR1(rho$along, d[2]-valid))),
         v = as.rip(outer(var.AR1(rho$along, d[1]-valid),
                          var.AR1(rho$across, d[2]))))
}

smoothFourierMag <- function(X, span = 2/3)
{
    L <- round(dim(X) * span)
    S <- outer(dnorm(seq(-2, 2, length = L[1])),
               dnorm(seq(-2, 2, length = L[2])))
    S <- as.rip(S/sum(S))
    ## pdf(file = basename(tempfile(fileext=".pdf")))
    ## image(rip.filter(X, S))
    ## dev.off()
    rip.filter(X, S, borderType = "replicate") # what should be borderType?
}

## FIXME: allow for adjustment by noise variance?
g.nonpar.withk <- function(x, valid = TRUE, span = 2/3,
                           k = NULL, smoothK = TRUE)
{
    ## if k != NULL, use it to adjust 
    dx.h <- if (valid) conv2valid(x, rip.flip(rip.grad$x))
            else rip.filter(x, rip.grad$x)
    dx.v <- if (valid) conv2valid(x, rip.flip(rip.grad$y))
            else rip.filter(x, rip.grad$y)
    Xh <- Mod(rip.ndft(dx.h))
    Xv <- Mod(rip.ndft(dx.v))
    if (!is.null(k))
    {
        Kh <- Mod(rip.dft(k, pad = dim(Xh)))
        Kv <- Mod(rip.dft(k, pad = dim(Xv)))
        if (smoothK) { ## divide by smoothed K instead of raw K
            Kh <- sqrt(smoothFourierMag(Kh^2, span = span * smoothK))
            Kv <- sqrt(smoothFourierMag(Kv^2, span = span * smoothK))
        }
        Xh[] <- Xh / Kh
        Xv[] <- Xv / Kv
    }
    list(h = smoothFourierMag(Xh^2, span = span),
         v = smoothFourierMag(Xv^2, span = span))
}



g.nonpar <- function(x, valid = TRUE, span = 2/3,
                     k = NULL, smoothK = TRUE, niter = 0, ...)
{
    ## initial estimate
    g <- g.nonpar.withk(x, valid = valid, span = span, k = k, smoothK = smoothK)
    for (i in seq_len(niter))
    {
        k <- symmetric.blur(x, resize = 1, g = g, ...)
        k[] <- k / sum(k)
        g <- g.nonpar.withk(x, valid = valid, span = span, k = k, smoothK = smoothK)
    }
    g
}



h.1d <- function(N, omega = seq_len(N) - 1)
{
    Mod(complex(modulus = 1, argument= -2 * pi * omega / N) - 1)
}


h.theoretical <- function(d)
{
    list(h = as.rip(outer(rep(1, d[1]), h.1d(d[2]))),
         v = as.rip(outer(h.1d(d[1]), rep(1, d[2]))))
}



