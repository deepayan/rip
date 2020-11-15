

posterior.variance <-
    function(y, k, x, lambda, alpha = 2,
             rho = list(along = 0, across = 0),
             yerror = c("normal", "huber", "bisquare", "poisson"),
             huber.k = 1.345, bisquare.c = 4.685, # must be > 1.548 but not checked
             wt.thres = 0.01)
{
    ## For now, don't split, but splitting won't be difficult if needed
    yerror <- match.arg(yerror)
    force.full.rank <- FALSE
    rank.offset <- 1 - force.full.rank
    rank.convtype <- if (force.full.rank) "same" else "valid"
    dlatent <- dim(x) # size of latent image
    Tk <- conv2sparse(k, dlatent[1], dlatent[2])
    Td.h <- conv2sparse(rip.flip(rip.grad$x), dlatent[1], dlatent[2],
                        conv.type = rank.convtype)
    Td.v <- conv2sparse(rip.flip(rip.grad$y), dlatent[1], dlatent[2],
                        conv.type = rank.convtype)
    ## Update with AR model if applicable
    if (rho$across != 0 || rho$along != 0)
    {
        hconv <- with(rho, as.rip(matrix(c(along * across, -along, -across, 1), 2, 2)))
        vconv <- with(rho, as.rip(matrix(c(along * across, -across, -along, 1), 2, 2)))
        hconv[] <- hconv / sqrt(sum(hconv^2)) # normalize variance
        vconv[] <- vconv / sqrt(sum(vconv^2)) # normalize variance
        Td.h <- conv2sparse(rip.flip(hconv), dlatent[1], dlatent[2] - rank.offset,
                            conv.type = "same") %*% Td.h
        Td.v <- conv2sparse(rip.flip(vconv), dlatent[1] - rank.offset, dlatent[2],
                            conv.type = "same") %*% Td.v
    }
    ## We only use the diagonal, and getting the full crossproduct is
    ## expensive.
    dcrossprod <- function(x) colSums(x^2)
    if (alpha == 2 && yerror == "normal")
    {
        A <- dcrossprod(Tk) + 2 * lambda * (dcrossprod(Td.h) + dcrossprod(Td.v))
    }
    else
    {
        yy <- as.vector(y)
        xx <- as.vector(x)
        u.h <- Td.h %*% xx
        u.v <- Td.v %*% xx
        sparsewts.h <- pmax(abs(u.h), wt.thres)^(alpha-2)
        sparsewts.v <- pmax(abs(u.v), wt.thres)^(alpha-2)
        Wh <- Diagonal(nrow(Td.h), sqrt(as.vector(sparsewts.h)))
        Wv <- Diagonal(nrow(Td.v), sqrt(as.vector(sparsewts.v)))
        ## Weights for robust loss
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
                         Diagonal(x = sqrt(pmax(ww, wt.thres)))
                     },
                     huber = {
                         ee <- as.vector(yy - Tk %*% xx) # may contain NAs
                         ## s <- 0.005 OR
                         s <- median(abs(ee), na.rm = TRUE) / 0.6745 # MAD
                         ww <- wt.huber(ee / s)
                         Diagonal(x = sqrt(pmax(ww, wt.thres)))
                     },
                     bisquare = {
                         ee <- as.vector(yy - Tk %*% xx) # may contain NAs
                         ## s <- 0.005 # OR
                         s <- median(abs(ee), na.rm = TRUE) / 0.6745 # MAD
                         ww <- wt.bisquare(ee / s)
                         Diagonal(x = sqrt(pmax(ww, wt.thres)))
                     })
        Tkk <- Wy %*% Tk
        A <- dcrossprod(Tkk) + alpha * lambda * (
            dcrossprod(Wh %*% Td.h) + dcrossprod(Wv %*% Td.v)
        )
    }
    ## NOTE: A is now just the diagonal, not the full
    S <- x
    S[] <- 1 / as.vector(A)
    S
}
