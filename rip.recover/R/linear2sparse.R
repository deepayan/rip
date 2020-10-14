
## Utilities to construct sparse matrices from linear operations
## specified in various specialized forms, with the option of cacheing
## results. The most important use-case (convolution / filter + valid
## / same) is optimized (conv2sparse.fast). FIXME: Improve public API.


sparseRDAPath <- function(prefix, dim, cache.dir = ".")
{
    dir <- file.path(cache.dir, sprintf("scache-%d-%d", dim[1], dim[2]))
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    file.path(dir, sprintf("sparse-%s.rda", prefix))
}



## The AR tranformation (though not any general g) is a
## convolution. We could use that fact and use a variant of
## conv2sparse (that retains size), but this version explicitly
## applies the Fourier-domain AR transform on each elementary vector

ar2sparse <- function(rho, M, N, match.sd = FALSE,
                      ..., # extra arguments passed on to stationary2iid
                      cache = NULL, cache.dir = NULL,
                      giveCsparse = TRUE, zap.digits = 3)
{
    zero <- as.rip(matrix(0, M, N))
    zero <- conv2valid(zero, ...)
    ## The output matrix is NxN
    ## Want to create i,j,x for sparseMatrix(i,j,x).
    ## There will be one such set for each column of the output
    g <- as.rip(outer(var.AR1(rho[1], nrow(zero)),
                      var.AR1(rho[2], ncol(zero))))
    linear2sparse(stationary2iid, M, N,
                  g = g, match.sd = match.sd, ...,
                  cache = cache, cache.dir = cache.dir, 
                  giveCsparse = giveCsparse, zap.digits = zap.digits)
}


## special case of linear2sparse, where FUN is convolution. The
## difference is that the threshold for converting small values to 0
## is global (not per-column), determined as a fraction of max(k).

## To save this in a cache, we need an unique identifier, which we
## cannot really expect the end-user to provide. Unless we can think
## of something better, we take the somewhat convoluted approach of
## saving the kernel in a temporary file and use its md5sum as the
## unique identifier.

kernel2md5 <- function(k)
{
    ftmp <- tempfile()
    save(k, file = ftmp)
    on.exit(unlink(ftmp))
    tools::md5sum(ftmp)
}

conv2sparse <- function(k, M, N, convolution = TRUE, fast = TRUE,
                        cache = NULL, cache.dir = NULL,
                        zap.digits = 4,
                        conv.type = c("valid", "same", "full"))
{
    if (fast)
        return(conv2sparse.fast(k = k, M = M, N = N,
                                convolution = convolution,
                                zap.digits = zap.digits,
                                conv.type = conv.type))
    conv.type <- match.arg(conv.type)
    convFun <- switch(conv.type,
                      valid = conv2valid,
                      same = conv2same,
                      full = conv2full)
    if (!is.null(cache.dir))
    {
        if (is.null(cache)) cache <- kernel2md5(k)
        f <- sparseRDAPath(cache, c(M, N), cache.dir = cache.dir)
        if (file.exists(f)) {
            load(f)
            return(T)
        }
        else cat("             \t[ cache file: ", f, " ]")
    }
    if (!convolution) # flip row / column order (otehrwise filtering)
        k[] <- k[rev(seq_len(nrow(k))), rev(seq_len(ncol(k)))]
    MN <- M * N
    zero <- as.rip(matrix(0, M, N))
    mn <- prod(dim(convFun(zero, k)))
    ## The output matrix is mn x MN, where mn is usually not MN,
    ## because only "valid" entries are retained in the output.  Want
    ## to create i,j,x for sparseMatrix(i,j,x).  There will be one
    ## such set for each of the MN columns of the output
    zap.threshold <- max(abs(k)) / 10^zap.digits
    ijColumn <- function(I, J)
    {
        e <- zero
        e[I, J] <- 1
        zapsmallp(convFun(e, k), threshold = zap.threshold)
    }
    ilist <- jlist <- xlist <- matrix(list(), M, N)
    for (I in 1:M)
    {
        cat(sprintf("\r %5d / %d", I, M))
        for (J in 1:N)
        {
            IJC <- ijColumn(I, J) # dense (dim may be less than M * N)
            wnz <- as.vector(which(IJC != 0))
            wcol <- (J-1L)*M + (I-1L) + 1L # column index in sparse matrix
            ilist[[I, J]] <- wnz
            jlist[[I, J]] <- rep(wcol, length(wnz))
            xlist[[I, J]] <- IJC[wnz]
        }
    }
    cat("\n")
    ## print(list(d = dim(IJC), M = M, N = N))
    T <- sparseMatrix(unlist(ilist), unlist(jlist), x = unlist(xlist),
                      dims = c(mn, MN), giveCsparse = TRUE)
    if (!is.null(cache.dir)) save(T, file = f)
    T
}


## Specialized version of conv2sparse, using kronecker products to
## build up block Toeplitz matrix, with each block a Toeplitz matrix
## corresponding to one column


## unsymmetric centered toeplitz. If x has odd length, middle element
## becomes diagonal. Otherwise, round to left/right? to match with conv2sparse()
ctoeplitz <- function(x, size)
    ## x contains part corresponding to convolution
    ## size is dimension of resulting Toeplitz matrix
{
    n <- length(x)
    shift <- floor((n-1)/2)
    ## create matrix of dimension size+shift, then drop rows/columns
    trow <- toeplitz(c(x, rep(0, size + shift - n)))
    trow[lower.tri(trow)] <- 0
    if (shift == 0)
        Matrix(trow, sparse = TRUE)
    else
        Matrix(trow[-(size+seq_len(shift)), -seq_len(shift)], sparse = TRUE)
}

## Ignore valid for now, just try to reproduce "same" first. The rest
## should be easy, hopefully.

conv2sparseFast.same <-
    function(k, M, N, zap.digits = 4)
{
    zap.threshold <- max(abs(k)) / 10^zap.digits
    k[] <- zapsmallp(k, threshold = zap.threshold)
    ## output is MN x MN matrix
    ## columns of k must be centered as in ctoeplitz()
    n <- ncol(k)
    shift <- floor((n-1)/2)
    eye <- Diagonal(N)
    ## each column of k becomes one block
    DIAG <- 1 + shift
    ## first column for now
    kcol <- c(k[,DIAG,drop=TRUE], rep(0, M-nrow(k)))
    ans <- kronecker(eye, ctoeplitz(k[,DIAG,drop=TRUE], M))
    ## other columns: first above diagonal
    for (i in seq_len(n))
    {
        if (i != DIAG)
        {
            eye[] <- bandSparse(N, N, i-DIAG, list(rep(1, N)))
            ans <- ans + kronecker(eye, ctoeplitz(k[,i,drop=TRUE], M))
        }
    }
    ans
}

conv2sparse.fast <-
    function(k, M, N, convolution = TRUE, zap.digits = 4,
             conv.type = c("valid", "same", "full"))
{
    conv.type <- match.arg(conv.type)
    if (conv.type == "full") stop("conv.type='full' not supported yet.")
    ## FIXME: check but should be just t(conv2sparse.fast(rip.flip(k),
    ## ?, ?, conv.type = "valid"))
    if (convolution) k <- rip.flip(k)
    ans <- conv2sparseFast.same(k, M = M, N = N, zap.digits = zap.digits)
    if (conv.type == "same") return(ans)
    if (conv.type == "valid") # drop incomplete rows
    {
        dummy <- k
        dummy[] <- 1 ## all non-zero, to detect and drop incomplete rows
        dummy.ans <- conv2sparseFast.same(dummy, M = M, N = N)
        rnnz <- rowSums(dummy.ans != 0) # number of non-zero entries
        rmax <- max(rnnz)
        return(ans[rnnz == rmax, ])
    }
}


##' Use a function that operates linearly on matrices to generate a
##' sparse matrix that operates on corresponding vectorized forms.
##'
##' Suppose the function 'FUN' operates on M-by-N real-valued matrices
##' to produce m-by-n real-valued matrices (a typical example would be
##' a 2-D convolution). _If_ this function is linear when considered
##' as a function from the space of M-by-N real matrices R^(MxN) to
##' the space of m-by-n real matrices R^(mxn), then it can in
##' principle also be represented as multiplication by a matrix T on
##' the correponding vectorized forms (i.e., from R^(M*N) to
##' R^(m*n)). In other words, 'FUN(x)' is equivalent to 'T %*%
##' as.vector(x)'.
##'
##' In general, the matrix 'T' may be difficult to describe. This
##' function assumes (without checking) that 'FUN' is indeed linear,
##' and constructs 'T' by applying 'FUN' on each canonical vector in
##' R^(M*N), i.e., all M-by-N matrices with one element set to 1, and
##' all other elements set to 0.
##'
##' As this conversion is not useful unless the resulting 'T' is
##' sparse, the result is represented as a sparse matrix (using the
##' 'Matrix' package).
##' 
##' @title Convert linear matrix-to-matrix function to sparse matrix
##'     form
##' @param FUN A function taking a real-valued matrix as input and
##'     another real-valued matrix as output.
##' @param M Row dimension of the input matrix
##' @param N Column dimension of the input matrix
##' @param ... Additional arguments passed on to 'FUN'
##' @param cache 'NULL', or a string that will be used to generate the
##'     name of a file where the result will be cached. If this file
##'     already exists, then the cache will be used and the actual
##'     computation will be skipped. Note that it is not checked
##'     whether 'FUN' is the same as the function used to generate the
##'     cache. However, the cache file location does encode the value
##'     of 'M' and 'N'.
##' @param cache.dir A directory giving the root of the
##'     cache. Defaults to working directory.
##' @param giveCsparse Passed on to the
##'     \code{\link[Matrix]{sparseMatrix}} constructor.
##' @param zap.digits To control the level of sparsity, results of
##'     'FUN' are modified to set "small" values to zero. This
##'     argument controls which values are considered small enough;
##'     specifically, the threshold is the maximum (per call to 'FUN')
##'     divided by '10^zap.digits'.
##' @return A sparse matrix produced by \code{\link[Matrix]{sparseMatrix}}.
linear2sparse <- function(FUN, M, N, ..., cache, cache.dir = NULL, 
                          giveCsparse = TRUE, zap.digits = 3)
    ## FIXME: add option to suppress progress indicator (e.g., show.progress=FALSE)
{
    if (!is.null(cache.dir))
    {
        f <- sparseRDAPath(cache, c(M, N), cache.dir = cache.dir)
        if (file.exists(f)) {
            load(f)
            return(T)
        }
        else cat("             \t[ cache file: ", f, " ]")
    }
    MN <- M * N
    zero <- as.rip(matrix(0, M, N))
    ## The output matrix is mxn
    ## Want to create i,j,x for sparseMatrix(i,j,x).
    ## There will be one such set for each column of the output
    ilist <- jlist <- xlist <- vector(mode = "list", length = MN)
    process_ij <- function(I, J)
    {
        e <- zero
        e[I, J] <- 1
        IJC <- as.vector(zapsmallp(FUN(e, ...), digits = zap.digits)) # dense
        wnz <- which(IJC != 0)
        wcol <- (J-1L)*M + (I-1L) + 1L # column index in sparse matrix
        ilist[[wcol]] <<- wnz
        jlist[[wcol]] <<- rep(wcol, length(wnz))
        xlist[[wcol]] <<- IJC[wnz]
    }
    mn <- prod(dim(FUN(zero, ...)))
    for (I in 1:M)
    {
        cat(sprintf("\r %5d / %d", I, M))
        for (J in 1:N) process_ij(I, J)
    }
    cat("\n")
    ## pop_sparse_comps(process_ij, M, N)
    ## str(list(ilist[1:5], jlist[1:5], xlist[1:5]))
    T <- sparseMatrix(unlist(ilist), unlist(jlist),
                      x = unlist(xlist), dims = c(mn, MN),
                      giveCsparse = giveCsparse)
    if (!is.null(cache.dir)) save(T, file = f)
    T
}
