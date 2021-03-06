
## Tools to split and merge images to make individual problems smaller
## (and possibly parallelizable). We refer to the individual
## sub-images as image "patches".

## This only anticipates two use cases: the recovered image (to be
## created by merging recovered image patches)
##
## - will be of the same size as the input image (to be split
##   into patches), as in deconvolution or denoising, or
##
## - they will be larger by a constant factor on both dimensions (as
##   in super-resolution or upscaling)


## Computes the sub-image row / column extents (one dimension) given
## total length (n), size of each patch (patch), and amount of overlap
## (overlap). Currently, if the last patch exceeds the full size, then
## it is shifted back, which will make the overlap between the last
## and second-to-last pacthes larger than the other overlaps. FIXME:
## It may be better to instead pre-adjust 'overlap' so that the amount
## of excess in the last image is minimized, keeping the same number
## of resulting patches.

computeSplitRanges <- function(n, patch, overlap)
{
    if (patch > n) {
        patch <- n
        overlap <- 0
    }
    gap <- patch - overlap
    ends <- seq(patch, n + gap - 1, by = gap)
    starts <- ends - patch + 1
    ## last patch may exceed boundary: need to shift back
    k <- length(ends)
    if (ends[k] > n) {
        adjust <- ends[k] - n
        starts[k] <- starts[k] - adjust
        ends[k] <- ends[k] - adjust
    }
    list(starts = starts, ends = ends)
}

## old version

## splitImage <- function(x, patch = 100, overlap = 22)
##     ## x is expected to be a "rip" object
##     ## FIXME: allow for multi-channel objects
## {
##     x <- as.rip(x)
##     if (nchannel(x) != 1) stop("Multi-channel images are not yet supported")
##     ## str(x[], give.attr = FALSE)
##     patch <- rep(patch, length.out = 2)
##     overlap <- rep(overlap, length.out = 2)
##     srow <- computeSplitRanges(nrow(x), patch[1], overlap[1])
##     scol <- computeSplitRanges(ncol(x), patch[2], overlap[2])
##     res <- matrix(list(), nrow = length(srow$starts), ncol = length(scol$starts))
##     for (i in seq_along(srow$starts))
##         for (j in seq_along(scol$starts))
##             res[[i, j]] <-
##                 x[srow$starts[i]:srow$ends[i],
##                   scol$starts[j]:scol$ends[j]]
##     attr(res, "split") <- list(row = srow, col = scol)
##     res
## }


splitImage <- function(x, patch = 100, overlap = 22)
    rip.split(x, patch = patch, overlap = overlap, as.array = FALSE)


enlargeSplitSizes <- function(s, factor = 1)
{
    if (length(factor) != 1 || factor != round(factor))
        stop("enlargement factor must be scalar integer.")
    ## s$row and s$col define split. They need to be adjusted for an
    ## image double the size.
    ## Procedure:
    ##    starts |--> 2 * (starts-1) + 1
    ##      ends |--> 2 * ends
    within(s,
    {
        row$starts <- factor * (row$starts-1) + 1
        col$starts <- factor * (col$starts-1) + 1
        row$ends <- factor * row$ends
        col$ends <- factor * col$ends
    })
}


adjustSplitSizesFull <- function(s, h = c(0, 0))
{
    ## s$row and s$col define split. They need to be adjusted so that
    ## all the starts are decreased by h and all the ends are
    ## increased by h. To make sure we start at 1, finally shift
    ## everything forward by h. Effectively, this means (1) leave
    ## starts alone, and (2) increase ends by 2*h.
    within(s,
    {
        row$ends <- row$ends + 2 * h[1]
        col$ends <- col$ends + 2 * h[2]
    })
}


unsplitImage <- function(xsplit, enlarge.factor = 1, full = FALSE)
{
    ## For deconvolution problems, xsplit[[i]] are slightly larger
    ## than corresponding ysplit[[i]] (due to inclusion of
    ## outside-boundary elements implicitly included via a blur
    ## kernel). By convention, we trim on all sides (assuming the
    ## kernel was centered) to fit the image back in.
    ##
    ## The amount of trimming required can be inferred from the
    ## mismatch in the size of the images in 'xsplit' and the extents
    ## recorded in the "split" attribute (corresponding to the split
    ## of the original image).
    ## 
    ## FIXME: Do we need edge-distance weighted recombination? Seems
    ## OK as it is
    s <- attr(xsplit, "split")
    if (is.null(s)) stop("Information about split is missing.")
    if (!missing(enlarge.factor)) s <- enlargeSplitSizes(s, enlarge.factor)
    ## latent and output (to be retained) dimensions of each
    ## patch. Although we do not check, these should be same for all
    ## patches.
    latent.dim <- dim(xsplit[[1]])
    output.dim <- c(s$row$ends[1], s$col$ends[1])
    ## number of rows / columns to remove on each side
    h <- (latent.dim - output.dim) / 2
    stopifnot(all(h == round(h)))
    ## By convention, we usually want the same size as input, but
    ## sometimes it is useful to do further work with the "full"
    ## latent image.
    if (full) s <- adjustSplitSizesFull(s, h = h)
    ## str(list(s = s, h = h))
    srow <- s$row
    scol <- s$col
    x <- matrix(NA_real_, max(srow$ends), max(scol$ends))
    for (i in seq_along(srow$starts))
        for (j in seq_along(scol$starts))
        {
            xpart <- xsplit[[i, j]]
            if (!full)
            {
                ## need to trim borders to remove h rows/columns on all sides
                xdim <- dim(xpart) 
                xpart <- xpart[seq(h[1] + 1, xdim[1] - h[1]),
                               seq(h[2] + 1, xdim[2] - h[2])]
            }
            ## image(as.rip(xpart), main = sprintf("%d, %d", i, j))
            rrange <- srow$starts[i]:srow$ends[i]
            crange <- scol$starts[j]:scol$ends[j]
            rpart <- 1:nrow(xpart)
            cpart <- 1:ncol(xpart)
            ## crude weighting: remove half of overlap with previous segment
            if (i > 1)
            {
                orow <- srow$ends[i-1] - srow$starts[i] + 1
                if (orow > 1) {
                    skiprow <- round(orow / 2)
                    rrange <- tail(rrange, -skiprow)
                    rpart <- tail(rpart, -skiprow)
                }
            }
            if (j > 1)
            {
                ocol <- scol$ends[j-1] - scol$starts[j] + 1
                if (ocol > 1) {
                    skipcol <- round(ocol / 2)
                    crange <- tail(crange, -skipcol)
                    cpart <- tail(cpart, -skipcol)
                }
            }
            x[rrange, crange] <- as.matrix(xpart[rpart, cpart])
        }
    return(as.rip(x))
}


## splitting / unsplitting is generally useful for other purposes as
## well, e.g., to extract patches, process them, and restore
## them. Assuming that all patches are of the same dimension (possibly
## overlapping), we can either store them as a list (as above) or a
## 3-dimensional array. Both have their uses (arrays can do matrix ops
## easily, lists are easier to resize), and we could think of
## converting between them.

## For now, the (exported) interface below only implements the 3-D
## array representation


rip.split <- function(x, patch = 8, overlap = 4, as.array = TRUE)
    ## x is expected to be a "rip" object
    ## FIXME: allow for multi-channel objects
{
    x <- as.rip(x)
    if (nchannel(x) != 1) stop("Multi-channel images are not yet supported")
    ## str(x[], give.attr = FALSE)
    patch <- rep(patch, length.out = 2)
    overlap <- rep(overlap, length.out = 2)
    srow <- computeSplitRanges(nrow(x), patch[1], overlap[1])
    scol <- computeSplitRanges(ncol(x), patch[2], overlap[2])
    pdim <- c(length(srow$starts), length(scol$starts))
    if (as.array)
    {
        ans <- array(0, dim = c(patch, prod(pdim)))
        for (j in seq_len(pdim[2]))
            for (i in seq_len(pdim[1]))
                ans[, , i + (j - 1L) * pdim[1]] <-
                    x[srow$starts[i]:srow$ends[i],
                      scol$starts[j]:scol$ends[j]]
        class(ans) <- c("rip.split", "array")
    }
    else # store as list with dimensions
    {
        ans <- matrix(list(), nrow = pdim[1], ncol = pdim[2])
        for (j in seq_len(pdim[2]))
            for (i in seq_len(pdim[1]))
                ans[[i, j]] <-
                    x[srow$starts[i]:srow$ends[i],
                      scol$starts[j]:scol$ends[j]]
        class(ans) <- c("rip.split", "list")
    }
    attr(ans, "split") <- list(row = srow, col = scol)
    ans
}


rip.unsplit <- function(xsplit, ...)
{
    if (is.list(xsplit)) unsplitImage(xsplit)
    else rip.unsplit.array(xsplit, ...)
}

rip.unsplit.array <- function(xsplit, merge = c("sum", "avg", "weight"))
{
    merge <- match.arg(merge) # sum is transpose of the rip.split() operator
    s <- attr(xsplit, "split")
    srow <- s$row
    scol <- s$col
    if (is.null(s)) stop("Information about split is missing.")
    pdim <- c(length(srow$starts), length(scol$starts))
    outdim <- c(tail(srow$ends, 1), tail(scol$ends, 1))
    x <- matrix(0, outdim[1], outdim[2])
    if (merge == "weight") stop("weighted merging not implemented yet")
    if (merge == "avg") # nothing to do for "sum"
    {
        ## need to figure out how many values get summed over for each
        ## destination, and average. We do this by keeping track using
        ## an auxiliary array
        xcounts <- matrix(0L, outdim[1], outdim[2])
    }
    for (j in seq_len(pdim[2]))
        for (i in seq_len(pdim[1]))
        {
            x[srow$starts[i]:srow$ends[i],
              scol$starts[j]:scol$ends[j]] <-
                x[srow$starts[i]:srow$ends[i],
                  scol$starts[j]:scol$ends[j]] + xsplit[, , i + (j - 1L) * pdim[1]]
            if (merge == "avg")
            {
                xcounts[srow$starts[i]:srow$ends[i],
                        scol$starts[j]:scol$ends[j]] <-
                    xcounts[srow$starts[i]:srow$ends[i],
                            scol$starts[j]:scol$ends[j]] + 1L
            }
        }
    if (merge == "avg") {
        nzero <- xcounts > 0 # should be all usually
        x[nzero] <- x[nzero] / xcounts[nzero]
    }
    as.rip(x)
}





