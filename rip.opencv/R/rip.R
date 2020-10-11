
## A rip object represents an image (typically) or other matrix-like
## object in a format similar to OpenCV's native format, to enable
## efficient conversion in C++ code.

## A rip object is always a matrix, with an attribute giving the
## number of rows, columns, and "channels". From the end-user
## perspective, it is more or less enough to know how images are
## represented. Grayscale images have only one channel, and are stored
## in the obvious manner. Multi-channel images have a somewhat
## non-standard representation: An 2x2 RGB image is stored as
##
## BGRBGR
## BGRBGR
##
## That is, for a given image column, channels are stored in
## successive columns, and the convention of ordering is also
## different (BGR instead of RGB).


##' Convert compatible data to a rip object
##'
##' .. content for \details{} ..
##' @title Convert to rip object
##' @param x Matrix containing data of rip object. The default method
##'     handles matrix-like input, where data for multiple channels
##'     should be in successive columns
##' @param channel number of channels in the rip object
##' @param type type of the rip object
##' @param depth depth of the rip object
##' @return converted rip object
##' @author Kaustav Nandy
as.rip <- function(x, ...)
{
    UseMethod("as.rip")
}

as.rip.rip <- function(x, ...) x

as.rip.default <- function(x, channel = 1, ...)
{
    if (!is.matrix(x)) stop("x must be a matrix")
    if (length(channel) != 1) stop("channel must have length one")
    d <- dim(x)
    attr(x, "cvdim") <- c(nrow = d[1], ncol = d[2] / channel,
                          nchannel = unname(channel))
    class(x) <- c("rip", "matrix")
    x
}

as.rip.matrix <- function(x, channel = 1, ...)
{
    as.rip.default(x, channel = channel, ...)
}


## Convert a R 'raster' image to corresponding 'natural' 3-dimensional
## array representation. This maybe be useful to have available at the
## user-level, but for now it's just an internal utility. No alpha
## channel support for now. "raster" objects are stored internally in
## row-major order (see Note in ?as.raster), so we need to apply
## as.matrix() first.

raster2array <- function(x)
{
    x <- as.matrix(x)
    ans <- array(NA_real_, dim = c(dim(x), 3))
    rgb <- col2rgb(x, alpha = FALSE)
    for (i in 1:3) ans[,,i] <- rgb[i, ]
    ans
}


## Convert a 'rip' image to corresponding 'natural' 3-dimensional
## array representation. Note that the order of RGB layers in "rip"
## objects differs from that of "raster" objects, so one should be
## careful when converting. By default, we reverse the ordering (BGR
## <-> RGB) when converting to / from arrays with 3 or 4 channels.
rip2array <- function(x, reverse.rgb = TRUE)
{
    stopifnot(inherits(x, "rip"))
    nc <- nchannel(x)
    ans <- array(NA_real_, dim = c(nrow(x), ncol(x) / nc, nc))
    if (nc == 1) ans[] <- x
    else
    {
        csub <- logical(nc)
        for (i in seq_len(nc))
        {
            isub <- csub
            isub[i] <- TRUE
            ans[,,i] <- x[, isub]
        }
    }
    if (nc == 3 && reverse.rgb) ans[,,3:1] 
    else if (nc == 4 && reverse.rgb) ans[,,c(3:1,4)] else ans
}

as.array.rip <- function(x, reverse.rgb = TRUE, ...) rip2array(x, reverse.rgb)

as.rip.array <- function(x, reverse.rgb = TRUE)
{
    odim <- dim(x)
    if (length(odim) != 3)
        stop("Input array must have exactly 3 dimensions")
    nc <- odim[3] # 3 and 4 (alpha channel) are handled specially
    a <- if (nc == 3 && reverse.rgb) x[, , 3:1, drop = FALSE] # switch RGB to BGR
         else if (nc == 4 && reverse.rgb) x[, , c(3:1, 4), drop = FALSE] # switch RGB to BGR
         else x
    a <- aperm(a, c(1, 3, 2))
    dim(a) <- c(odim[1], prod(odim[c(2, 3)]))
    as.rip.matrix(a, channel = nc)
}

as.rip.raster <- function(x) # ignore alpha channel for now
{
    ## First obtain 'natural' RGB array representation in R
    a <- raster2array(x) # must have exactly 3 channels
    as.rip.array(a)
}

as.rip.nativeRaster <- function(x)
{
    ## First, convert to usual column-major matrix form: usually done
    ## by grDevices:::as.matrix.raster() (see Note in ?as.raster), but
    ## this doesn't work on "nativeRaster" objects.
    x <- matrix(x, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    ## x values are interpreted as 32-bit unsigned integers, with each
    ## byte giving one channel (RGBA).
    R <- bitwAnd(x, 255L); x[] <- bitwShiftR(x, 8)
    G <- bitwAnd(x, 255L); x[] <- bitwShiftR(x, 8)
    B <- bitwAnd(x, 255L); x[] <- bitwShiftR(x, 8)
    A <- bitwAnd(x, 255L)
    a <- array(0, dim = c(dim(x), 4))
    a[,,1] <- B
    a[,,2] <- G
    a[,,3] <- R
    a[,,4] <- A
    as.rip.array(a, reverse.rgb = FALSE)
}


## if rescale=FALSE, values must be interpreted as between 0 and 1
as.raster.rip <- function(x, ..., rescale = TRUE, restrict = TRUE) 
{
    cvdim <- attr(x, "cvdim")
    if (cvdim[3] != 1) {
        dim(x) <- cvdim[c("nrow", "nchannel", "ncol")]
        x <- aperm(x, c(1, 3, 2))
        x <- x[, , 3:1]
    }
    x[] <- if (rescale) 255 * (x - min(x))/(max(x) - min(x))
           else 255 * x
    if (restrict) {
        x[x < 0] <- 0
        x[x > 255] <- 255
    }
    as.raster(unclass(x), max = 255)
}

nchannel <- function(x) attr(x, "cvdim")["nchannel"]



##' Print method for rip objects
##'
##' .. content for \details{} ..
##' @title 
##' @param x An object of class 'rip'
##' @param ... 
##' @return 
##' @author Deepayan Sarkar
print.rip <- function(x, ...)
{
    d <- attr(x, "cvdim")
    cat(sprintf("[%g x %g] image with %g channel(s)\n", d[1], d[2], d[3]))
}


## FIXME: use grid instead

plot.rip <- function(x, rescale = TRUE, ...)
    plot(as.raster(x, rescale = rescale), ...)


image.rip <- function(x, ..., rescale = TRUE, restrict = !rescale,
                      shift = FALSE, xlab = "", ylab = "", asp = 1)
{
    cvdim <- attr(x, "cvdim")
    if (cvdim[3] != 1)
    {
        dim(x) <- cvdim[c("nrow", "nchannel", "ncol")]
        x <- aperm(x, c(1, 3, 2))
        x <- x[,,3:1]
    }
    if (rescale) x[] <- 255 * (x - min(x)) / (max(x) - min(x))
    if (restrict)
    {
        x[x < 0] <- 0
        x[x > 255] <- 255
    }
    if (shift && nrow(x) > 2 && ncol(x) > 2)
    {
        i1 <- floor(nrow(x)/2)
        i2 <- floor(ncol(x)/2)
        x[] <- x[c((i1 + 1):nrow(x), 1:i1),
                 c((i2 + 1):ncol(x), 1:i2)]
    }
    x <- as.raster(unclass(x), max = 255)
    plot(c(0, ncol(x) + 1), c(0, nrow(x) + 1),
         type = "n", xaxs = "i", yaxs = "i",
         xlab = xlab, ylab = ylab, asp = asp, ...)
    rasterImage(x, 1, 1, ncol(x), nrow(x), interpolate = FALSE)
}


