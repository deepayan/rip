### File import and export

##' Import image data from an external file or write image data to an
##' external file.
##'
##' User-friendly interface to the OpenCV \code{cv::imread} and
##' \code{cv::imwrite} functions via the \code{rip.cv$IO} module. The
##' file format is automatically determined by OpenCV.
##' 
##' @title Import from or export to an image file
##' @param file Path of the image file.
##' @param x An object of class \code{"rip"} containing data that can
##'     be interpreted as an image by OpenCV (usually with 1, 3, or 4
##'     channels with values between 0 and 255).
##' @param type Character string that determines whether the imported
##'     data will be stored as grayscale (default) or color (including
##'     possibly an alpha channel). If \code{"original"}, the choice
##'     is taken from the image.
##' @return image as a rip object.
##' @author Kaustav Nandy
rip.import <- function(file, type = c("grayscale", "color", "original"))
{
    type <- match.arg(type)
    itype <- switch(type, original = -1L, grayscale = 0L, color = 1L)
    file <- path.expand(file)
    if (!file.exists(file))
        stop("Specified file does not exist or has inadequate permissions.")
    rip.cv$IO$imread(file, itype)
}

rip.export <- function(x, file = "")
{
    rip.cv$IO$imwrite(x, file)
}


## Filtering, padding, resizing, etc.

border2enum <-
    function(type = c("constant", "replicate", "reflect", "wrap", "reflect_101"))
{
    ## Suppose the image is abcdef. Then
    ## "constant": 000000|abcdefgh|000000 
    ## "replicate": aaaaaa|abcdefgh|hhhhhhh
    ## "reflect": fedcba|abcdefgh|hgfedcb
    ## "wrap": cdefgh|abcdefgh|abcdefg
    ## "reflect_101": gfedcb|abcdefgh|gfedcba
    type <- match.arg(type)
    switch(type,
           constant = rip.cv$enums$BorderTypes["BORDER_CONSTANT"],
           replicate = rip.cv$enums$BorderTypes["BORDER_REPLICATE"],
           reflect = rip.cv$enums$BorderTypes["BORDER_REFLECT"],
           wrap = rip.cv$enums$BorderTypes["BORDER_WRAP"],
           reflect_101 = rip.cv$enums$BorderTypes["BORDER_REFLECT_101"])
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Forms a border around the rip image object
##' @param x rip object 
##' @param top border size at top direction
##' @param bottom border size at bottom direction
##' @param left border size at left direction
##' @param right border size at right direction
##' @param borderType border type
##' @return bordered image
##' @author kaustav nandy
rip.pad <-
    function(x, pad = optDFTlength(dim(x)),
             offset = round((pad-dim(x))/2),
             value = 0,
             top = offset[1],
             left = offset[2],
             bottom = pad[1] - nrow(x) - top,
             right = pad[2] - ncol(x) - left,
             borderType = c("constant", "replicate", "reflect", "wrap", "reflect_101"))
{
    x <- as.rip(x)
    if (FALSE) # R version
    {
        y <- matrix(value, pad[1], pad[2])
        y[offset[1] + seq_len(nrow(x)), offset[2] + seq_len(ncol(x))] <- x
        as.rip(y)
    }
    bcode <- border2enum(match.arg(borderType))
    rip.cv$core$copyMakeBorder(x, top, bottom, left, right, bcode, value)
}


rip.resize <-
    function(x, d = c(0, 0), fx = 1, fy = fx,
             method = c("nearest", "linear", "area", "cubic", "lanczos4"))
{
    x <- as.rip(x)
    ## In geomtrans$resize(), either d or both fx and fy must be 0
    if (!missing(d)) {
        if (all(d == dim(x))) return(x)
        d <- rev(d) # (row,column) in R -> (x,y) in opencv
        fy <- fx <- 0 # ignored without warning if specified
        mscale <- min(d / dim(x))
    }
    else mscale <- min(fx, fy)
    if (missing(method))
        method <- if (mscale < 1) "area" else "linear"
    else
        method <- match.arg(method)
    iflags <- rip.cv$enums$InterpolationFlags
    icode <- switch(method,
                    nearest = iflags["INTER_NEAREST"],
                    linear = iflags["INTER_LINEAR"],
                    area = iflags["INTER_AREA"],
                    cubic = iflags["INTER_CUBIC"],
                    lanczos4 = iflags["INTER_LANCZOS4"])
    rip.cv$imgproc$resize(x, d, fx, fy, icode)
}


rip.blur <-
    function(x, method = c("mean", "median", "gaussian"),
             ksize = c(1L, 1L), sigma.x = 1, sigma.y = 1,
             anchor = c(-1L, -1L),
             borderType = c("reflect_101", "constant", "replicate", "reflect", "wrap"))
{
    x <- as.rip(x)
    ## Note: medianBlur() may not work for large ksize with double data
    ## TODO: add checks on valid ksize (odd etc)
    method <- match.arg(method)
    ksize <- rep(ksize, length.out = 2)
    anchor <- rep(anchor, length.out = 2)
    bcode <- border2enum(match.arg(borderType))
    switch(method,
           mean =
               {
                   rip.cv$filter$blur(x, ksize, anchor, bcode)
               },
           median =
               {
                   if (!(ksize[1] %in% c(1L, 3L, 5L)))
                       stop("medianBlur() is supported for ksize=3 and 5 only.")
                   rip.cv$filter$medianBlur(x, ksize[1])
               },
           gaussian =
               {
                   if (!all(ksize == 0L) && any(ksize %%2 == 0L))
                       stop("ksize must be 0 or odd in gaussianBlur().")
                   rip.cv$filter$GaussianBlur(x, ksize, sigma.x, sigma.y, bcode)
               })
}

## Filtering (flip=TRUE for convolution)

## borderType WRAP should match standard DFT, but that's not allowed
## by cv::filter2D(). Using BORDER_REPLICATE as default for
## now. FIXME: Not sure how OpenCV uses DFT for 'large' k (more than
## 11x11) - maybe should just use our own that always uses DFT. Need
## to figure out how anchor affects this.

rip.filter <- 
    function(x, k, flip = FALSE, anchor = c(-1L, -1L), delta = 0,
             borderType = c("replicate", "constant", "reflect", "wrap", "reflect_101"))
{
    x <- as.rip(x)
    k <- as.rip(k)
    if (nchannel(k) != 1) stop("filter 'k' must have a single channel")
    if (flip) k[] <- rip.flip(k)
    bcode <- border2enum(match.arg(borderType))
    rip.cv$filter$filter2D(x, k, anchor, delta, bcode)
}


## conv2() in MATLAB has the option of doing "full", "same", or
## "valid" convolution, and this distinction is important for
## us. cv::filter2D only does "same", so we have our own wrapper
## around that, rip.conv(). Use flip = FALSE for correlation filter.



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Matlab style 2-D convolution
##' @param x 'rip' object to be convolved
##' @param k 'rip' object to be convolved with
##' @param type shape parameter which determines the part of the
##' 2-dimensional convolution to be returned. If \code{"full"}, the
##' full convolution matrix is returned, for \code{"same"}, 'rip'
##' object of same size as 'x' is returned, and in case of
##' \code{"valid"}, only the part of the convolution that are without
##' zero padded edges are returned.
##' @param flip 
##' @return Convolution of 'x' and 'k'
##' @author kaustav nandy
rip.conv <- function(x, k, type = c("full", "valid", "same"), flip = TRUE)
{
    ## NOTE: Although cv::filter2D docs don't mention this explicitly,
    ## borderType = "BORDER_CONSTANT" uses value = 0, which is what we
    ## want.
    type <- match.arg(type)
    x <- as.rip(x)
    k <- as.rip(k)
    d <- dim(k)
    ## if (flip) k[] <- rip.flip(k)
    ## flip <- FALSE
    pad.extra <- d - 1 # extra rows / columns for type = "full"
    ## Filtering / convolution happens with kernel centered (anchor =
    ## -1), so the breakup is obvious when the number of rows /
    ## columns in k are odd. When even, "centering" is ambiguous, so
    ## we need to fix a convention; we choose to be consistent with MATLAB.
    pad.bottomright <- floor(pad.extra / 2)
    pad.topleft <- pad.extra - pad.bottomright # == floor((pad.extra + 1) / 2)
    ## Two options: do completely case-by-case, or always do full and
    ## then drop rows / columns if necessary. First option is easier
    ## to understand, so go with that.
    anchor <- ceiling(d / 2 - 1) # == d - floor(d / 2) - 1
    anchor <- rev(anchor) # Need in (x,y) order, not (row,column)
    switch(type,
           same =
           {
               rip.filter(x, k, flip = flip,
                          anchor = anchor,
                          borderType = "constant")
           },
           full =
           {
               rip.filter(rip.pad(x, value = 0,
                                  top = pad.topleft[1],
                                  bottom = pad.bottomright[1],
                                  left = pad.topleft[2],
                                  right = pad.bottomright[2],
                                  borderType = "constant"),
                          k, flip = flip,
                          anchor = anchor,
                          borderType = "constant")
           },
           valid =
           {
               y <- rip.filter(x, k, flip = flip,
                               anchor = anchor,
                               borderType = "constant")
               ## floor((n-1)/2) : N - floor(n/2)
               keep.row <- seq(1 + floor((d[1]-1) / 2), nrow(x) - floor(d[1]/2))
               keep.col <- seq(1 + floor((d[2]-1) / 2), ncol(x) - floor(d[2]/2))
               as.rip(y[keep.row, keep.col])
           })
}


### Discrete Fourier Transform and related utilities

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param x \code{"rip"} object to transform
##'
##' @param inverse Logical flag. If \code{TRUE}, inverse DFT is
##'     computed, assuming real output.
##'
##' @param scale Logical flag. If \code{TRUE}, divides transformed
##'     values by number of elements.
##'
##' @param rowwise Logical flag. If \code{TRUE}, computes DFT rowwise.
##'
##' @param complex Logical flag. If \code{TRUE}, stores output as
##'     complex numbers, otherwise imaginary part is stored in
##'     subdiagonal exploiting symmetry. Ignored when
##'     \code{inverse=TRUE}, in which case it is inferred from the
##'     input.
##'
##' @param nonzerorows Integer. Assumes that only first nonzerorows
##'     have non-zero data. When 'inverse=TRUE', only first
##'     nonzerorows rows in the output will have non-zero data.
##'
##' @return 
rip.dft <-
    function(x, pad = NULL, inverse = FALSE, scale = inverse,
             rowwise = FALSE, complex = TRUE,
             nonzerorows = 0)
{
    x <- as.rip(x)
    if (!is.null(pad)) x <- rip.pad(x, pad = pad, value = 0, borderType = "constant")
    eflags <- rip.cv$enums$DftFlags
    flags <- if (inverse) bitwOr(eflags["DFT_INVERSE"], eflags["DFT_REAL_OUTPUT"])
             else 0L
    if (scale) flags <- bitwOr(flags, eflags["DFT_SCALE"])
    if (rowwise) flags <- bitwOr(flags, eflags["DFT_ROWS"])
    if (inverse && !missing(complex))
        warning("Explicit specification of 'complex' flag ignored.")
    if (inverse) complex <- is.complex(x)
    if (complex) flags <- bitwOr(flags, eflags["DFT_COMPLEX_OUTPUT"])
    if (inverse && complex)
    {
        X <- matrix(0, nrow = nrow(x), ncol = 2L * ncol(x))
        X[, c(TRUE, FALSE)] <- Re(x)
        X[, c(FALSE, TRUE)] <- Im(x)
        x <- as.rip(X, channel = 2)
    }
    d <- rip.cv$transforms$dft(x, flags, nonzerorows)
    if (!inverse && complex)
    {
        d <-
            as.rip(matrix(complex(real = d[, c(TRUE, FALSE), drop = FALSE], 
                                  imaginary = d[, c(FALSE, TRUE), drop = FALSE]),
                          nrow = nrow(d)))
    }
    d
}

rip.ndft <- function(x, pad = NULL, inverse = FALSE)
    ## 'normalized' version with fewer options, to avoid scaling
    ## issues
{
    L <- if (is.null(pad)) prod(dim(x)) else prod(pad)
    f <- if (inverse) sqrt(L) else 1/sqrt(L)
    rip.dft(x, pad = pad, inverse = inverse) * f
}


rip.dct <- function(x, pad = NULL, inverse = FALSE, rowwise = FALSE)
{
    x <- as.rip(x)
    if (!is.null(pad)) x <- rip.pad(x, pad = pad, value = 0, borderType = "constant")
    eflags <- rip.cv$enums$DftFlags
    flags <- 0L
    if (inverse) flags <- bitwOr(flags, eflags["DCT_INVERSE"])
    if (rowwise) flags <- bitwOr(flags, eflags["DCT_ROWS"])
    rip.cv$transforms$dct(x, flags)
}


## Color to grayscale

rip.desaturate <- function(x, method = c("simple", "decolor", "convert"))
{
    x <- as.rip(x)
    nc <- nchannel(x)
    if (nc == 1) return(x)
    if (nc != 3) stop("Input image must have exactly 3 channels")
    ## FIXME: handle alpha channel?
    method <- match.arg(method)
    switch(method,
           simple =
           {
               a <- as.array(x)
               as.rip(0.299 * a[,,1] + 0.587 * a[,,2]+ 0.114 * a[,,3])
           },
           decolor = rip.cv$photo$decolor(x),
           convert =
           {
               col2gray <- rip.cv$enums$ColorConversionCodes["COLOR_BGR2GRAY"]
               rip.cv$imgproc$cvtColor(x, col2gray)
           })
}

## Miscellaneous utilites

rip.flip <- function(k, cv = nchannel(k) > 1)
{
    k <- as.rip(k)
    if (cv)
        rip.cv$core$flip(k)
    else
    {
        if (nchannel(k) != 1) stop("Input must have a single channel")
        k[] <- k[rev(seq_len(nrow(k))), rev(seq_len(ncol(k)))]
        k
    }
}

rip.shift <- function(x, inverse = FALSE,
                      orow = 1 + ceiling(nrow(x)/2), ocol = 1 + ceiling(ncol(x)/2))
{
    x <- as.rip(x)
    if (nchannel(x) != 1) stop("Input must have a single channel")
    ## The [orow, ocol] entry is shifted to the 'origin' [1,1]
    if (nrow(x) > 2 && ncol(x) > 2) {
        if (!inverse)
            x[] <- x[c(orow:nrow(x), 1:(orow - 1)),
                     c(ocol:ncol(x), 1:(ocol - 1))]
        else {
            orow <- 1 + floor(nrow(x)/2)
            ocol <- 1 + floor(ncol(x)/2)
            x[] <- x[c(orow:nrow(x), 1:(orow - 1)),
                     c(ocol:ncol(x), 1:(ocol - 1))]
        }
    }
    x
}


