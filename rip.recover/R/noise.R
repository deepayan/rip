
##' Add noise to an image using a Poisson noise model 
##'
##' .. content for \details{} ..
##' @title
##' @param x Input image matrix, with values expected to be between 0
##'     and 255.
##' @param S Amount of noise to be added, 'NA' or 'Inf' to skip. Each
##'     value 'v' in the image is independently replaced by a Possion
##'     variate with mean 'S * v', which is then divided by
##'     'S'. Higher values of 'S' thus mean less noise.
##' @param restrict Logical flag indicating whether to restrict the
##'     result to lie between 0 and 255.
##' @param round Logical flag indicating whether to round the result
##'     to the nearest integer.
##' @param rescale Logical flag indicating whether to scale the
##'     result, by dividing by 255.
##' @return
addPoissonNoise <- function(x, S = 1, restrict = TRUE,
                            round = TRUE, rescale = TRUE)
{
    ## x assumed to be values between 0 to 255
    if (max(x) < 2) warning("Image appears to be scaled: values range between ",
                            round(min(x), 3), " and ", round(max(x), 3))
    if (is.finite(S)) x[] <- rpois(length(x), lambda = S * x) / S
    if (restrict)
    {
        x[x < 0] <- 0
        x[x > 255] <- 255
    }
    if (round) x[] <- round(x)
    if (rescale) x[] <- x / 255
    x
}

