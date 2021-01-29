
## See for parameter meanings https://en.wikipedia.org/wiki/Structural_similarity

rip.ssim <- function(x, y, k1 = 0.01, k2 = 0.03, L = 255,
                     ksize = 11L, sigma = 1.5, 
                     scale = diff(range(x)) < 1.01,
                     summarize = mean)
{
    if (scale)
    {
        x[] <- 255 * x
        y[] <- 255 * y
    }
    c1 <- (k1 * L)^2
    c2 <- (k2 * L)^2
    mu.x <- rip.blur(x, "gaussian", c(ksize, ksize), sigma, sigma)
    mu.y <- rip.blur(y, "gaussian", c(ksize, ksize), sigma, sigma)
    mu.xy <- mu.x * mu.y
    var.x <- rip.blur(x^2, "gaussian", c(ksize, ksize), sigma, sigma) - mu.x^2
    var.y <- rip.blur(y^2, "gaussian", c(ksize, ksize), sigma, sigma) - mu.y^2
    cov.xy <- rip.blur(x * y, "gaussian", c(ksize, ksize), sigma, sigma) - mu.xy
    n  <-  (2 * mu.xy + c1) * (2 * cov.xy + c2)
    d <- (mu.x^2 + mu.y^2 + c1) * (var.x + var.y + c2)
    ssim <- n / d
    if (is.null(summarize)) ssim else summarize(ssim)
}

