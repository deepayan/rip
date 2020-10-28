# rip: Image Processing in R

This repository contains R packages intended to provide image
processing tools in R.

The `rip.opencv` package provides an interface to a subset of the
[OpenCV](http://opencv.org/) computer vision library. See [this
introduction](https://deepayan.github.io/rip/opencv-intro.html) for
details.

The `rip.recover` package provides tools for image reconstruction, and
currently implements a Bayesian approach for image deconvolution and
super-resolution using a "natural" gradient image prior. A short
introduction is available
[here](https://deepayan.github.io/rip/recover-intro.html).


## Installation

These packages are not available on CRAN yet, but can be installed
from Github. Installation of `rip.opencv` requires OpenCV development
libraries, and you will need a working toolchain to compile from
source. This should work seamlessly on most platforms, or at least
give hints on how to proceed. For example:

```
library(remotes) # install first if necessary
remotes::install_github("deepayan/rip/rip.opencv")
remotes::install_github("deepayan/rip/rip.recover")
```



