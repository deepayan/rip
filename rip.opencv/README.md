
# rip.opencv

The goal of this package is to provide access to selected routines in
the [OpenCV](http://opencv.org/) computer vision library from R. In
some ways it is similar to the
[opencv](https://github.com/ropensci/opencv) R package, which may be
more useful depending on your use case.

The package provides a low-level interface to several OpenCV routines
through Rcpp modules. High-level R wrappers are available for some
routines, and more may be added in future.

For more details, see the package vignette or its longer version
[here](https://deepayan.github.io/rip/opencv-intro.html).

## Installation

The package is not available on CRAN yet (no man pages yet!), but can
be installed from Github. Installation requires OpenCV development
libraries. Thanks to a configure script adapted from the `opencv`
package, installing from source should work seamlessly on most
platforms, or at least give hints on how to proceed. For example:

```
library(remotes) # install first if necessary
remotes::install_github("deepayan/rip/rip.opencv")
```

## Overview

Rather than exposing OpenCV classes and methods directly through
external pointers, the package uses standard R objects to represent
corresponding OpenCV objects, and explicitly converts between the two
forms as necessary. The idea behind this design is that most
operations will be performed in R, and OpenCV is used to simply make
an additional suite of operations available.

The package defines a "rip" class in R that parallels the `cv::Mat`
class in OpenCV, essentially an analogue of three-dimensional arrays
in R. Conversion functions written in C++ map "rip" objects to and
from `cv::Mat` objects as necessary. 

Wrappers to OpenCV functions are exposed as `Rcpp` modules,
specifically, as the environment `rip.cv`:

```r
ls(rip.cv)
```
```
[1] "enums"      "feature"    "filter"     "geomtrans"  "IO"        
[6] "misc"       "photo"      "transforms" "utils"     
```

Each of these modules consist of related functions (except `enums`,
which is a list of mnemonic codes), for example

```r
rip.cv$IO # after being initialized
```
```
Rcpp module 'IO' 
        3 functions: 
         imread : 2 arguments
        imwrite : 2 arguments
         vfread : 2 arguments

        0 classes : 
```

More routines may be added in future; an incomplete list of potential
candidates are available [here](opencv-functions.md).

Some high-level convenience functions are also available to cover
common use cases.


