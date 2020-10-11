#ifndef __COMMON_H__
#define __COMMON_H__

#include <RcppCommon.h>
#include <opencv2/opencv.hpp>

// FIXME: Figure out how to add custom as() and wrap() methods for
// cv::Mat <-> SEXP conversion. This should make it easier to expose
// some opencv functions, although most will need explicit wrappers
// because they use in-place modifications (return type is void). See
// https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-extending.pdf

#include <Rcpp.h>
#include "conversion.h"

#endif
