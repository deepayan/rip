#ifndef __CONVERSION_H__
#define __CONVERSION_H__

/* conversion functions R matrix <-> cv::Mat */

Rcpp::NumericMatrix CV2RCPP(cv::Mat);
cv::Mat RCPP2CV(Rcpp::NumericMatrix, int);

#endif
