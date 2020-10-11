#include "common.h"

using namespace Rcpp;

/* Border Padding */
Rcpp::NumericMatrix cv_copyMakeBorder(Rcpp::NumericMatrix imgMat, 
				      int top, int bottom, 
				      int left, int right,
				      int borderType,
				      double value)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat outImg;
    cv::copyMakeBorder(M, outImg, top, bottom, left, right,
		       borderType, cv::Scalar(value));
    return CV2RCPP(outImg);
}

// PSNR 
double cv_PSNR(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y)
{
    cv::Mat M1 = RCPP2CV(x, 5), M2 = RCPP2CV(y, 5);
    return PSNR(M1, M2);
}


Rcpp::NumericMatrix cv_flip(Rcpp::NumericMatrix imgMat, int code)
{
    cv::Mat out, M = RCPP2CV(imgMat, 5);
    flip(M, out, code);
    return CV2RCPP(out);
}

// FIXME: add rotateFlags enum
Rcpp::NumericMatrix cv_rotate(Rcpp::NumericMatrix imgMat, int code)
{
    cv::Mat out, M = RCPP2CV(imgMat, 5);
    rotate(M, out, code);
    return CV2RCPP(out);
}

RCPP_MODULE(core) 
{
    function("copyMakeBorder", &cv_copyMakeBorder,
	     List::create(_["x"],
			  _["top"] = 0, _["bottom"] = 0,
			  _["left"] = 0, _["right"] = 0,
			  _["borderType"] = 0, _["value"] = 0.0),
	     "Pad a matrix by adding borders.");
    function("PSNR", &cv_PSNR, List::create(_["x"], _["y"]),
	     "Peak Signal-to-Noise Ratio (PSNR)");
    function("flip", &cv_flip, List::create(_["x"], _["code"] = -1L),
	     "Flip (reverse) rows (code = 0), columns (code = 1) or both (code = -1)");
    function("rotate", &cv_rotate, List::create(_["x"], _["code"] = 0L),
	     "Rotate image clockwise: code 0=90, 1=180, 2=270");
}
