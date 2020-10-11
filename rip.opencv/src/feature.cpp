#include "common.h"

using namespace Rcpp;
using namespace cv;

/* Canny Edge Detector */
Rcpp::NumericMatrix
cv_Canny(Rcpp::NumericMatrix imgMat, 
	 double threshold1, double threshold2,
	 int appertureSize, bool L2Gradient) 
{
    cv::Mat M = RCPP2CV(imgMat, 0);
    cv::Mat outImg;
    cv::Canny(M, outImg, threshold1, threshold2, appertureSize, L2Gradient);
    return (CV2RCPP(outImg));
}

/* Corner Eigen Values and Vectors */
/* imgMat has to be single channel 8-bit or floating-point image */
Rcpp::NumericMatrix
cv_cornerEigenValsAndVecs(Rcpp::NumericMatrix imgMat,
			  int blockSize, int kSize)
{
    cv::Mat M = RCPP2CV(imgMat, 0);
    cv::Mat outImg;
    cv::cornerEigenValsAndVecs(M, outImg, blockSize, kSize, BORDER_DEFAULT);
    return (CV2RCPP(outImg));
}

/* Harris Corner Detection */
Rcpp::NumericMatrix
cv_cornerHarris(Rcpp::NumericMatrix imgMat,
		int blockSize, int kSize, double k)
{
    cv::Mat M = RCPP2CV(imgMat, 0);
    cv::Mat outImg;
    cv::cornerHarris(M, outImg, blockSize, kSize, k);
    return (CV2RCPP(outImg));
}


/* Minimum eigen values of gradient matrix for corner detection*/
Rcpp::NumericMatrix
cv_cornerMinEigenVal(Rcpp::NumericMatrix imgMat, int blockSize, int ksize)
{
    cv::Mat M = RCPP2CV(imgMat, 0);
    cv::Mat outImg(M.rows, M.cols, CV_32FC1);
    cv::cornerMinEigenVal(M, outImg, blockSize, ksize);
    return (CV2RCPP(outImg));
}


/* Refines the corner location */

RCPP_MODULE(feature) 
{
    function("Canny", &cv_Canny, "Canny Edge Detection");
    function("cornerEigenValsAndVecs", &cv_cornerEigenValsAndVecs, "Calculate eigenvalues and eigenvectors of image blocks for corner detection");
    function("cornerHarris", &cv_cornerHarris, "Harris Corner Detection");
    function("cornerMinEigenVal", &cv_cornerMinEigenVal, "Minimum eigen values of gradient matrices");
}

