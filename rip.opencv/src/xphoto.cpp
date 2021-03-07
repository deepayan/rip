#include "common.h"

#include <opencv2/xphoto.hpp>

using namespace cv;
using namespace Rcpp;


// Functions from cv::xphoto

Rcpp::NumericMatrix
cv_xphoto_inpaint(Rcpp::NumericMatrix imgMat, Rcpp::NumericMatrix maskMat,
		  int algorithmType)
{
    // maskMat == zero => missing (different from cv::inpaint)
    cv::Mat outImg, M = RCPP2CV(imgMat, 0), mask = RCPP2CV(maskMat, 0);
    cv::xphoto::inpaint(M, mask, outImg, algorithmType);
    return CV2RCPP(outImg);
}


// bm3dDenoising() is not available unless OPENCV_ENABLE_NONFREE is
// set when building OpenCV. This is usually not the case in binary
// distributions such as Debian.

// Fortunately, this only gives a run-time error, not give a
// compile-time error, so not bothering to figure out whether it is
// available or not.


Rcpp::NumericMatrix
cv_xphoto_bm3dDenoising(Rcpp::NumericMatrix imgMat,
			double h,
			int templateWindowSize,
			int searchWindowSize)
{
    cv::Mat outImg, M = RCPP2CV(imgMat, 0);
    cv::xphoto::bm3dDenoising(M, outImg, (float) h,
			      templateWindowSize, searchWindowSize);
 			       // float    h = 1,
			       // int    templateWindowSize = 4,
			       // int    searchWindowSize = 16,
			       // int    blockMatchingStep1 = 2500,
			       // int    blockMatchingStep2 = 400,
			       // int    groupSize = 8,
			       // int    slidingStep = 1,
			       // float    beta = 2.0f,
			       // int    normType = cv::NORM_L2,
			       // int    step = cv::xphoto::BM3D_STEPALL,
			       // int    transformType = cv::xphoto::HAAR 
    return CV2RCPP(outImg);
}


Rcpp::NumericMatrix
cv_xphoto_dctDenoising(Rcpp::NumericMatrix imgMat,
		       double sigma = 1.0, int psize = 16)
{
    cv::Mat outImg, M = RCPP2CV(imgMat, 0);
    cv::xphoto::dctDenoising(M, outImg, sigma, psize);
    return CV2RCPP(outImg);
}


/*

void cv::xphoto::inpaint(const Mat &src,
			 const Mat &mask,
			 Mat & 	dst,
			 const int algorithmType 
			 )

mask: non-zero pixels indicate valid image area, zero pixels indicate area to be inpainted

algorithmType = INPAINT_SHIFTMAP: Any number of channels from 1 to
4. In case of 3- and 4-channels images the function expect them in
CIELab colorspace or similar one, where first color component shows
intensity, while second and third shows colors. Nonetheless you can
try any colorspaces.

algorithmType = INPAINT_FSR_BEST or INPAINT_FSR_FAST: 1-channel
grayscale or 3-channel BGR image.

void cv::xphoto::bm3dDenoising(InputArray    src,
			       OutputArray    dst,
			       float    h = 1,
			       int    templateWindowSize = 4,
			       int    searchWindowSize = 16,
			       int    blockMatchingStep1 = 2500,
			       int    blockMatchingStep2 = 400,
			       int    groupSize = 8,
			       int    slidingStep = 1,
			       float    beta = 2.0f,
			       int    normType = cv::NORM_L2,
			       int    step = cv::xphoto::BM3D_STEPALL,
			       int    transformType = cv::xphoto::HAAR 
			       )   

void cv::xphoto::dctDenoising(const Mat &src,
			      Mat &    dst,
			      const double    sigma,
			      const int    psize = 16 
			      )   

void cv::xphoto::oilPainting(InputArray 	src,
			     OutputArray 	dst,
			     int 	size,
			     int 	dynRatio,
			     int 	code 
			     )	

*/


RCPP_MODULE(xphoto)
{
    function("inpaint", &cv_xphoto_inpaint,
	     List::create(_["x"], _["mask"],
			  _["algorithmType"] = 0),
	     "Impute missing pixels (xphoto).");
    function("bm3dDenoising", &cv_xphoto_bm3dDenoising,
	     List::create(_["x"], _["h"] = 1.0,
			  _["templateWindowSize"] = 4,
			  _["searchWindowSize"] = 16),
	     "BM3D Denoising (xphoto). May not be available unless built with OPENCV_ENABLE_NONFREE enabled");
    function("dctDenoising", &cv_xphoto_dctDenoising,
	     List::create(_["x"], _["sigma"] = 1.0, _["psize"] = 16),
	     "DCT Denoising (xphoto).");
}

