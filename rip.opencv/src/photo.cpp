#include "common.h"

using namespace cv;
using namespace Rcpp;

Rcpp::NumericMatrix
cv_fastNlMeansDenoising(Rcpp::NumericMatrix imgMat, float h = 3.0, 
			int templateWindowSize = 7, int searchWindowSize = 21)
{
    cv::Mat outImg, M = RCPP2CV(imgMat, 0);
    cv::fastNlMeansDenoising(M, outImg, h, templateWindowSize, searchWindowSize);
    return CV2RCPP(outImg);
}


Rcpp::NumericMatrix
cv_fastNlMeansDenoisingColored(Rcpp::NumericMatrix imgMat, float h, 
			       int templateWindowSize, int searchWindowSize)
{
    cv::Mat outImg, M = RCPP2CV(imgMat, 0);
    cv::fastNlMeansDenoisingColored(M, outImg, h, templateWindowSize, searchWindowSize);
    return CV2RCPP(outImg);
}


Rcpp::NumericMatrix
cv_inpaint(Rcpp::NumericMatrix imgMat, Rcpp::NumericMatrix maskMat, 
	   double inpaintRadius = 5.0, int flags = cv::INPAINT_TELEA)
{
    // maskMat should be 8-bit single channel image (check?).
    // non-zero => missing
    cv::Mat outImg, M = RCPP2CV(imgMat, 0), mask = RCPP2CV(maskMat, 0);
    cv::inpaint(M, mask, outImg, inpaintRadius, flags);
    return CV2RCPP(outImg);
}

Rcpp::NumericMatrix
cv_decolor(Rcpp::NumericMatrix imgMat, int color = 0)
{
    cv::Mat outImg, M = RCPP2CV(imgMat, 0), outImg_color;
    // not clear what the third argument is, but it's a color image
    // different from the input
    cv::decolor(M, outImg, outImg_color);
    if (color) return CV2RCPP(outImg_color);
    else return CV2RCPP(outImg);
}

Rcpp::NumericMatrix
cv_pencilSketch(Rcpp::NumericMatrix imgMat, int color = 0,
		float sigma_s = 60,
		float sigma_r = 0.07f,
		float shade_factor = 0.02f)
{
    cv::Mat outImg, outImg_color, M = RCPP2CV(imgMat, 0);
    cv::pencilSketch(M, outImg, outImg_color, sigma_s, sigma_r, shade_factor);
    if (color) return CV2RCPP(outImg_color);
    else return CV2RCPP(outImg);
}

Rcpp::NumericMatrix
cv_stylization(Rcpp::NumericMatrix imgMat,
	       float sigma_s = 60.0,
	       float sigma_r = 0.45f)
{
    cv::Mat outImg, M = RCPP2CV(imgMat, 0);
    cv::stylization(M, outImg, sigma_s, sigma_r);
    return CV2RCPP(outImg);
}

Rcpp::NumericMatrix
cv_edgePreservingFilter(Rcpp::NumericMatrix imgMat, int flags,
			float sigma_s = 60.0,
			float sigma_r = 0.4f)
{
    cv::Mat outImg, M = RCPP2CV(imgMat, 0);
    cv::edgePreservingFilter(M, outImg, flags, sigma_s, sigma_r);
    return CV2RCPP(outImg);
}



RCPP_MODULE(photo)
{
    function("fastNlMeansDenoising", &cv_fastNlMeansDenoising,
	     List::create(_["x"], _["h"] = 3.0,
			  _["templateWindowSize"] = 7,
			  _["searchWindowSize"] = 21),
	     "NL means denoising.");
    function("fastNlMeansDenoisingColored", &cv_fastNlMeansDenoisingColored,
	     List::create(_["x"], _["h"] = 3.0,
			  _["templateWindowSize"] = 7,
			  _["searchWindowSize"] = 21),
	     "NL means denoising (colored).");
    function("inpaint", &cv_inpaint,
	     List::create(_["x"], _["mask"],
			  _["radius"] = 5, _["flags"] = 1),
	     "Impute missing pixels.");
    function("decolor", &cv_decolor,
	     List::create(_["x"], _["color"] = 0),
	     "Contrast-preserving transformation of a color image to a grayscale image.");
    function("pencilSketch", &cv_pencilSketch,
	     List::create(_["x"], _["color"] = 0,
			  _["sigma_s"] = 60.0, _["sigma_r"] = 0.07,
			  _["shade_factor"] = 0.02),
	     "Pencil-like non-photorealistic line drawing.");
    function("stylization", &cv_stylization,
	     List::create(_["x"], _["sigma_s"] = 60.0, _["sigma_r"] = 0.45),
	     "Non-photorealistic stylization.");
    function("edgePreservingFilter", &cv_edgePreservingFilter,
	     List::create(_["x"], _["flags"] = 1, _["sigma_s"] = 60.0, _["sigma_r"] = 0.4),
	     "Edge preserving filter.");
}
