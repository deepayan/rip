#include "common.h"

using namespace cv;
using namespace Rcpp;

// Border types: values=?
// * BORDER_REPLICATE:     aaaaaa|abcdefgh|hhhhhhh
// * BORDER_REFLECT:       fedcba|abcdefgh|hgfedcb
// * BORDER_REFLECT_101:   gfedcb|abcdefgh|gfedcba
// * BORDER_WRAP:          cdefgh|abcdefgh|abcdefg
// * BORDER_CONSTANT:      iiiiii|abcdefgh|iiiiiii  with some specified 'i'



/* Bilateral filter of an image */
Rcpp::NumericMatrix
cv_bilateralFilter(Rcpp::NumericMatrix x, int d,
		   double sigmaColor, double sigmaSpace,
		   int borderType = BORDER_REFLECT_101)
{
    cv::Mat out, M = RCPP2CV(x, 5);
    cv::bilateralFilter(M, out, d, sigmaColor, sigmaSpace, borderType);
    return CV2RCPP(out);
}

// adaptiveBilateralFilter(M, out, ksize, sigmaSpace, maxSigmaColor=20.0,
// 			anchor=Point(-1, -1), int borderType )
// ???


Rcpp::NumericMatrix
cv_blur(Rcpp::NumericMatrix x,
	Rcpp::NumericVector size_,
	Rcpp::NumericVector anchor_,
	int borderType = BORDER_REFLECT_101) 
{
    cv::Mat M = RCPP2CV(x, 5);
    cv::Mat out;
    std::vector<int> size = as< std::vector<int> >(size_);
    std::vector<int> anchor = as< std::vector<int> >(anchor_);
    cv::blur(M, out, Size(size.at(0), size.at(1)),
	     Point(anchor.at(0), anchor.at(1)), borderType);
    return CV2RCPP(out);
}

Rcpp::NumericMatrix
cv_GaussianBlur(Rcpp::NumericMatrix imgMat,
		Rcpp::NumericVector size_, 
		double sigmaX, double sigmaY,
		int borderType = BORDER_REFLECT_101) 
{
    std::vector<int> size = as< std::vector<int> >(size_);
    cv::Mat outImg, M = RCPP2CV(imgMat, 5);
    cv::GaussianBlur(M, outImg, Size(size.at(0), size.at(1)), sigmaX, sigmaY, borderType);
    return CV2RCPP(outImg);
}

Rcpp::NumericMatrix
cv_boxFilter(Rcpp::NumericMatrix imgMat,
	     Rcpp::NumericVector size_, Rcpp::NumericVector anchor_,
	     bool normalize,
	     int borderType = BORDER_REFLECT_101) 
{
    std::vector<int> size = as< std::vector<int> >(size_);
    std::vector<int> anchor = as< std::vector<int> >(anchor_);
    cv::Mat outImg, M = RCPP2CV(imgMat, 5);
    cv::boxFilter(M, outImg, -1,
		  cv::Size(size.at(0), size.at(1)),
		  cv::Point(anchor.at(0), anchor.at(1)),
		  normalize, borderType);
    return CV2RCPP(outImg);
}


/* Apply 2d filter to an image given a kerel */
Rcpp::NumericMatrix
cv_filter2D(Rcpp::NumericMatrix imgMat,
	    Rcpp::NumericMatrix kernel_,
	    IntegerVector anchor_, double delta,
	    int borderType = BORDER_REFLECT_101)
{
    cv::Mat out, M = RCPP2CV(imgMat, 5), kernel = RCPP2CV(kernel_, 5);
    std::vector<int> anchor = as< std::vector<int> >(anchor_);
    cv::filter2D(M, out, -1, kernel, cv::Point(anchor[0], anchor[1]), delta, borderType);
    return CV2RCPP(out);
}



/* Get derivative kernel */

Rcpp::List
cv_getDerivKernels(int dx, int dy, int ksize, bool normalize, int ktype)
{
    cv::Mat kx, ky;
    cv::getDerivKernels(kx, ky, dx, dy, ksize, normalize, ktype);
    return Rcpp::List::create(Rcpp::Named("x") = CV2RCPP(kx), 
			      Rcpp::Named("y") = CV2RCPP(ky));   
}

/* Get Gaussian Kernel */

Rcpp::NumericMatrix
cv_getGaussianKernel(int ksize, double sigma, int ktype)
{
    cv::Mat out = cv::getGaussianKernel(ksize, sigma, ktype);
    return CV2RCPP(out);
}


Rcpp::NumericMatrix
cv_getGaborKernel(Rcpp::NumericVector ksize_, double sigma, double theta,
		  double lambd, double gamma, double psi, int ktype)
{
    std::vector<int> ksize = as< std::vector<int> >(ksize_);
    cv::Mat out = cv::getGaborKernel(Size(ksize.at(0), ksize.at(1)),
				     sigma, theta, lambd, gamma, psi, ktype); 
    return CV2RCPP(out);
}

Rcpp::NumericMatrix
cv_medianBlur(Rcpp::NumericMatrix x, int ksize)
{
    cv::Mat out, M = RCPP2CV(x, 5);
    cv::medianBlur(M, out, ksize);
    return CV2RCPP(out);
}


Rcpp::NumericMatrix
cv_Scharr(Rcpp::NumericMatrix imgMat,
	  int dx, int dy, double scale = 1,
	  double delta = 0.0,
	  int bordertype = BORDER_REFLECT_101)
{
    cv::Mat out, M = RCPP2CV(imgMat, 5);
    cv::Scharr(M, out, -1, dx, dy, scale, delta);
    return CV2RCPP(out);
}




Rcpp::NumericMatrix
cv_Sobel(Rcpp::NumericMatrix imgMat, int dx, int dy, 
	 int ksize, double scale = 1,
	 double delta = 0.0,
	 int bordertype = BORDER_REFLECT_101)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat outImg;
    cv::Sobel(M, outImg, -1, dx, dy, ksize, scale, delta);
    return (CV2RCPP(outImg));
}

Rcpp::NumericMatrix
cv_sepFilter2D(Rcpp::NumericMatrix imgMat, 
	       Rcpp::NumericMatrix kernelX_, Rcpp::NumericMatrix kernelY_, 
	       Rcpp::NumericVector anchor_, double delta = 0.0,
	       int bordertype = BORDER_REFLECT_101)
{
    // applies kernelX row-wise, then kernelY column-wise
    cv::Mat
	out,
	M = RCPP2CV(imgMat, 5),
	kernelX = RCPP2CV(kernelX_, 5),
	kernelY = RCPP2CV(kernelY_, 5);
    std::vector<int> anchor = as< std::vector<int> >(anchor_);
    cv::sepFilter2D(M, out, -1, kernelX, kernelY, 
		    Point(anchor.at(0), anchor.at(1)), delta);
    return (CV2RCPP(out));
}

    
RCPP_MODULE(filter) 
{
    function("getGaussianKernel", &cv_getGaussianKernel, "Gaussian Kernel");
    function("getDerivKernels", &cv_getDerivKernels, "Derivative Kernel");
    function("getGaborKernel", &cv_getGaborKernel, "Gabor Kernel");
    function("medianBlur", &cv_medianBlur, "Median blur");
    function("GaussianBlur", &cv_GaussianBlur, "Gaussian Blurring");
    function("bilateralFilter", &cv_bilateralFilter, "Bilateral Filtering to an image");
    function("boxFilter", &cv_boxFilter, "Box Blurring");
    function("blur", &cv_blur, "Normalized Box Blurring");
    function("filter2D", &cv_filter2D, "Convolves an image with a kernel");
    function("sepFilter2D", &cv_sepFilter2D, "Applies separable linear filter to image.");
    function("Sobel", &cv_Sobel, "Sobel derivative");
    function("Scharr", &cv_Scharr, "Scharr derivative");
}


