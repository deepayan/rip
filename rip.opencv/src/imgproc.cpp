#include "common.h"

using namespace Rcpp;
using namespace cv;


Rcpp::NumericMatrix
cv_getStructuringElement(int shape, Rcpp::IntegerVector ksize, Rcpp::IntegerVector anchor)
{
    // std::vector<int> anchor = as< std::vector<int> >(anchor_);
    return CV2RCPP(getStructuringElement(shape,
					 Size(ksize.at(0), ksize.at(1)),
					 Point(anchor.at(0), anchor.at(1))));
}

Rcpp::NumericMatrix
cv_dilate(Rcpp::NumericMatrix imgMat,
	  Rcpp::NumericMatrix kern, 
	  Rcpp::IntegerVector anchor_, int iteration)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat K = RCPP2CV(kern, 5);
    std::vector<int> anchor = as< std::vector<int> >(anchor_);
    cv::Mat outImg;
    cv::dilate(M, outImg, K, Point(anchor[0], anchor[1]));
    return CV2RCPP(outImg);
}

Rcpp::NumericMatrix
cv_erode(Rcpp::NumericMatrix imgMat,
	 Rcpp::NumericMatrix kern, 
	 Rcpp::IntegerVector anchor_, int iteration)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat K = RCPP2CV(kern, 5);
    std::vector<int> anchor = as< std::vector<int> >(anchor_);
    cv::Mat outImg;
    cv::erode(M, outImg, K, Point(anchor[0], anchor[1]));
    return CV2RCPP(outImg);
}


/* Resizing generic geometrical transformation to an image */
/* Interpolation methods are given by the following constants */
/****** INTER_NEAREST == 0 ******/
/****** INTER_LINEAR == 1 ******/
/****** INTER_CUBIC == 2 ******/
/****** INTER_AREA == 3 ******/
/****** INTER_LANCZOS4 == 4  ******/
/*******************************************************/



Rcpp::NumericMatrix
cv_resize(Rcpp::NumericMatrix imgMat, Rcpp::NumericVector dsize_, 
	  double fx, double fy, int interpolation)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    std::vector<int> dsize = as< std::vector<int> >(dsize_);
    cv::Mat outImg;
    cv::resize(M, outImg, Size(dsize.at(0), dsize.at(1)), fx, fy, interpolation);
    return CV2RCPP(outImg);
}


Rcpp::NumericMatrix
cv_warpAffine(Rcpp::NumericMatrix imgMat,
	      Rcpp::NumericMatrix tmat_,
	      Rcpp::NumericVector dsize_, 
	      int flags, int borderMode)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat tmat = RCPP2CV(tmat_, 5);
    std::vector<int> dsize = as< std::vector<int> >(dsize_);
    cv::Mat outImg;
    cv::warpAffine(M, outImg, tmat, Size(dsize.at(0), dsize.at(1)), flags, borderMode);
    return CV2RCPP(outImg);
}


Rcpp::NumericMatrix
cv_invertAffineTransform(Rcpp::NumericMatrix tmat_)
{
    cv::Mat tmat = RCPP2CV(tmat_, 5);
    cv::Mat itmat;
    cv::invertAffineTransform(tmat, itmat);
    return(CV2RCPP(itmat));
}

Rcpp::NumericMatrix
cv_warpPerspective(Rcpp::NumericMatrix imgMat,
		   Rcpp::NumericMatrix tmat_, 
		   Rcpp::NumericVector size_)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat tmat = RCPP2CV(tmat_, 5);
    std::vector<int> size = as< std::vector<int> >(size_);
    cv::Mat outImg;
    cv::warpPerspective(M, outImg, tmat, Size(size.at(0), size.at(1)));
    return CV2RCPP(outImg);
}


Rcpp::NumericMatrix
cv_getRectSubPix(Rcpp::NumericMatrix imgMat, Rcpp::NumericVector patchsize_, 
		 Rcpp::NumericVector center_)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    std::vector<int> patchsize = as< std::vector<int> >(patchsize_);
    std::vector<int> center = as< std::vector<int> >(center_);
    cv::Mat outImg;
    cv::Point2f centerP(center.at(0), center.at(1));
    cv::getRectSubPix(M, Size(patchsize.at(0), patchsize.at(1)), centerP, outImg, -1);
    return CV2RCPP(outImg);
}

Rcpp::NumericMatrix
cv_getRotationMatrix2D(Rcpp::NumericVector center_, double angle, double scale)
{
    std::vector<int> center = as< std::vector<int> >(center_);
    Point2f centerP(center.at(0), center.at(1));
    cv::Mat outMat = cv::getRotationMatrix2D(centerP, angle, scale);
    return CV2RCPP(outMat);
}

/* Upsampling an image */
Rcpp::NumericMatrix
cv_pyrUp(Rcpp::NumericMatrix imgMat,
	 Rcpp::NumericVector size_,
	 int borderType = BORDER_REFLECT_101)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat outImg;
    std::vector<int> size = as< std::vector<int> >(size_);
    cv::pyrUp(M, outImg, Size(size.at(0), size.at(1)), borderType);
    return CV2RCPP(outImg);
}


/* Downsampling an image */
Rcpp::NumericMatrix
cv_pyrDown(Rcpp::NumericMatrix imgMat,
	   Rcpp::NumericVector size_,
	   int borderType = BORDER_REFLECT_101)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat outImg;
    std::vector<int> size = as< std::vector<int> >(size_);
    cv::pyrDown(M, outImg, Size(size.at(0), size.at(1)), borderType);
    return CV2RCPP(outImg);
}


Rcpp::NumericMatrix
cv_pyrMeanShiftFiltering(Rcpp::NumericMatrix imgMat, double sp, double sr, int maxlevel)
{
    cv::Mat M = RCPP2CV(imgMat, 0);
    cv::Mat outImg;
    cv::pyrMeanShiftFiltering(M, outImg, sp, sr, maxlevel,
                              TermCriteria(TermCriteria::MAX_ITER+TermCriteria::EPS,5,1));
    return CV2RCPP(outImg);
}


Rcpp::NumericMatrix
cv_cvtColor(Rcpp::NumericMatrix imgMat, int code, int dstCn)
{
    cv::Mat outImg, M = RCPP2CV(imgMat, 0);
    cv::cvtColor(M, outImg, code, dstCn);
    return CV2RCPP(outImg);
}

Rcpp::NumericMatrix
cv_matchTemplate(Rcpp::NumericMatrix imgMat,
		 Rcpp::NumericMatrix temp, int method)
{
    cv::Mat M = RCPP2CV(imgMat, 5);
    cv::Mat outImg;
    cv::Mat T = RCPP2CV(temp, 5);
    cv::matchTemplate(M, T, outImg, method);
    return (CV2RCPP(outImg));
}


Rcpp::NumericMatrix cv_equalizeHist(Rcpp::NumericMatrix x)
{
    cv::Mat out, M = RCPP2CV(x, 5);
    cv::equalizeHist(M, out);
    return CV2RCPP(out);
}


RCPP_MODULE(imgproc) 
{
    function("getStructuringElement", &cv_getStructuringElement,
	     List::create(_["shape"] = 2L,
			  _["ksize"] = IntegerVector::create(3, 3),
			  _["anchor"] = IntegerVector::create(-1, -1)),
	     "Generate structuring element ");
    function("erode", &cv_erode, "Erodes an image");
    function("dilate", &cv_dilate, "Dilates an image");
    function("resize", &cv_resize, "Image resizing");
    function("warpAffine", &cv_warpAffine, "Affine warpping");
    function("warpPerspective", &cv_warpPerspective, "Perspective Warpping");
    function("getRotationMatrix2D", &cv_getRotationMatrix2D, "get rotation matrix");
    function("invertAffineTransform", &cv_invertAffineTransform, "inverts an affine transformation");
    function("getRectSubPix", &cv_getRectSubPix, "Retrives a pixel rectangle");
    function("pyrUp", &cv_pyrUp, "Pyramid upsample");
    function("pyrDown", &cv_pyrDown, "Pyramid downsample");
    function("pyrMeanShiftFiltering", &cv_pyrMeanShiftFiltering,
	     "Pyramid mean shift filtering");
    function("cvtColor", &cv_cvtColor,
	     List::create(_["x"], _["code"] = 6, _["dstCn"] = 0),
	     "Convert between colorspaces (default: BGR to grayscale).");
    function("matchTemplate", &cv_matchTemplate, "Template matching");
    function("equalizeHist", &cv_equalizeHist, "Histogram Equalization");
}


