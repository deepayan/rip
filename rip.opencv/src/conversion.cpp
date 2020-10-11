#include <Rcpp.h>
#include <opencv2/opencv.hpp>

#include "conversion.h"

using namespace Rcpp;
using namespace cv;

/* helper functions */
Rcpp::NumericMatrix convertUCMat_CV2RCPP(cv::Mat);
Rcpp::NumericMatrix convertCMat_CV2RCPP(cv::Mat);
Rcpp::NumericMatrix convertUSMat_CV2RCPP(cv::Mat);
Rcpp::NumericMatrix convertSMat_CV2RCPP(cv::Mat);
Rcpp::NumericMatrix convertIMat_CV2RCPP(cv::Mat);
Rcpp::NumericMatrix convertFMat_CV2RCPP(cv::Mat);
Rcpp::NumericMatrix convertDMat_CV2RCPP(cv::Mat);

void convertUCMat_RCPP2CV(Rcpp::NumericMatrix, cv::Mat, int, int);
void convertFMat_RCPP2CV(Rcpp::NumericMatrix, cv::Mat, int, int);
void convertDMat_RCPP2CV(Rcpp::NumericMatrix, cv::Mat, int, int);


/* Conversion of OpenCV matrix to RCPP NumericMatrix (with some attributes) */
Rcpp::NumericMatrix CV2RCPP(cv::Mat M)
{
    Rcpp::NumericMatrix ans;
    switch (M.depth()) {
    case 0: 
	ans = convertUCMat_CV2RCPP(M);
	break;
    case 5:
	ans = convertFMat_CV2RCPP(M);
	break;
    case 6:
	ans = convertDMat_CV2RCPP(M);
	break;
    case 1:
	ans = convertCMat_CV2RCPP(M);
	break;
    case 2:
	ans = convertUSMat_CV2RCPP(M);
	break;
    case 3:
	ans = convertSMat_CV2RCPP(M);
	break;
    case 4:
	ans = convertIMat_CV2RCPP(M);
	break;
    default:
	::Rf_error("Cannot convert");
    }
    ans.attr("cvdim") = NumericVector::create(Rcpp::Named("nrow") = M.rows,
					      Rcpp::Named("ncol") = M.cols,
					      Rcpp::Named("nchannel") = M.channels());
    ans.attr("class") = StringVector::create("rip", "matrix");
    return ans;
}

/*********************************************/
/**** This part includes separate functions that  
      convert OpenCV matrix of different type to NumericMatrix ****/
/********************************************/

/* Converting OpenCV uchar Mat to NumericMatrix */
/* image depth == 0 */
/***********************************/

Rcpp::NumericMatrix convertUCMat_CV2RCPP(cv::Mat M)
{
    NumericMatrix x(M.rows, M.cols * M.channels());
    for (int i = 0; i < M.rows; i++)
	for (int j = 0; j < M.cols * M.channels(); j++) {
	    x(i, j) = M.at<uchar>(i,j);
	}
    return x;
}


/* Converting OpenCV signed char Mat to NumericMatrix */
/* image depth == 1 */
/***********************************/

Rcpp::NumericMatrix convertCMat_CV2RCPP(cv::Mat M)
{
    NumericMatrix x(M.rows, M.cols * M.channels());
    for (int i = 0; i < M.rows; i++)
	for (int j = 0; j < M.cols * M.channels(); j++) {
	    x(i, j) = M.at<char>(i,j);
	}
    return x;
}



/* Converting OpenCV unsigned int Mat to NumericMatrix */
/******** image depth == 2 *********/
/***********************************/

Rcpp::NumericMatrix convertUSMat_CV2RCPP(cv::Mat M)
{
    NumericMatrix x(M.rows, M.cols * M.channels());
    for (int i = 0; i < M.rows; i++)
	for (int j = 0; j < M.cols * M.channels(); j++) {
	    x(i, j) = M.at<ushort>(i,j);
	}
    return x;
}

/* Converting OpenCV signed short int Mat to NumericMatrix */
/******** image depth == 3 *********/
/***********************************/

Rcpp::NumericMatrix convertSMat_CV2RCPP(cv::Mat M)
{
    NumericMatrix x(M.rows, M.cols * M.channels());
    for (int i = 0; i < M.rows; i++)
	for (int j = 0; j < M.cols * M.channels(); j++) {
	    x(i, j) = M.at<short>(i,j);
	}
    return x;
}


/* Converting OpenCV signed int Mat to NumericMatrix */
/******** image depth == 4 *********/
/***********************************/

Rcpp::NumericMatrix convertIMat_CV2RCPP(cv::Mat M)
{
    NumericMatrix x(M.rows, M.cols * M.channels());
    for (int i = 0; i < M.rows; i++)
	for (int j = 0; j < M.cols * M.channels(); j++) {
	    x(i, j) = M.at<int>(i,j);
	}
    return x;
}


/* Converting OpenCV float Mat to NumericMatrix */
/******** image depth == 5 *********/
/***********************************/

Rcpp::NumericMatrix convertFMat_CV2RCPP(cv::Mat M)
{
    NumericMatrix x(M.rows, M.cols * M.channels());
    for (int i = 0; i < M.rows; i++)
	for (int j = 0; j < M.cols * M.channels(); j++) {
	    x(i, j) = M.at<float>(i,j);
	}
    return x;
}


/* Converting OpenCV double Mat to NumericMatrix */
/******** image depth == 6 *********/
/***********************************/
Rcpp::NumericMatrix convertDMat_CV2RCPP(cv::Mat M)
{
    NumericMatrix x(M.rows, M.cols * M.channels());
    for (int i = 0; i < M.rows; i++)
	for (int j = 0; j < M.cols * M.channels(); j++) {
	    x(i, j) = M.at<double>(i,j);
	}
    return x;
}


/* Conversion of RCPP NumericMatrix containing allowed type of data to
   OpenCV matrix */


// depth: depth is the depth of the data in the Rcpp matrix. Usually
// we store it as Numericmatrix (double) as long as we are inside R,
// which corresponds to depth=5 in OpenCV. But some functions in
// OpenCV are not valid for depth=5, so we might need to use a
// different depth in those cases.

cv::Mat RCPP2CV(Rcpp::NumericMatrix x, int depth = 5)
{
    std::vector<int> dims = as< std::vector<int> > (x.attr("cvdim"));
    int nRow = dims[0];
    int nCol = dims[1];
    int nChannel = dims[2];
    int effCol = nCol * nChannel;
    int type = CV_MAKETYPE(depth, nChannel);
    cv::Mat M(nRow, nCol, type); // FIXME: What is type when we know depth and nChannels?
    if (depth == 0)
	convertUCMat_RCPP2CV(x, M, nRow, effCol);
    else if (depth == 5)
	convertFMat_RCPP2CV(x, M, nRow, effCol);
    else if (depth == 6)
	convertDMat_RCPP2CV(x, M, nRow, effCol);
    else 
	stop("Cannot convert matrix to cv::Mat");
    return M;
}

/* Converting NumericMatrix of various types to OpenCV Mat*/

void convertUCMat_RCPP2CV(Rcpp::NumericMatrix lx, cv::Mat M, int nRow, int effCol)
{
    int i, j;
    for(i = 0; i < nRow; i++) {
	for(j = 0; j < effCol; j++) {
	    M.at<uchar>(i, j) = lx(i, j);
	}
    }
}

void convertFMat_RCPP2CV(Rcpp::NumericMatrix lx, cv::Mat M, int nRow, int effCol)
{
    int i, j;
    for(i = 0; i < nRow; i++) {
	for(j = 0; j < effCol; j++) {
	    M.at<float>(i, j) = lx(i, j);
	}
    }
}

void convertDMat_RCPP2CV(Rcpp::NumericMatrix lx, cv::Mat M, int nRow, int effCol)
{
    int i, j;
    for(i = 0; i < nRow; i++) {
	for(j = 0; j < effCol; j++) {
	    M.at<double>(i, j) = lx(i, j);
	}
    }
}


