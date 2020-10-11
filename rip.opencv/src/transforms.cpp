#include "common.h"

// FIXME: add dct and idct ?

using namespace Rcpp;

int cv_getOptimalDFTSize(int vecSize)
{
    return cv::getOptimalDFTSize(vecSize);
}

Rcpp::NumericMatrix cv_dft(Rcpp::NumericMatrix x, int flags = 0, int nonzerorows = 0)
{
    // flags: OR of 1  = DFT_INVERSE,
    //              2  = DFT_SCALE (divide by number of elements)
    //              4  = DFT_ROWS  (transforms row-wise)
    //              16 = DFT_COMPLEX_OUTPUT (default is to store
    //                   imaginary part in subdiagonal exploiting symmetry)
    //              32 = DFT_REAL_OUTPUT (for inverse, assumes symmetry, produced real)
    // nonzerorows: only first nonzerorows have non-zero data (output
    //              when inverting)
    cv::Mat out, f, M = RCPP2CV(x, 5);
    M.convertTo(f, CV_32F);
    cv::dft(f, out, flags, nonzerorows);
    return CV2RCPP(out);
}

Rcpp::NumericMatrix cv_idft(Rcpp::NumericMatrix imgMat, int flags, int nonzerorows)
{
    cv::Mat out, f, M = RCPP2CV(imgMat, 5);
    M.convertTo(f, CV_32F);
    idft(f, out, flags);
    return (CV2RCPP(out));
}

// useful mostly when using DFT_REAL_OUTPUT in dft()
Rcpp::NumericMatrix cv_mulSpectrums(Rcpp::NumericMatrix x1, 
				    Rcpp::NumericMatrix x2, bool conj)
{
    // conj: whether to conjugate x2 before mult (gives correlation
    //        rather than convolution)
    cv::Mat M1 = RCPP2CV(x1, 5);
    cv::Mat M2 = RCPP2CV(x2, 5);
    cv::Mat out;
    cv::mulSpectrums(M1, M2, out, 0, conj);
    return CV2RCPP(out);
}

RCPP_MODULE(transforms) 
{
    function("dft", &cv_dft, "Discrete Fourier Transformation");
    function("idft", &cv_idft, "Inverse DFT");
    function("mulSpectrums", &cv_mulSpectrums, "Multiply DFT Spectrums");
    function("getOptimalDFTSize", &cv_getOptimalDFTSize, "get optimal DFT size for Fourier Transformation");
}

