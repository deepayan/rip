
#include "common.h"

using namespace cv;
using namespace Rcpp;

// FIXME TODO For better maintenance, generate these from some kind of
// more easily editable database, e.g. a CSV file


Rcpp::IntegerVector cv_enum_BorderTypes()
{
    IntegerVector
	e = IntegerVector::create(_["BORDER_CONSTANT"] = (int) BORDER_CONSTANT,
				  _["BORDER_REPLICATE"] = (int) BORDER_REPLICATE,
				  _["BORDER_REFLECT"] = (int) BORDER_REFLECT,
				  _["BORDER_WRAP"] = (int) BORDER_WRAP,
				  _["BORDER_REFLECT_101"] = (int) BORDER_REFLECT_101);
    return e;
}


Rcpp::IntegerVector cv_enum_DftFlags()
{
    IntegerVector
	e = IntegerVector::create(_["DCT_INVERSE"] = (int) DCT_INVERSE,
				  _["DCT_ROWS"] = (int) DCT_ROWS,
				  _["DFT_INVERSE"] = (int) DFT_INVERSE,
				  _["DFT_SCALE"] = (int) DFT_SCALE,
				  _["DFT_ROWS"] = (int) DFT_ROWS,
				  _["DFT_COMPLEX_OUTPUT"] = (int) DFT_COMPLEX_OUTPUT,
				  _["DFT_REAL_OUTPUT"] = (int) DFT_REAL_OUTPUT);
    return e;
}


// INTER_NEAREST - a nearest-neighbor interpolation
// INTER_LINEAR - a bilinear interpolation (used by default)
// INTER_AREA - resampling using pixel area relation. It may be a
//              preferred method for image decimation, as it gives
//              moire-free results. But when the image is zoomed, it
//              is similar to the INTER_NEAREST method.
// INTER_CUBIC - a bicubic interpolation over 4x4 pixel neighborhood
// INTER_LANCZOS4 - a Lanczos interpolation over 8x8 pixel neighborhood

Rcpp::IntegerVector cv_enum_InterpolationFlags()
{
    IntegerVector
	e = IntegerVector::create(_["INTER_NEAREST"] = (int) INTER_NEAREST,
				  _["INTER_LINEAR"] = (int) INTER_LINEAR,
				  _["INTER_CUBIC"] = (int) INTER_CUBIC,
				  _["INTER_AREA"] = (int) INTER_AREA,
				  _["INTER_LANCZOS4"] = (int) INTER_LANCZOS4,
				  _["INTER_MAX"] = (int) INTER_MAX,
				  _["WARP_FILL_OUTLIERS"] = (int) WARP_FILL_OUTLIERS,
				  _["WARP_INVERSE_MAP"] = (int) WARP_INVERSE_MAP);
    return e;
}


Rcpp::IntegerVector cv_enum_ImreadModes()
{
    // There are more, but not making them available
    IntegerVector
	e = IntegerVector::create(_["IMREAD_UNCHANGED"] = (int) IMREAD_UNCHANGED,
				  _["IMREAD_GRAYSCALE"] = (int) IMREAD_GRAYSCALE,
				  _["IMREAD_COLOR"] = (int) IMREAD_COLOR);
    return e;
}





Rcpp::List cv_enum_ColorConversionCodes()
{
    // There are some more very specialized conversion omitted here;
    // see opencv docs if needed
    IntegerVector
	e1 = IntegerVector::create(_["COLOR_BGR2BGRA"] = (int) COLOR_BGR2BGRA,
				   _["COLOR_RGB2RGBA"] = (int) COLOR_RGB2RGBA,
				   _["COLOR_BGRA2BGR"] = (int) COLOR_BGRA2BGR,
				   _["COLOR_RGBA2RGB"] = (int) COLOR_RGBA2RGB,
				   _["COLOR_BGR2RGBA"] = (int) COLOR_BGR2RGBA,
				   _["COLOR_RGB2BGRA"] = (int) COLOR_RGB2BGRA,
				   _["COLOR_RGBA2BGR"] = (int) COLOR_RGBA2BGR,
				   _["COLOR_BGRA2RGB"] = (int) COLOR_BGRA2RGB,
				   _["COLOR_BGR2RGB"] = (int) COLOR_BGR2RGB,
				   _["COLOR_RGB2BGR"] = (int) COLOR_RGB2BGR,
				   _["COLOR_BGRA2RGBA"] = (int) COLOR_BGRA2RGBA,
				   _["COLOR_RGBA2BGRA"] = (int) COLOR_RGBA2BGRA,
				   _["COLOR_BGR2GRAY"] = (int) COLOR_BGR2GRAY,
				   _["COLOR_RGB2GRAY"] = (int) COLOR_RGB2GRAY,
				   _["COLOR_GRAY2BGR"] = (int) COLOR_GRAY2BGR,
				   _["COLOR_GRAY2RGB"] = (int) COLOR_GRAY2RGB,
				   _["COLOR_GRAY2BGRA"] = (int) COLOR_GRAY2BGRA,
				   _["COLOR_GRAY2RGBA"] = (int) COLOR_GRAY2RGBA,
				   _["COLOR_BGRA2GRAY"] = (int) COLOR_BGRA2GRAY,
				   _["COLOR_RGBA2GRAY"] = (int) COLOR_RGBA2GRAY),
	e2 = IntegerVector::create(_["COLOR_BGR2XYZ"] = (int) COLOR_BGR2XYZ,
				   _["COLOR_RGB2XYZ"] = (int) COLOR_RGB2XYZ,
				   _["COLOR_XYZ2BGR"] = (int) COLOR_XYZ2BGR,
				   _["COLOR_XYZ2RGB"] = (int) COLOR_XYZ2RGB,
				   _["COLOR_BGR2YCrCb"] = (int) COLOR_BGR2YCrCb,
				   _["COLOR_RGB2YCrCb"] = (int) COLOR_RGB2YCrCb,
				   _["COLOR_YCrCb2BGR"] = (int) COLOR_YCrCb2BGR,
				   _["COLOR_YCrCb2RGB"] = (int) COLOR_YCrCb2RGB,
				   _["COLOR_BGR2HSV"] = (int) COLOR_BGR2HSV,
				   _["COLOR_RGB2HSV"] = (int) COLOR_RGB2HSV,
				   _["COLOR_BGR2Lab"] = (int) COLOR_BGR2Lab,
				   _["COLOR_RGB2Lab"] = (int) COLOR_RGB2Lab,
				   _["COLOR_BGR2Luv"] = (int) COLOR_BGR2Luv,
				   _["COLOR_RGB2Luv"] = (int) COLOR_RGB2Luv,
				   _["COLOR_BGR2HLS"] = (int) COLOR_BGR2HLS,
				   _["COLOR_RGB2HLS"] = (int) COLOR_RGB2HLS,
				   _["COLOR_HSV2BGR"] = (int) COLOR_HSV2BGR,
				   _["COLOR_HSV2RGB"] = (int) COLOR_HSV2RGB,
				   _["COLOR_Lab2BGR"] = (int) COLOR_Lab2BGR,
				   _["COLOR_Lab2RGB"] = (int) COLOR_Lab2RGB),
	e3 = IntegerVector::create(_["COLOR_Luv2BGR"] = (int) COLOR_Luv2BGR,
				   _["COLOR_Luv2RGB"] = (int) COLOR_Luv2RGB,
				   _["COLOR_HLS2BGR"] = (int) COLOR_HLS2BGR,
				   _["COLOR_HLS2RGB"] = (int) COLOR_HLS2RGB,
				   _["COLOR_LBGR2Lab"] = (int) COLOR_LBGR2Lab,
				   _["COLOR_LRGB2Lab"] = (int) COLOR_LRGB2Lab,
				   _["COLOR_LBGR2Luv"] = (int) COLOR_LBGR2Luv,
				   _["COLOR_LRGB2Luv"] = (int) COLOR_LRGB2Luv,
				   _["COLOR_Lab2LBGR"] = (int) COLOR_Lab2LBGR,
				   _["COLOR_Lab2LRGB"] = (int) COLOR_Lab2LRGB,
				   _["COLOR_Luv2LBGR"] = (int) COLOR_Luv2LBGR,
				   _["COLOR_Luv2LRGB"] = (int) COLOR_Luv2LRGB,
				   _["COLOR_BGR2YUV"] = (int) COLOR_BGR2YUV,
				   _["COLOR_RGB2YUV"] = (int) COLOR_RGB2YUV,
				   _["COLOR_YUV2BGR"] = (int) COLOR_YUV2BGR,
				   _["COLOR_YUV2RGB"] = (int) COLOR_YUV2RGB);
    // This is a crude workaround for the limit of 20 in create(). Or
    // we could just keep only the conversions to/from BGR
    return List::create(e1, e2, e3); // unlist() in R
}

// enums from https://docs.opencv.org/3.4.20/d7/d1b/group__imgproc__misc.html

Rcpp::IntegerVector cv_enum_AdaptiveThresholdTypes() {
    IntegerVector
	e = IntegerVector::create(_["ADAPTIVE_THRESH_MEAN_C"] = (int) ADAPTIVE_THRESH_MEAN_C,
				  _["ADAPTIVE_THRESH_GAUSSIAN_C"] = (int) ADAPTIVE_THRESH_GAUSSIAN_C);
    return e;
}
 
Rcpp::IntegerVector cv_enum_DistanceTransformLabelTypes() {
    IntegerVector
	e = IntegerVector::create(_["DIST_LABEL_CCOMP"] = (int) DIST_LABEL_CCOMP,
				  _["DIST_LABEL_PIXEL"] = (int) DIST_LABEL_PIXEL);
    return e;
}
 
Rcpp::IntegerVector cv_enum_DistanceTransformMasks() {
    IntegerVector
	e = IntegerVector::create(_["DIST_MASK_3"] = (int) DIST_MASK_3,
				  _["DIST_MASK_5"] = (int) DIST_MASK_5,
				  _["DIST_MASK_PRECISE"] = (int) DIST_MASK_PRECISE);
    return e;
}

Rcpp::IntegerVector cv_enum_DistanceTypes() {
    IntegerVector
	e = IntegerVector::create(_["DIST_USER"] = (int) DIST_USER,
				  _["DIST_L1"] = (int) DIST_L1,
				  _["DIST_L2"] = (int) DIST_L2,
				  _["DIST_C"] = (int) DIST_C,
				  _["DIST_L12"] = (int) DIST_L12,
				  _["DIST_FAIR"] = (int) DIST_FAIR,
				  _["DIST_WELSCH"] = (int) DIST_WELSCH,
				  _["DIST_HUBER"] = (int) DIST_HUBER);
    return e;
}
 
Rcpp::IntegerVector cv_enum_FloodFillFlags() {
    IntegerVector
	e = IntegerVector::create(_["FLOODFILL_FIXED_RANGE"] = (int) FLOODFILL_FIXED_RANGE,
				  _["FLOODFILL_MASK_ONLY"] = (int) FLOODFILL_MASK_ONLY);
    return e;
}

Rcpp::IntegerVector cv_enum_GrabCutClasses() {
    IntegerVector
	e = IntegerVector::create(_["GC_BGD"] = (int) GC_BGD,
				  _["GC_FGD"] = (int) GC_FGD,
				  _["GC_PR_BGD"] = (int) GC_PR_BGD,
				  _["GC_PR_FGD"] = (int) GC_PR_FGD);
    return e;
}
 
Rcpp::IntegerVector cv_enum_GrabCutModes() {
    IntegerVector
	e = IntegerVector::create(_["GC_INIT_WITH_RECT"] = (int) GC_INIT_WITH_RECT,
				  _["GC_INIT_WITH_MASK"] = (int) GC_INIT_WITH_MASK,
				  _["GC_EVAL"] = (int) GC_EVAL,
				  _["GC_EVAL_FREEZE_MODEL"] = (int) GC_EVAL_FREEZE_MODEL);
    return e;
}
 
Rcpp::IntegerVector cv_enum_ThresholdTypes() {
    IntegerVector
	e = IntegerVector::create(_["THRESH_BINARY"] = (int) THRESH_BINARY,
				  _["THRESH_BINARY_INV"] = (int) THRESH_BINARY_INV,
				  _["THRESH_TRUNC"] = (int) THRESH_TRUNC,
				  _["THRESH_TOZERO"] = (int) THRESH_TOZERO,
				  _["THRESH_TOZERO_INV"] = (int) THRESH_TOZERO_INV,
				  _["THRESH_MASK"] = (int) THRESH_MASK,
				  _["THRESH_OTSU"] = (int) THRESH_OTSU,
				  _["THRESH_TRIANGLE"] = (int) THRESH_TRIANGLE);
    return e;
}




// unnamed enum types

Rcpp::IntegerVector cv_enum_Misc()
{
    IntegerVector
	e = IntegerVector::create(_["RECURS_FILTER"] = (int) RECURS_FILTER,
				  _["NORMCONV_FILTER"] = (int) NORMCONV_FILTER,
				  _["INPAINT_NS"] = (int) INPAINT_NS,
				  _["INPAINT_TELEA"] = (int) INPAINT_TELEA,
				  _["NORMAL_CLONE"] = (int) NORMAL_CLONE,
				  _["MIXED_CLONE"] = (int) MIXED_CLONE,
				  _["MONOCHROME_TRANSFER"] = (int) MONOCHROME_TRANSFER);
    return e;
}



RCPP_MODULE(enums) 
{
    function("BorderTypes", &cv_enum_BorderTypes, "OpenCV BorderType enum values");
    function("DftFlags", &cv_enum_DftFlags, "OpenCV DftFlags enum values");
    function("InterpolationFlags", &cv_enum_InterpolationFlags, "OpenCV Interpolation Types enum values");
    function("ImreadModes", &cv_enum_ImreadModes, "OpenCV ImreadModes enum values");
    function("ColorConversionCodes", &cv_enum_ColorConversionCodes,
	     "OpenCV color conversion enum values");
    function("AdaptiveThresholdTypes", &cv_enum_AdaptiveThresholdTypes, "OpenCV AdaptiveThresholdTypes enum values");
    function("DistanceTransformLabelTypes", &cv_enum_DistanceTransformLabelTypes, "OpenCV DistanceTransformLabelTypes enum values");
    function("DistanceTransformMasks", &cv_enum_DistanceTransformMasks, "OpenCV DistanceTransformMasks enum values");
    function("DistanceTypes", &cv_enum_DistanceTypes, "OpenCV DistanceTypes enum values");
    function("FloodFillFlags", &cv_enum_FloodFillFlags, "OpenCV FloodFillFlags enum values");
    function("GrabCutClasses", &cv_enum_GrabCutClasses, "OpenCV GrabCutClasses enum values");
    function("GrabCutModes", &cv_enum_GrabCutModes, "OpenCV GrabCutModes enum values");
    function("ThresholdTypes", &cv_enum_ThresholdTypes, "OpenCV ThresholdTypes enum values");
    function("Misc", &cv_enum_Misc, "OpenCV miscellaneous enum values");
}
