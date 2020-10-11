
#include "common.h"

using namespace cv;
using namespace Rcpp;


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
	e = IntegerVector::create(_["DFT_INVERSE"] = (int) DFT_INVERSE,
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
    function("Misc", &cv_enum_Misc, "OpenCV miscellaneous enum values");
}
