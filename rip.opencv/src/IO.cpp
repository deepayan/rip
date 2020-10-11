#include "common.h"

// Rcpp::NumericMatrix CV2RCPP(cv::Mat);
// cv::Mat RCPP2CV(Rcpp::NumericMatrix, int);

using namespace cv;
using namespace Rcpp;

// flags: opencv 2:
// cv::IMREAD_UNCHANGED = -1, 
// cv::IMREAD_GRAYSCALE = 0, 
// cv::IMREAD_COLOR = 1, 
// ...
// (opencv3 has more flags).

Rcpp::NumericMatrix cv_imread(std::vector<std::string> file, int flags = 0)
{
    std::string f = file[0];
    cv::Mat image = cv::imread(f, flags);
    if (!image.data) stop("cv::imread() could not read image data from file.");
    return CV2RCPP(image);
}

// cv::imwrite() 
void cv_imwrite(Rcpp::NumericMatrix x, Rcpp::CharacterVector file) 
{
    cv::Mat img = RCPP2CV(x, 5);
    std::string outfile = as<std::string>(file);
    cv::imwrite(outfile, img);
    return;
}


// Read frame from video file, but very crudely

Rcpp::NumericMatrix readVideoFrame(std::string file, int index)
{
    cv::VideoCapture cap(file);
    int i = 0;
    if (!cap.isOpened())
        stop("Could not open the file");
    else {
        for( ; ; ) {
	    cv::Mat frame;
            cap >> frame;
            i++;
            if (i == index)
                return (CV2RCPP(frame));
        }
    }
}


RCPP_MODULE(IO) {
    function("imread", &cv_imread,
	     List::create(_["file"], _["type"] = -1L),
	     "Read image from file: type 0=grayscale, 1=color");
    function("imwrite", &cv_imwrite, "Write image to file");
    function("vfread", &readVideoFrame, "Read frame from video file");
}


