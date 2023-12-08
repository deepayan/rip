# OpenCV functions by header file

An incomplete list of potential candidate function for which R
wrappers can be useful. This is obtained (from opencv 3.2.0) by

```
grep CV_EXPORTS_W /usr/include/opencv2/*.hpp | grep -v class
```

with further manual editing to retain only "interesting" ones. In
particular, only functions are considered here, and not classes. Note
that many interesting uses of OpenCV involve classes. For example,
video capture is done using the `VideoCapture` class in
`opencv2/videoio.hpp`, and these can be in principle exposed either
using an external pointer (which we have not tried yet), or write
custom functions (see `rip.cv$IO$vfread` for a simple example).



## /usr/include/opencv2/imgcodecs.hpp:

```
CV_EXPORTS_W Mat imread( const String& filename, int flags = IMREAD_COLOR );
CV_EXPORTS_W bool imreadmulti(const String& filename, std::vector<Mat>& mats, int flags = IMREAD_ANYCOLOR);
CV_EXPORTS_W bool imwrite( const String& filename, InputArray img,
CV_EXPORTS_W Mat imdecode( InputArray buf, int flags );
CV_EXPORTS_W bool imencode( const String& ext, InputArray img,
```

## /usr/include/opencv2/photo.hpp:

```
CV_EXPORTS_W void inpaint( InputArray src, InputArray inpaintMask,
CV_EXPORTS_W void fastNlMeansDenoising( InputArray src, OutputArray dst, float h = 3,
CV_EXPORTS_W void fastNlMeansDenoising( InputArray src, OutputArray dst,
CV_EXPORTS_W void fastNlMeansDenoisingColored( InputArray src, OutputArray dst,
CV_EXPORTS_W void fastNlMeansDenoisingMulti( InputArrayOfArrays srcImgs, OutputArray dst,
CV_EXPORTS_W void fastNlMeansDenoisingMulti( InputArrayOfArrays srcImgs, OutputArray dst,
CV_EXPORTS_W void fastNlMeansDenoisingColoredMulti( InputArrayOfArrays srcImgs, OutputArray dst,
CV_EXPORTS_W void denoise_TVL1(const std::vector<Mat>& observations,Mat& result, double lambda=1.0, int niters=30);
CV_EXPORTS_W void decolor( InputArray src, OutputArray grayscale, OutputArray color_boost);
CV_EXPORTS_W void seamlessClone( InputArray src, InputArray dst, InputArray mask, Point p,
CV_EXPORTS_W void colorChange(InputArray src, InputArray mask, OutputArray dst, float red_mul = 1.0f,
CV_EXPORTS_W void illuminationChange(InputArray src, InputArray mask, OutputArray dst,
CV_EXPORTS_W void textureFlattening(InputArray src, InputArray mask, OutputArray dst,
CV_EXPORTS_W void edgePreservingFilter(InputArray src, OutputArray dst, int flags = 1,
CV_EXPORTS_W void detailEnhance(InputArray src, OutputArray dst, float sigma_s = 10,
CV_EXPORTS_W void pencilSketch(InputArray src, OutputArray dst1, OutputArray dst2,
CV_EXPORTS_W void stylization(InputArray src, OutputArray dst, float sigma_s = 60,
```

## /usr/include/opencv2/imgproc.hpp:

In module `imgproc`:

```
CV_EXPORTS_W Mat getStructuringElement(int shape, Size ksize, Point anchor = Point(-1,-1));
CV_EXPORTS_W void erode( InputArray src, OutputArray dst, InputArray kernel,
CV_EXPORTS_W void dilate( InputArray src, OutputArray dst, InputArray kernel,
CV_EXPORTS_W void resize( InputArray src, OutputArray dst,
CV_EXPORTS_W void warpAffine( InputArray src, OutputArray dst,
CV_EXPORTS_W void warpPerspective( InputArray src, OutputArray dst,
CV_EXPORTS_W Mat getRotationMatrix2D( Point2f center, double angle, double scale );
CV_EXPORTS_W void invertAffineTransform( InputArray M, OutputArray iM );
CV_EXPORTS_W void getRectSubPix( InputArray image, Size patchSize,
CV_EXPORTS_W void pyrDown( InputArray src, OutputArray dst,
CV_EXPORTS_W void pyrUp( InputArray src, OutputArray dst,
CV_EXPORTS_W void pyrMeanShiftFiltering( InputArray src, OutputArray dst,
CV_EXPORTS_W void cvtColor( InputArray src, OutputArray dst, int code, int dstCn = 0 );
CV_EXPORTS_W void matchTemplate( InputArray image, InputArray templ,
CV_EXPORTS_W void equalizeHist( InputArray src, OutputArray dst );
```


In module `feature`:

```
CV_EXPORTS_W void Canny( InputArray image, OutputArray edges,
CV_EXPORTS_W void cornerMinEigenVal( InputArray src, OutputArray dst,
CV_EXPORTS_W void cornerHarris( InputArray src, OutputArray dst, int blockSize,
CV_EXPORTS_W void cornerEigenValsAndVecs( InputArray src, OutputArray dst,
```


In module `filter`:

```
CV_EXPORTS_W Mat getGaussianKernel( int ksize, double sigma, int ktype = CV_64F );
CV_EXPORTS_W void getDerivKernels( OutputArray kx, OutputArray ky,
CV_EXPORTS_W Mat getGaborKernel( Size ksize, double sigma, double theta, double lambd,
CV_EXPORTS_W void medianBlur( InputArray src, OutputArray dst, int ksize );
CV_EXPORTS_W void GaussianBlur( InputArray src, OutputArray dst, Size ksize,
CV_EXPORTS_W void bilateralFilter( InputArray src, OutputArray dst, int d,
CV_EXPORTS_W void boxFilter( InputArray src, OutputArray dst, int ddepth,
CV_EXPORTS_W void blur( InputArray src, OutputArray dst,
CV_EXPORTS_W void filter2D( InputArray src, OutputArray dst, int ddepth,
CV_EXPORTS_W void sepFilter2D( InputArray src, OutputArray dst, int ddepth,
CV_EXPORTS_W void Sobel( InputArray src, OutputArray dst, int ddepth,
CV_EXPORTS_W void Scharr( InputArray src, OutputArray dst, int ddepth,
```


TODO next:

```
CV_EXPORTS_W void Laplacian( InputArray src, OutputArray dst, int ddepth,
CV_EXPORTS_W void goodFeaturesToTrack( InputArray image, OutputArray corners,
CV_EXPORTS_W void HoughLines( InputArray image, OutputArray lines,
CV_EXPORTS_W void HoughLinesP( InputArray image, OutputArray lines,
CV_EXPORTS_W void HoughCircles( InputArray image, OutputArray circles,

CV_EXPORTS_W double threshold( InputArray src, OutputArray dst,
CV_EXPORTS_W void adaptiveThreshold( InputArray src, OutputArray dst,

CV_EXPORTS_W void watershed( InputArray image, InputOutputArray markers );
CV_EXPORTS_W int floodFill( InputOutputArray image, InputOutputArray mask,
```




Not available:


```
CV_EXPORTS_W Ptr<LineSegmentDetector> createLineSegmentDetector(
CV_EXPORTS_W void sqrBoxFilter( InputArray _src, OutputArray _dst, int ddepth,
CV_EXPORTS_W void spatialGradient( InputArray src, OutputArray dx,
CV_EXPORTS_W void preCornerDetect( InputArray src, OutputArray dst, int ksize,
CV_EXPORTS_W void cornerSubPix( InputArray image, InputOutputArray corners,
CV_EXPORTS_W void morphologyEx( InputArray src, OutputArray dst,
CV_EXPORTS_W void remap( InputArray src, OutputArray dst,
CV_EXPORTS_W void convertMaps( InputArray map1, InputArray map2,
CV_EXPORTS_W Mat getPerspectiveTransform( InputArray src, InputArray dst );
CV_EXPORTS_W Mat getAffineTransform( InputArray src, InputArray dst );
CV_EXPORTS_W void logPolar( InputArray src, OutputArray dst,
CV_EXPORTS_W void linearPolar( InputArray src, OutputArray dst,
CV_EXPORTS_W void integral( InputArray src, OutputArray sum, int sdepth = -1 );
CV_EXPORTS_W void accumulate( InputArray src, InputOutputArray dst,
CV_EXPORTS_W void accumulateSquare( InputArray src, InputOutputArray dst,
CV_EXPORTS_W void accumulateProduct( InputArray src1, InputArray src2,
CV_EXPORTS_W void accumulateWeighted( InputArray src, InputOutputArray dst,
CV_EXPORTS_W Point2d phaseCorrelate(InputArray src1, InputArray src2,
CV_EXPORTS_W void createHanningWindow(OutputArray dst, Size winSize, int type);
CV_EXPORTS_W void undistort( InputArray src, OutputArray dst,
CV_EXPORTS_W void initUndistortRectifyMap( InputArray cameraMatrix, InputArray distCoeffs,
CV_EXPORTS_W float initWideAngleProjMap( InputArray cameraMatrix, InputArray distCoeffs,
CV_EXPORTS_W Mat getDefaultNewCameraMatrix( InputArray cameraMatrix, Size imgsize = Size(),
CV_EXPORTS_W void undistortPoints( InputArray src, OutputArray dst,
CV_EXPORTS_W void calcHist( InputArrayOfArrays images,
CV_EXPORTS_W void calcBackProject( InputArrayOfArrays images, const std::vector<int>& channels,
CV_EXPORTS_W double compareHist( InputArray H1, InputArray H2, int method );

CV_EXPORTS_W void grabCut( InputArray img, InputOutputArray mask, Rect rect,
CV_EXPORTS_W void distanceTransform( InputArray src, OutputArray dst,
CV_EXPORTS_W void demosaicing(InputArray _src, OutputArray _dst, int code, int dcn = 0);
CV_EXPORTS_W Moments moments( InputArray array, bool binaryImage = false );
CV_EXPORTS_W void HuMoments( const Moments& m, OutputArray hu );
CV_EXPORTS_W int connectedComponents(InputArray image, OutputArray labels,
CV_EXPORTS_W int connectedComponentsWithStats(InputArray image, OutputArray labels,
CV_EXPORTS_W void findContours( InputOutputArray image, OutputArrayOfArrays contours,
CV_EXPORTS_W void approxPolyDP( InputArray curve,
CV_EXPORTS_W double arcLength( InputArray curve, bool closed );
CV_EXPORTS_W Rect boundingRect( InputArray points );
CV_EXPORTS_W double contourArea( InputArray contour, bool oriented = false );
CV_EXPORTS_W RotatedRect minAreaRect( InputArray points );
CV_EXPORTS_W void boxPoints(RotatedRect box, OutputArray points);
CV_EXPORTS_W void minEnclosingCircle( InputArray points,
CV_EXPORTS_W double minEnclosingTriangle( InputArray points, CV_OUT OutputArray triangle );
CV_EXPORTS_W double matchShapes( InputArray contour1, InputArray contour2,
CV_EXPORTS_W void convexHull( InputArray points, OutputArray hull,
CV_EXPORTS_W void convexityDefects( InputArray contour, InputArray convexhull, OutputArray convexityDefects );
CV_EXPORTS_W bool isContourConvex( InputArray contour );
CV_EXPORTS_W float intersectConvexConvex( InputArray _p1, InputArray _p2,
CV_EXPORTS_W RotatedRect fitEllipse( InputArray points );
CV_EXPORTS_W void fitLine( InputArray points, OutputArray line, int distType,
CV_EXPORTS_W double pointPolygonTest( InputArray contour, Point2f pt, bool measureDist );
CV_EXPORTS_W int rotatedRectangleIntersection( const RotatedRect& rect1, const RotatedRect& rect2, OutputArray intersectingRegion  );
CV_EXPORTS_W Ptr<CLAHE> createCLAHE(double clipLimit = 40.0, Size tileGridSize = Size(8, 8));
CV_EXPORTS_W void applyColorMap(InputArray src, OutputArray dst, int colormap);
CV_EXPORTS_W void line(InputOutputArray img, Point pt1, Point pt2, const Scalar& color,
CV_EXPORTS_W void arrowedLine(InputOutputArray img, Point pt1, Point pt2, const Scalar& color,
CV_EXPORTS_W void rectangle(InputOutputArray img, Point pt1, Point pt2,
CV_EXPORTS_W void circle(InputOutputArray img, Point center, int radius,
CV_EXPORTS_W void ellipse(InputOutputArray img, Point center, Size axes,
CV_EXPORTS_W void ellipse(InputOutputArray img, const RotatedRect& box, const Scalar& color,
CV_EXPORTS_W void drawMarker(CV_IN_OUT Mat& img, Point position, const Scalar& color,
CV_EXPORTS_W void fillConvexPoly(InputOutputArray img, InputArray points,
CV_EXPORTS_W void fillPoly(InputOutputArray img, InputArrayOfArrays pts,
CV_EXPORTS_W void polylines(InputOutputArray img, InputArrayOfArrays pts,
CV_EXPORTS_W void drawContours( InputOutputArray image, InputArrayOfArrays contours,
CV_EXPORTS_W bool clipLine(Rect imgRect, CV_OUT CV_IN_OUT Point& pt1, CV_OUT CV_IN_OUT Point& pt2);
CV_EXPORTS_W void ellipse2Poly( Point center, Size axes, int angle,
CV_EXPORTS_W void putText( InputOutputArray img, const String& text, Point org,
CV_EXPORTS_W Size getTextSize(const String& text, int fontFace,
```


## /usr/include/opencv2/aruco.hpp

```
CV_EXPORTS_W void detectMarkers(InputArray image, const Ptr<Dictionary> &dictionary, OutputArrayOfArrays corners,
CV_EXPORTS_W void estimatePoseSingleMarkers(InputArrayOfArrays corners, float markerLength,
CV_EXPORTS_W int estimatePoseBoard(InputArrayOfArrays corners, InputArray ids, const Ptr<Board> &board,
CV_EXPORTS_W void refineDetectedMarkers(
CV_EXPORTS_W void drawDetectedMarkers(InputOutputArray image, InputArrayOfArrays corners,
CV_EXPORTS_W void drawAxis(InputOutputArray image, InputArray cameraMatrix, InputArray distCoeffs,
CV_EXPORTS_W void drawMarker(const Ptr<Dictionary> &dictionary, int id, int sidePixels, OutputArray img,
CV_EXPORTS_W void drawPlanarBoard(const Ptr<Board> &board, Size outSize, OutputArray img,
CV_EXPORTS_W double calibrateCameraAruco(
```

## /usr/include/opencv2/core.hpp:


In module `transforms`:

```
CV_EXPORTS_W void dft(InputArray src, OutputArray dst, int flags = 0, int nonzeroRows = 0);
CV_EXPORTS_W void idft(InputArray src, OutputArray dst, int flags = 0, int nonzeroRows = 0);
CV_EXPORTS_W void mulSpectrums(InputArray a, InputArray b, OutputArray c,
CV_EXPORTS_W int getOptimalDFTSize(int vecsize);
```

In module `core`:

```
CV_EXPORTS_W void copyMakeBorder(InputArray src, OutputArray dst,
CV_EXPORTS_W double PSNR(InputArray src1, InputArray src2);
CV_EXPORTS_W void flip(InputArray src, OutputArray dst, int flipCode);
CV_EXPORTS_W void rotate(InputArray src, OutputArray dst, int rotateCode);
```

Not available


```
CV_EXPORTS_W int borderInterpolate(int p, int len, int borderType);
CV_EXPORTS_W void add(InputArray src1, InputArray src2, OutputArray dst,
CV_EXPORTS_W void subtract(InputArray src1, InputArray src2, OutputArray dst,
CV_EXPORTS_W void multiply(InputArray src1, InputArray src2,
CV_EXPORTS_W void divide(InputArray src1, InputArray src2, OutputArray dst,
CV_EXPORTS_W void divide(double scale, InputArray src2,
CV_EXPORTS_W void scaleAdd(InputArray src1, double alpha, InputArray src2, OutputArray dst);
CV_EXPORTS_W void addWeighted(InputArray src1, double alpha, InputArray src2,
CV_EXPORTS_W void convertScaleAbs(InputArray src, OutputArray dst,
CV_EXPORTS_W void convertFp16(InputArray src, OutputArray dst);
CV_EXPORTS_W void LUT(InputArray src, InputArray lut, OutputArray dst);
CV_EXPORTS_W int countNonZero( InputArray src );
CV_EXPORTS_W void findNonZero( InputArray src, OutputArray idx );
CV_EXPORTS_W Scalar mean(InputArray src, InputArray mask = noArray());
CV_EXPORTS_W void meanStdDev(InputArray src, OutputArray mean, OutputArray stddev,
CV_EXPORTS_W double norm(InputArray src1, int normType = NORM_L2, InputArray mask = noArray());
CV_EXPORTS_W double norm(InputArray src1, InputArray src2,
CV_EXPORTS_W void batchDistance(InputArray src1, InputArray src2,
CV_EXPORTS_W void normalize( InputArray src, InputOutputArray dst, double alpha = 1, double beta = 0,
CV_EXPORTS_W void minMaxLoc(InputArray src, CV_OUT double* minVal,
CV_EXPORTS_W void reduce(InputArray src, OutputArray dst, int dim, int rtype, int dtype = -1);
CV_EXPORTS_W void merge(InputArrayOfArrays mv, OutputArray dst);
CV_EXPORTS_W void split(InputArray m, OutputArrayOfArrays mv);
CV_EXPORTS_W void mixChannels(InputArrayOfArrays src, InputOutputArrayOfArrays dst,
CV_EXPORTS_W void extractChannel(InputArray src, OutputArray dst, int coi);
CV_EXPORTS_W void insertChannel(InputArray src, InputOutputArray dst, int coi);
CV_EXPORTS_W void repeat(InputArray src, int ny, int nx, OutputArray dst);
CV_EXPORTS_W void hconcat(InputArrayOfArrays src, OutputArray dst);
CV_EXPORTS_W void vconcat(InputArrayOfArrays src, OutputArray dst);
CV_EXPORTS_W void bitwise_and(InputArray src1, InputArray src2,
CV_EXPORTS_W void bitwise_or(InputArray src1, InputArray src2,
CV_EXPORTS_W void bitwise_xor(InputArray src1, InputArray src2,
CV_EXPORTS_W void bitwise_not(InputArray src, OutputArray dst,
CV_EXPORTS_W void absdiff(InputArray src1, InputArray src2, OutputArray dst);
CV_EXPORTS_W void inRange(InputArray src, InputArray lowerb,
CV_EXPORTS_W void compare(InputArray src1, InputArray src2, OutputArray dst, int cmpop);
CV_EXPORTS_W void min(InputArray src1, InputArray src2, OutputArray dst);
CV_EXPORTS_W void max(InputArray src1, InputArray src2, OutputArray dst);
CV_EXPORTS_W void sqrt(InputArray src, OutputArray dst);
CV_EXPORTS_W void pow(InputArray src, double power, OutputArray dst);
CV_EXPORTS_W void exp(InputArray src, OutputArray dst);
CV_EXPORTS_W void log(InputArray src, OutputArray dst);
CV_EXPORTS_W void polarToCart(InputArray magnitude, InputArray angle,
CV_EXPORTS_W void cartToPolar(InputArray x, InputArray y,
CV_EXPORTS_W void phase(InputArray x, InputArray y, OutputArray angle,
CV_EXPORTS_W void magnitude(InputArray x, InputArray y, OutputArray magnitude);
CV_EXPORTS_W bool checkRange(InputArray a, bool quiet = true, CV_OUT Point* pos = 0,
CV_EXPORTS_W void patchNaNs(InputOutputArray a, double val = 0);
CV_EXPORTS_W void gemm(InputArray src1, InputArray src2, double alpha,
CV_EXPORTS_W void mulTransposed( InputArray src, OutputArray dst, bool aTa,
CV_EXPORTS_W void transpose(InputArray src, OutputArray dst);
CV_EXPORTS_W void transform(InputArray src, OutputArray dst, InputArray m );
CV_EXPORTS_W void perspectiveTransform(InputArray src, OutputArray dst, InputArray m );
CV_EXPORTS_W void completeSymm(InputOutputArray mtx, bool lowerToUpper = false);
CV_EXPORTS_W void setIdentity(InputOutputArray mtx, const Scalar& s = Scalar(1));
CV_EXPORTS_W double determinant(InputArray mtx);
CV_EXPORTS_W Scalar trace(InputArray mtx);
CV_EXPORTS_W double invert(InputArray src, OutputArray dst, int flags = DECOMP_LU);
CV_EXPORTS_W bool solve(InputArray src1, InputArray src2,
CV_EXPORTS_W void sort(InputArray src, OutputArray dst, int flags);
CV_EXPORTS_W void sortIdx(InputArray src, OutputArray dst, int flags);
CV_EXPORTS_W int solveCubic(InputArray coeffs, OutputArray roots);
CV_EXPORTS_W double solvePoly(InputArray coeffs, OutputArray roots, int maxIters = 300);
CV_EXPORTS_W bool eigen(InputArray src, OutputArray eigenvalues,
CV_EXPORTS_W void calcCovarMatrix( InputArray samples, OutputArray covar,
CV_EXPORTS_W void PCACompute(InputArray data, InputOutputArray mean,
CV_EXPORTS_W void PCACompute(InputArray data, InputOutputArray mean,
CV_EXPORTS_W void PCAProject(InputArray data, InputArray mean,
CV_EXPORTS_W void PCABackProject(InputArray data, InputArray mean,
CV_EXPORTS_W void SVDecomp( InputArray src, OutputArray w, OutputArray u, OutputArray vt, int flags = 0 );
CV_EXPORTS_W void SVBackSubst( InputArray w, InputArray u, InputArray vt,
CV_EXPORTS_W double Mahalanobis(InputArray v1, InputArray v2, InputArray icovar);
CV_EXPORTS_W void dct(InputArray src, OutputArray dst, int flags = 0);
CV_EXPORTS_W void idct(InputArray src, OutputArray dst, int flags = 0);
CV_EXPORTS_W void setRNGSeed(int seed);
CV_EXPORTS_W void randu(InputOutputArray dst, InputArray low, InputArray high);
CV_EXPORTS_W void randn(InputOutputArray dst, InputArray mean, InputArray stddev);
CV_EXPORTS_W void randShuffle(InputOutputArray dst, double iterFactor = 1., RNG* rng = 0);
CV_EXPORTS_W double kmeans( InputArray data, int K, InputOutputArray bestLabels,
```


## /usr/include/opencv2/calib3d.hpp:

```
CV_EXPORTS_W void Rodrigues( InputArray src, OutputArray dst, OutputArray jacobian = noArray() );
CV_EXPORTS_W Mat findHomography( InputArray srcPoints, InputArray dstPoints,
CV_EXPORTS_W Vec3d RQDecomp3x3( InputArray src, OutputArray mtxR, OutputArray mtxQ,
CV_EXPORTS_W void decomposeProjectionMatrix( InputArray projMatrix, OutputArray cameraMatrix,
CV_EXPORTS_W void matMulDeriv( InputArray A, InputArray B, OutputArray dABdA, OutputArray dABdB );
CV_EXPORTS_W void composeRT( InputArray rvec1, InputArray tvec1,
CV_EXPORTS_W void projectPoints( InputArray objectPoints,
CV_EXPORTS_W bool solvePnP( InputArray objectPoints, InputArray imagePoints,
CV_EXPORTS_W bool solvePnPRansac( InputArray objectPoints, InputArray imagePoints,
CV_EXPORTS_W Mat initCameraMatrix2D( InputArrayOfArrays objectPoints,
CV_EXPORTS_W bool findChessboardCorners( InputArray image, Size patternSize, OutputArray corners,
CV_EXPORTS_W void drawChessboardCorners( InputOutputArray image, Size patternSize,
CV_EXPORTS_W bool findCirclesGrid( InputArray image, Size patternSize,
CV_EXPORTS_W double calibrateCamera( InputArrayOfArrays objectPoints,
CV_EXPORTS_W void calibrationMatrixValues( InputArray cameraMatrix, Size imageSize,
CV_EXPORTS_W double stereoCalibrate( InputArrayOfArrays objectPoints,
CV_EXPORTS_W void stereoRectify( InputArray cameraMatrix1, InputArray distCoeffs1,
CV_EXPORTS_W bool stereoRectifyUncalibrated( InputArray points1, InputArray points2,
CV_EXPORTS_W float rectify3Collinear( InputArray cameraMatrix1, InputArray distCoeffs1,
CV_EXPORTS_W Mat getOptimalNewCameraMatrix( InputArray cameraMatrix, InputArray distCoeffs,
CV_EXPORTS_W void convertPointsToHomogeneous( InputArray src, OutputArray dst );
CV_EXPORTS_W void convertPointsFromHomogeneous( InputArray src, OutputArray dst );
CV_EXPORTS_W Mat findFundamentalMat( InputArray points1, InputArray points2,
CV_EXPORTS_W Mat findEssentialMat( InputArray points1, InputArray points2,
CV_EXPORTS_W Mat findEssentialMat( InputArray points1, InputArray points2,
CV_EXPORTS_W void decomposeEssentialMat( InputArray E, OutputArray R1, OutputArray R2, OutputArray t );
CV_EXPORTS_W int recoverPose( InputArray E, InputArray points1, InputArray points2,
CV_EXPORTS_W int recoverPose( InputArray E, InputArray points1, InputArray points2,
CV_EXPORTS_W void computeCorrespondEpilines( InputArray points, int whichImage,
CV_EXPORTS_W void triangulatePoints( InputArray projMatr1, InputArray projMatr2,
CV_EXPORTS_W void correctMatches( InputArray F, InputArray points1, InputArray points2,
CV_EXPORTS_W void filterSpeckles( InputOutputArray img, double newVal,
CV_EXPORTS_W Rect getValidDisparityROI( Rect roi1, Rect roi2,
CV_EXPORTS_W void validateDisparity( InputOutputArray disparity, InputArray cost,
CV_EXPORTS_W void reprojectImageTo3D( InputArray disparity,
CV_EXPORTS_W double sampsonDistance(InputArray pt1, InputArray pt2, InputArray F);
CV_EXPORTS_W  int estimateAffine3D(InputArray src, InputArray dst,
CV_EXPORTS_W cv::Mat estimateAffine2D(InputArray from, InputArray to, OutputArray inliers = noArray(),
CV_EXPORTS_W cv::Mat estimateAffinePartial2D(InputArray from, InputArray to, OutputArray inliers = noArray(),
CV_EXPORTS_W int decomposeHomographyMat(InputArray H,
CV_EXPORTS_W void projectPoints(InputArray objectPoints, OutputArray imagePoints, InputArray rvec, InputArray tvec,
CV_EXPORTS_W void distortPoints(InputArray undistorted, OutputArray distorted, InputArray K, InputArray D, double alpha = 0);
CV_EXPORTS_W void undistortPoints(InputArray distorted, OutputArray undistorted,
CV_EXPORTS_W void initUndistortRectifyMap(InputArray K, InputArray D, InputArray R, InputArray P,
CV_EXPORTS_W void undistortImage(InputArray distorted, OutputArray undistorted,
CV_EXPORTS_W void estimateNewCameraMatrixForUndistortRectify(InputArray K, InputArray D, const Size &image_size, InputArray R,
CV_EXPORTS_W double calibrate(InputArrayOfArrays objectPoints, InputArrayOfArrays imagePoints, const Size& image_size,
CV_EXPORTS_W void stereoRectify(InputArray K1, InputArray D1, InputArray K2, InputArray D2, const Size &imageSize, InputArray R, InputArray tvec,
CV_EXPORTS_W double stereoCalibrate(InputArrayOfArrays objectPoints, InputArrayOfArrays imagePoints1, InputArrayOfArrays imagePoints2,
```

## /usr/include/opencv2/features2d.hpp:

```
CV_EXPORTS_W void drawKeypoints( InputArray image, const std::vector<KeyPoint>& keypoints, InputOutputArray outImage,
CV_EXPORTS_W void drawMatches( InputArray img1, const std::vector<KeyPoint>& keypoints1,
```

## /usr/include/opencv2/objdetect.hpp:

```
CV_EXPORTS_W void groupRectangles(CV_IN_OUT std::vector<Rect>& rectList, CV_OUT std::vector<int>& weights,
```

## /usr/include/opencv2/optflow.hpp:

```
CV_EXPORTS_W void calcOpticalFlowSF( InputArray from, InputArray to, OutputArray flow,
CV_EXPORTS_W void calcOpticalFlowSF( InputArray from, InputArray to, OutputArray flow, int layers,
CV_EXPORTS_W void calcOpticalFlowSparseToDense ( InputArray from, InputArray to, OutputArray flow,
CV_EXPORTS_W Mat readOpticalFlow( const String& path );
CV_EXPORTS_W bool writeOpticalFlow( const String& path, InputArray flow );
CV_EXPORTS_W Ptr<VariationalRefinement> createVariationalFlowRefinement();
CV_EXPORTS_W Ptr<DenseOpticalFlow> createOptFlow_DeepFlow();
CV_EXPORTS_W Ptr<DenseOpticalFlow> createOptFlow_SimpleFlow();
CV_EXPORTS_W Ptr<DenseOpticalFlow> createOptFlow_Farneback();
CV_EXPORTS_W Ptr<DenseOpticalFlow> createOptFlow_SparseToDense();
CV_EXPORTS_W Ptr<DISOpticalFlow> createOptFlow_DIS(int preset = DISOpticalFlow::PRESET_FAST);
```

## /usr/include/opencv2/plot.hpp:

```
CV_EXPORTS_W Ptr<Plot2d> createPlot2d(InputArray data);
CV_EXPORTS_W Ptr<Plot2d> createPlot2d(InputArray dataX, InputArray dataY);
```

## /usr/include/opencv2/stitching.hpp:

```
CV_EXPORTS_W Ptr<Stitcher> createStitcher(bool try_use_gpu = false);
```

## /usr/include/opencv2/ximgproc.hpp:

```
CV_EXPORTS_W void niBlackThreshold( InputArray _src, OutputArray _dst,
CV_EXPORTS_W void thinning( InputArray src, OutputArray dst, int thinningType = THINNING_ZHANGSUEN);
```

## /usr/include/opencv2/highgui.hpp:

```
CV_EXPORTS_W void namedWindow(const String& winname, int flags = WINDOW_AUTOSIZE);
CV_EXPORTS_W void destroyWindow(const String& winname);
CV_EXPORTS_W void destroyAllWindows();
CV_EXPORTS_W int startWindowThread();
CV_EXPORTS_W int waitKeyEx(int delay = 0);
CV_EXPORTS_W int waitKey(int delay = 0);
CV_EXPORTS_W void imshow(const String& winname, InputArray mat);
CV_EXPORTS_W void resizeWindow(const String& winname, int width, int height);
CV_EXPORTS_W void moveWindow(const String& winname, int x, int y);
CV_EXPORTS_W void setWindowProperty(const String& winname, int prop_id, double prop_value);
CV_EXPORTS_W void setWindowTitle(const String& winname, const String& title);
CV_EXPORTS_W double getWindowProperty(const String& winname, int prop_id);
CV_EXPORTS_W int getTrackbarPos(const String& trackbarname, const String& winname);
CV_EXPORTS_W void setTrackbarPos(const String& trackbarname, const String& winname, int pos);
CV_EXPORTS_W void setTrackbarMax(const String& trackbarname, const String& winname, int maxval);
CV_EXPORTS_W void setTrackbarMin(const String& trackbarname, const String& winname, int minval);
CV_EXPORTS_W void addText(const Mat& img, const String& text, Point org, const String& nameFont, int pointSize = -1, Scalar color = Scalar::all(0),
CV_EXPORTS_W void displayOverlay(const String& winname, const String& text, int delayms = 0);
CV_EXPORTS_W void displayStatusBar(const String& winname, const String& text, int delayms = 0);
```
