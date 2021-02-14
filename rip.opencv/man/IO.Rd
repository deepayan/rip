\name{rip.import}
\alias{rip.import}
\alias{rip.export}
\title{Import from or export to an image file}
\usage{
rip.import(file, type = c("grayscale", "color", "original"))
rip.export(x, file)
}
\arguments{
  \item{file}{Path of the image file.}

  \item{x}{An object of class \code{"rip"} containing data that can be
    interpreted as an image by OpenCV (usually with 1, 3, or 4 channels
    with values between 0 and 255).
  }

  \item{type}{Character string that determines whether the imported data
    will be stored as grayscale (default) or color (including possibly
    an alpha channel). If \code{"original"}, the choice is taken from
    the image.
  }
}
\value{

  If successful, \code{rip.import} returns the image data in \code{file}
  as an object of \code{\link{class}} \code{"rip"}. Such objects store
  data as a (typically numeric) matrix, regardless of the number of
  channels, in the same layout that is used natively by the
  \code{cv::Mat} type in OpenCV. The \code{"cvdim"} attribute gives the
  number or rows, columns, and channels as an integer vector of length 3
  with names \code{"nrow"}, \code{"ncol"}, and \code{"nchannel"}.

  Essentially, \code{"rip"} objects store data column-wise, but all
  channels for a given column in the image are stored in contiguous
  columns. For the purposes of importing and exporting color images, the
  first three channels are interpreted as Blue, Green, and Red (BGR)
  respectively (which is not the more typical Red, Green, and Blue, RGB
  ordering).
  
}
\description{
  Import image data from an external file or write image data to an
  external file.
}
\seealso{
  The \code{"rip"} class is a convenient format for manipulating data
  using OpenCV functions as it closely reflects the native OpenCV
  representation. It is also convenient for manipulation as a matrix
  using R functions as long it represents a single channel image. For
  multi-channel images, it may be more convenient to convert it either
  into a three-dimensional array or a
  \code{\link[graphics:rasterImage]{raster}} object using convenient
  \code{\link[as.rip]{conversion functions}}.
}
\details{
  User-friendly interface to the OpenCV \code{cv::imread} and
  \code{cv::imwrite} functions via the \code{rip.cv$IO} module. The file
  format is automatically determined by OpenCV.
}
\author{
  Kaustav Nandy
} 
