\name{zapsmallp}
\alias{zapsmallp}
\title{Set small values to zero}

\description{
  Set small values in a numeric vector to zero, similar to
  \code{\link{zapsmall}}, except that non-zero elements are retained
  as is instead of being rounded.
}

\usage{
zapsmallp(x, digits = 2, prop = 10^(-digits), threshold = max(abs(x)) * prop)
}
\arguments{
  
  \item{x}{Usually a numeric vector, matrix, or array.}

  \item{digits}{Numeric giving "number" of decimal digits to retain}

  \item{prop}{Proportion (fraction) of maximum below which elements are
    to be set to zero. Overrides \code{digits}.}

  \item{threshold}{Threshold absolute value below which elements
    are to be set to zero. Overrides \code{prop}.}
}

\value{
  Modified version of \code{x}, retaining attributes, with small
  elements set to zero.
}
