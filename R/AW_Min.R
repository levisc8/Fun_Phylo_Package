#' Calculate abundance weighted NND
#' 
#' @description Calculates abundance weighted minimum
#' 
#' @param x a vector or matrix of numerical values
#' @param w a vector or matrix of numerical weights whose positions
#' correspond to the values in \code{x}.
#' @param na.rm A logical indicating whether to ignore missing values.
#' 
#' @return A numerical value
#' 
#' @author Sam Levin
#' 
#' @seealso \code{\link[stats]{weighted.mean}}
#' 
#' @export


weighted.min <- function (x, w, na.rm = FALSE) 
{
  if (missing(w)) {
    if (na.rm) 
      x <- x[!is.na(x)]
    return(sum(x)/length(x))
  }
  if (length(w) != length(x)) 
    stop("'x' and 'w' must have the same length")
  w <- as.double(w)
  if (na.rm) {
    i <- !is.na(x)
    w <- w[i]
    x <- x[i]
  }
  min(((x * w)[w != 0])/sum(w))
}
