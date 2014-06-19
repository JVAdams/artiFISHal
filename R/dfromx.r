#' Determine Distance to Shore from Easting
#'
#' Determine the distance to the lake shore from a given easting.
#' @param x 		A numeric vector of the eastings to convert (in m).
#' @param d2shr.we	A numeric vector of length 2, distance (in the "x" direction) from west and east shores excluded from lake (in m).
#' @param eastr 	A numeric vector of length 2, easting range (minimum and maximum) corresponding to the west and east shores of the lake (in m).
#' @return 			A numeric vector, same length as \code{x}, with distances to shore (west or east, whichever is closer) for each easting provided (in m).
#' @export
#' @details			An internal function called by \code{\link{SimFish}}.
#' @seealso			\code{\link{SimFish}}, \code{\link{zfromx}}, \code{\link{xfromz}}.

dfromx <- function(x, d2shr.we, eastr) {
	mid.d <- ((eastr[1]-d2shr.we[1]) + (eastr[2]+d2shr.we[2]) ) / 2
	d <- rep(NA, length(x))
	d[!is.na(x) & x <= mid.d] <-  d2shr.we[1] + x[!is.na(x) & x <= mid.d]
	d[!is.na(x) & x >  mid.d] <- eastr[2] + d2shr.we[2] - x[!is.na(x) & x > mid.d]
	d
	}
