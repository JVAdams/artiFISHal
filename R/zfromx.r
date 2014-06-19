#' Determine Bottom Depth from Easting
#'
#' Determine the distance to the bottom from a given easting.
#' @param x 		A numeric vector of the eastings to convert (in m).
#' @param maxz		A numeric scalar, the maximum bottom depth of the lake (in m).
#' @param eastr 	A numeric vector of length 2, easting range (minimum and maximum) corresponding to the west and east shores of the lake (in m).
#' @param ints 		A numeric vector of length 2, the intercepts of the angled lake beds along the west and east shores of the lake (in m).
#' @param slopes 	A numeric vector of length 2, the slopes of the angled lake beds along the west and east shores of the lake (unitless).
#' @return 			A numeric vector, same length as \code{x}, with distances to bottom for each easting provided (in m).
#' @export
#' @details			An internal function called by \code{\link{SimFish}}.
#' @seealso			\code{\link{SimFish}}, \code{\link{dfromx}}, \code{\link{xfromz}}.

zfromx <- function(x, maxz, eastr, ints, slopes) {
	z <- rep(NA, length(x))
	z[!is.na(x) & x <= eastr[2]/3] <- slopes[1]*x[!is.na(x) & x <= eastr[2]/3] + ints[1]
	z[!is.na(x) & x >  eastr[2]/3] <- slopes[2]*x[!is.na(x) & x >  eastr[2]/3] + ints[2]
	ifelse(z > maxz, maxz, z)
	}
