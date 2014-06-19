#' Determine Easting from Bottom Depth
#'
#' Determine the easting from a given distance to lake bottom, incorporating random assignment to west or east shore if not specified.
#' @param z 		A numeric vector of the distances to lake bottom to convert (in m).
#' @param maxz		A numeric scalar, the maximum bottom depth of the lake (in m).
#' @param ints 		A numeric vector of length 2, the intercepts of the angled lake beds along the west and east shores of the lake (in m).
#' @param slopes 	A numeric vector of length 2, the slopes of the angled lake beds along the west and east shores of the lake (unitless).
#' @param shore 	A numeric or character scalar, indicating if the shore should be west (0), east (1) or randomly assigned ("random", default).
#' @return 			A numeric vector, same length as \code{x}, with distances to bottom for each easting provided (in m).
#' @export
#' @details			An internal function called by \code{\link{SimFish}}.
#' @seealso			\code{\link{SimFish}}, \code{\link{dfromx}}, \code{\link{zfromx}}.
#' The assignment to shore is necessary because the artificial lake will have two eastings for each bottom depth (except for the depth at the vertex).

xfromz <- function(z, maxz, ints, slopes, shore="random") {
	x <- rep(NA, length(z))
	side <- if(shore[1]=="random") sample(0:1, length(z), replace=TRUE) else shore
	x[!is.na(z) & side==0] <- (z[!is.na(z) & side==0] - ints[1]) / slopes[1]
	x[!is.na(z) & side==1] <- (z[!is.na(z) & side==1] - ints[2]) / slopes[2]
	x[z >= maxz] <- runif(sum(z >= maxz), (maxz - ints[1]) / slopes[1], (maxz - ints[2]) / slopes[2])
	x
	}
