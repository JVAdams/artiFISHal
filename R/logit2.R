#' Double Logistic Probability
#'
#' Probability distribution function of double logistic curve, with ascending
#' and descending limbs allowing for parabolic-shaped probabilities.
#'
#' @param x
#'   A numeric vector.
#' @param x50a
#'   A numeric scalar, the quantile at which small values of x have a
#'   50\% probability, the first "location" parameter.
#' @param slopea
#'   A numeric scalar, the slope at which the probability for small values of
#'   x increases, the first "scale" parameter.
#' @param x50b
#'   A numeric scalar, the quantile at which large values of x have a
#'   50\% probability, the second "location" parameter.
#' @param slopeb
#'   A numeric scalar, the slope at which the probability for large values of
#'   x decreases, the second "scale" parameter.
#'
#' @details
#'   The double logistic function mimics the single logistic function if either
#'   of the location parameters are set to the extremes, \code{x50a = -Inf} or
#'   \code{x50a = Inf}.
#'
#' @export
#' @references
#'   This function is based on a modification of MATLAB
#'   \href{http://www.mathworks.com/products/matlab/}{[link]} code provided by
#'   Kresimir Williams, NOAA-AFSC
#'   \href{http://www.afsc.noaa.gov}{afsc.noaa.gov}, at the Great Lakes Acoustic
#'   Users Group's Workshop on Trawl Performance,
#'   hosted by the Great Lakes Fishery Commission
#'   \href{http://www.glfc.org/}{glfc.org}, 22-24 April 2014, in Ann Arbor,
#'   Michigan, USA.
#' @examples
#'
#' x <- 1:400
#' y1 <- logit2(x=x, x50a=90, slopea=10, x50b=Inf, slopeb=-20)
#' plot(x, y1)
#'
#' y2 <- logit2(x=x, x50a=-Inf, slopea=10, x50b=300, slopeb=-20)
#' plot(x, y2)
#'
#' y3 <- logit2(x=x, x50a=90, slopea=10, x50b=300, slopeb=-20)
#' plot(x, y3)
#'

logit2 <- function(x, x50a, slopea, x50b, slopeb) {
	proba <- plogis(x, location=x50a, scale=slopea)
	probb <- 1 - plogis(x, location=x50b, scale=-slopeb)
	proba*probb
}
