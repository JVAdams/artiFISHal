#' Predict Node from a Fitted rpart Object
#'
#' Predict the node identification number from a fitted \code{rpart} object.
#'
#' @param object
#'   A fitted model object of class "rpart", the result of a call to the
#'   \code{\link[rpart]{rpart}} function.
#' @param newdata
#'   A data frame containing the values for which predictions will be made.
#'   The predictors referred to in the right side of \code{formula(object)}
#'   must be present by name in newdata.
#' @export
#' @import
#'   rpart
#' @seealso
#'   \code{\link[rpart]{rpart}}, \code{\link[rpart]{predict.rpart}}
#' @references
#'   This function is a modification of an approach suggested in a post
#'   \href{http://tolstoy.newcastle.edu.au/R/e4/help/08/07/17702.html}{[link]}
#'   to the R-help mailing list on 22 July 2008 by Brian D. Ripley
#'   \href{http://www.stats.ox.ac.uk/~ripley/}{[link]},
#'   Professor of Applied Statistics, University of Oxford, Oxford, UK.
#' @examples
#' \dontrun{
#'
#' z.auto <- rpart(Mileage ~ Weight, car.test.frame)
#' prednode(z.auto, car.test.frame)
#'
#' fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
#' prednode(fit, kyphosis)
#' }

prednode <- function(object, newdata) {
  library(rpart)
  fax <- sapply(newdata, class)=="factor"
  newdata[, fax] <- data.frame(lapply(newdata[, fax, drop=FALSE], as.numeric))
  newmat <- rpart:::rpart.matrix(newdata)
  rpart:::pred.rpart(object, newmat)
}
