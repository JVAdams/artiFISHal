#' View Midwater Trawl Selectivity Curves
#'
#' Visualization of midwater trawl selectivity curves based on parameters used
#' in \code{\link{AcMtEst}}.
#'
#' @param SelecParam
#'   A data frame with 6 columns in which each row provides the midwater trawl
#'   selectivity parameters for a given fish group and mesh panel zone.
#'   All columns must be completely filled in (no missing values).
#'   Selectivity is assumed to be 100\% for any group-zone combination not
#'   represented as a row in the data frame.
#'   For 100\% selectivity of small fish, use MtL50Small = -Inf and any slope.
#'   For 100\% selectivity of large fish, use MtL50Large = Inf and any slope.
#'   Column names and descriptions:
#'   \itemize{
#'     \item \code{G} = character, a one-letter nickname for the group
#'       (e.g., fish species and life stage) used in plotting
#'     \item \code{Zone} = character, mesh panel zone, one of "mouth",
#'       "middle", "aft", or "cod"
#'     \item \code{MtL50Small} = the length (in mm) at which small fish have a
#'       50\% probability of being captured by the trawl
#'     \item \code{MtSlopeSmall} = the (inverse) slope at which small fish
#'       probability of capture increases with length,
#'       smaller values are steeper
#'     \item \code{MtL50Large} = the length (in mm) at which large fish have a
#'       50\% probability of being captured by the trawl
#'     \item \code{MtSlopeLarge} = the (absolute value of the inverse) slope at
#'       which large fish probability of capture decreases with length,
#'       smaller values are steeper
#'   }
#'
#' @export
#' @import
#'   jvamisc
#' @seealso
#'   \code{\link{AcMtEst}}, \code{\link{logit2}}
#' @examples
#' \dontrun{
#'
#' selec <- data.frame(
#' 	G = c("A", "a", "A", "a", "A", "a"),
#' 	Zone = c("mouth", "mouth", "middle", "middle", "aft", "aft"),
#' 	MtL50Small = c(100, 90, 60, 50, 30, 2),
#' 	MtSlopeSmall = c(40, 40, 30, 30, 20, 20),
#' 	MtL50Large = c(180, 180, Inf, Inf, Inf, Inf),
#' 	MtSlopeLarge = c(20, 20, 100, 100, 100, 100)
#' )
#' ViewSelec(selec)
#' }

ViewSelec <- function(SelecParam) {

	# check validity of zones
	suz <- c("mouth", "middle", "aft", "cod")
	uz <- unique(SelecParam$Zone)
	badzones <- setdiff(uz, suz)
	if (length(badzones) > 0) {
    stop('Zones must be one of "mouth", "middle", "aft", or "cod".')
	}

	# check for missings
	missings <- sum(is.na(SelecParam))
	if (missings > 0) {
    stop("SelectParam data frame may not have any missing values.")
	}

	# calculate lengths at small probabilities for plotting limits
	smallp <- 0.01
	SelecParam$left <- SelecParam$MtL50Small -
    log((1-smallp)/smallp)*SelecParam$MtSlopeSmall
	SelecParam$left[SelecParam$left < 0] <- 0
	SelecParam$right <- SelecParam$MtL50Large -
    log((1-smallp)/smallp)*(-SelecParam$MtSlopeLarge)

	# fill in 100% selectivities for group-zones with no parameters
	sug <- sort(unique(SelecParam$G))
	full <- expand.grid(G=sug, Zone=suz)
	both <- merge(SelecParam, full, all=TRUE)
	sel100 <- is.na(both$MtL50Small)
	both$MtL50Small[sel100] <- -Inf
	both$MtSlopeSmall[sel100] <- 100
	both$MtL50Large[sel100] <- Inf
	both$MtSlopeLarge[sel100] <- 200

	x <- floor(min(both$left[is.finite(both$left)], na.rm=TRUE)):
    ceiling(max(both$right[is.finite(both$right)], na.rm=TRUE))

	dev.new()
	par(mfrow=c(2, 2), mar=c(3, 3, 2, 1), oma=c(2, 2, 0, 0))
	for(z in 1:length(suz)) {
		sel <- both$Zone==suz[z]
		plotblank(xlim=range(x), main=suz[z])
		abline(h=c(0, 0.5, 1), col="gray", lwd=2)
		for(g in 1:length(sug)) {
			sel2 <- both$Zone==suz[z] & both$G==sug[g]
			y <- logit2(x, both$MtL50Small[sel2], both$MtSlopeSmall[sel2],
        both$MtL50Large[sel2], -both$MtSlopeLarge[sel2])
			lines(spline(x, y, 1000), lwd=3, col=g)
			}
		}
	legend("bottomright", sug, col=seq(sug), lwd=3, bty="n")
	mtext("Fish length  (mm)", side=1, outer=TRUE)
	mtext("Midwater trawl selectivity", side=2, outer=TRUE)

}
