#' Visualize the Mesh Panel Zones of a Midwater Trawl
#'
#' Draw a diagram displaying the sizes of the mesh panel zones of a midwater trawl.
#'
#' @param PanelProps 	A numeric vector of length 4, size of the different mesh panel zones of the midwater trawl, 
#' mouth (outermost), middle, aft, and cod (inner), default c(0.4, 0.3, 0.2, 0.1).
#' Sizes are expressed as proportions of the distance from the outer edge of the trawl to the trawl center in both the vertical and horizontal directions,
#' and they should add up to 1.
#'
#' @details
#'
#' A diagram is produced giving a view of the midwater trawl mesh panel zones (drawn to scale) from the perspective a fish 
#' located directly in the center of the oncoming trawl path.
#'
#' @export
#' @import 				MASS
#' @seealso \code{\link{AcMtEst}}
#' @examples
#'
#' \dontrun{
#' ViewZones()
#' ViewZones(c(0.4, 0.4, 0.1, 0.1))
#' }
#'

ViewZones <- function(PanelProps=c(0.4, 0.3, 0.2, 0.1)) {
	# check validity of the trawl zone proportions that were input
	if(round(sum(PanelProps), 7) != 1) stop("Panel proportions should sum to 1.", fisherr)

	mouth.edge <- 1
	middle.edge <- 1 - PanelProps[1]
	aft.edge <- middle.edge - PanelProps[2]
	cod.edge <- aft.edge - PanelProps[3]

	dev.new(rescale="fit")
	par(mar=rep(0.1, 4), cex=2)
	eqscplot(0, 0, xlim=c(-1, 1), ylim=c(-1, 1), type="n", axes=FALSE, xlab="", ylab="")
	polygon(mouth.edge*c(-1, 1, 1, -1), mouth.edge*c(-1, -1, 1, 1), col=gray(0.4))
	polygon(middle.edge*c(-1, 1, 1, -1), middle.edge*c(-1, -1, 1, 1), col=gray(0.6))
	polygon(aft.edge*c(-1, 1, 1, -1), aft.edge*c(-1, -1, 1, 1), col=gray(0.8))
	polygon(cod.edge*c(-1, 1, 1, -1), cod.edge*c(-1, -1, 1, 1), col=gray(1))
	text(-middle.edge, middle.edge, "Mouth", pos=3, offset=0.1, col=gray(1))
	text(-aft.edge, aft.edge, "Middle", pos=3, offset=0.1, col=gray(0))
	text(-cod.edge, cod.edge, "Aft", pos=3, offset=0.1, col=gray(0.2))
	text(0, 0, "Cod", col=gray(0.4))
	}
