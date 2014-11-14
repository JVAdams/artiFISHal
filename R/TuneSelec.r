#' Tune Midwater Trawl Selectivity
#'
#' Interactive visualization of midwater trawl selectivity curve for aid in tuning selectivity parameters used in \code{\link{AcMtEst}}.
#'
#' @details
#'
#' Interactive sliders are provided for the four parameters of interest which are required for defining the midwater trawl selectivity
#' of the mesh panel zones as a double logistic curve:
#' \itemize{
#'   \item \code{MtL50Small} the length (in mm) at which small fish have a 50\% probability of being captured by the trawl
#'   \item \code{MtSlopeSmall} the slope at which small fish probability of capture increases with length
#'   \item \code{MtL50Large} the length (in mm) at which large fish have a 50\% probability of being captured by the trawl
#'   \item \code{MtSlopeLarge} the (absolute value of the) slope at which large fish probability of capture decreases with length
#'	}
#'
#' Note that the sliders pop up in a separate R window, which may be hidden if you click on another window.
#'
#' @export
#' @import 				rpanel
#' @seealso \code{\link{AcMtEst}}, \code{\link{logit2}}
#' @examples
#'
#' \dontrun{
#' TuneSelec()
#' }
#'

TuneSelec <- function() {

	# no idea why, but the function doesn't seem to work unless this package is attached
	library(rpanel)

	# probability graphing function
	double.draw <- function(panel) {
		y <- logit2(panel$x, panel$L501, panel$SR1, panel$L502, -panel$SR2)
		plotblank(xlim=range(panel$x[y>0.001]), ylim=0:1, xlab="Fish length  (mm)", ylab="Midwater trawl selectivity")
		abline(v=c(panel$L501, panel$L502), col="gray", lwd=2)
		abline(h=c(0, 0.5, 1), col="gray", lwd=2)
		lines(spline(panel$x, y, 1000), lwd=3)
		panel
		}

	# plot it, with a slider to adjust coeficients of the double logistic function
	dev.new()
	par(mar=c(4, 4, 1, 1))
	plot(1, 1)
	panel <- rp.control(x=1:2000, L501=10.1, SR1=10.1, L502=200.1, SR2=20.1)
	# draw an initial plot, so user isn't staring at a gray window before clicking on the slider
	double.draw(panel)
	rp.slider(panel, L501, 0.1, 1000, resolution=0.1, showvalue=T, action=double.draw, 
		title="MtL50Small:   length at 50% for small fish                                                                                                                 ")
	rp.slider(panel, SR1, 0.1, 200, resolution=0.1, showvalue=T, action=double.draw, title="MtSlopeSmall:   slope for small fish")
	rp.slider(panel, L502, 0.1, 1000, resolution=0.1, showvalue=T, action=double.draw, title="MtL50Large:   length at 50% for large fish")
	rp.slider(panel, SR2, 0.1, 200, resolution=0.1, showvalue=T, action=double.draw, title="MtSlopeLarge:   slope for large fish")

	}
