#' Determine the Largest Fish that can Pass through a Net
#'
#' Determine the largest height (or depth) of a fish that can pass through a single diamond-shaped mesh of a net.
#' @param BarMesh 		A numeric scalar, the length of one side of a square mesh (in inches).
#' @param H2WRatio		A numeric scalar, the ratio of the height (vertex to vertex) of a single square mesh (oriented as a diamond) 
#' to its width (vertex to vertex).
#' @param L2HRatio		A numeric scalar, > 1, the ratio of the length of a fish to its height (or depth).  If set to NA, the default, 
#' no length calculations will be done.
#' @param Plot			A logical scalar, indicating if a diagram (drawn to scale) of the mesh and the largest fish should be shown, default TRUE.
#'
#' @return 				A named vector with 2 elements: fish height (in inches) and fish length (in mm).
#'
#' @details The cross sectional shape of the fish (looking it right in the nose) is assumed to be a geometric shape similar to an ellipse,
#' two vertically-oriented tangential circles of equal size (like a figure 8).
#' 
#' @import 				MASS plotrix
#' @export
#' @references 
#' This function is based on a modification of ellipse-based \href{http://www.mathworks.com/products/matlab/}{MATLAB} code provided by Kresimir Williams, 
#' \href{http://www.afsc.noaa.gov}{NOAA-AFSC}, at the Great Lakes Acoustic Users Group's Workshop on Trawl Performance, 
#' hosted by the \href{http://www.glfc.org/}{Great Lakes Fishery Commission}, 22-24 April 2014, in Ann Arbor, Michigan, USA.
#'
#' @examples
#' MeshPass(BarMesh=2, H2WRatio=0.3)
#' MeshPass(BarMesh=2, H2WRatio=0.3, L2HRatio=4)
#' MeshPass(BarMesh=2, H2WRatio=1/0.3, L2HRatio=4)
#' 

MeshPass <- function(BarMesh, H2WRatio, L2HRatio=NA, Plot=TRUE) {
	if(length(BarMesh)!=1 | length(H2WRatio)!=1 | length(L2HRatio)!=1 | length(Plot)!=1) stop("All arguments to the MeshPass() function should be scalars of length 1.")
	hwr <- 1/H2WRatio
	# half width of partially open mesh
	W <- BarMesh / sqrt(1 + hwr^2)
	# half height of partially open mesh
	H <- W * hwr
	# half angle of mesh peak (should be > 90 degrees, or pi/2)
	theta <- atan(1/hwr)
	# radius of biggest circle that can fit inside right (or upper) half of partially open mesh
	r <- H * tan(theta/2)
	# fish body depth in inches
	fheight <- 4*r
	# fish length in mm
	flength <- fheight*L2HRatio*25.4

	if(Plot) {
		addfish <- function(radius) {
			# draw a "fish" on the mesh ... consisting of two tangent circles, with an overlaying square
			polygon(c(-radius, radius, radius, -radius), c(-radius, -radius, radius, radius), col="lightgray", border=NA)
			draw.ellipse(x=0, y= radius, a=radius, b=radius, angle=0, col="lightgray", border="darkgray")
			draw.ellipse(x=0, y=-radius, a=radius, b=radius, angle=0, col="lightgray", border="darkgray")
			arrows(0, -2*radius, 0, 2*radius, length=0.1, code=3, col=gray(0.3), lwd=2)
			}
		if(is.na(L2HRatio)) {
			title1 <- paste0('Fish height = ', signif(fheight, 3), '"')
			title2 <- paste0('Bar mesh = ', BarMesh, '", mesh ratio = ', signif(H2WRatio, 3))
			} else {
			title1 <- paste0('Fish height = ', signif(fheight, 3), '", length = ', signif(flength, 3), ' mm')
			title2 <- paste0('Bar mesh = ', BarMesh, '", mesh ratio = ', signif(H2WRatio, 3), ', fish length/height = ', signif(L2HRatio, 3))
			}
		windows(rescale="fit")
		eqscplot(0, 0, type="n", xlim=c(-H, H), ylim=c(-W, W), las=1, xlab="Inches", ylab="Inches", main=title1)
		addfish(r)
		polygon(c(0, H, 0, -H), c(-W, 0, W, 0), lwd=2, density=0)
		mtext(title2, side=3, cex=1.2, line=0.2)
		}

	c(FishHeight.in=fheight, FishLength.mm=flength)
	}
