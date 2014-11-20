#' Simulate a Fish Population
#'
#' Create a simulated population of pelagic fish in an artificial lake.
#' @param LakeName 		A character scalar, full name for artificial lake to be used in plot titles.
#' @param LkWidth		A numeric scalar, the width of the lake in the west-east direction (in m).
#' @param LkLength		A numeric scalar, the length of the lake in the south-north direction (in m).
#' @param BotDepMin		A numeric scalar, the minimum bottom depth of the lake, at both the west and east shorelines (in m).
#' @param BotDepMax		A numeric scalar, the maximum bottom depth of the lake (in m).
#' @param BotDepVertex	A numeric scalar, the vertical distance from the surface to the "vertex" of the lake bottom (in m), default \code{2*BotDepMax}.  
#' The "vertex" of the lake bottom is the point at which the angled lake beds along the west and east shores would intersect,
#' were they not cut off first by the specified \code{BotDepMax}.  
#' View this figure \href{https://raw.githubusercontent.com/JVAdams/artiFISHal/master/images/LakeFigures.JPG}{[link]} for a diagram of the artificial lake.
#' @param FishParam		A data frame with 18 columns in which each row describes a sub-population of fish to be placed in the artificial lake.
#' The first 11 columns must be completely filled in (no missing values).
#' The last 8 columns may have some missing values.  
#' However, in each row, \strong{either} water depth (\code{WD} and \code{WDE}) \strong{or} distance to bottom (\code{D2B} and \code{D2BE}) 
#' \strong{must} be filled in, but \strong{not both}.
#' Column names and descriptions:
#' \itemize{
#'   \item \code{G} = character, a one-letter nickname for the group (e.g., fish species and life stage) used in plotting
#'   \item \code{Z} = numeric, mean length (in mm)
#'   \item \code{ZE} = numeric, error around mean length, expressed as SD/mean
#'   \item \code{LWC1}, \code{LWC2} = numeric, length-weight regression coefficients, where wt = LWC1*len^LWC2, (wt in g, len in mm)
#'   \item \code{LWCE} = numeric, error around weight estimate, expressed as SD(estimate)/estimate
#'   \item \code{TSC1}, \code{TSC2} = numeric, target strength and length relation coefficients, ts = TSC1 + TSC2*log10(len/10), (ts in db, len in mm)
#'   \item \code{TSCE} = numeric, error around target strength estimate, expressed as SD(estimate)/estimate
#'   \item \code{PropN} = numeric, approximate proportion of population that the row represents (automatically adjusted to ensure they sum to 1)
#'   \item \code{E} = numeric, mean easting (m)
#'   \item \code{EE} = numeric, error around easting, expressed as SD/mean
#'   \item \code{N} = numeric, mean northing (m)
#'   \item \code{NE} = numeric, error around northing, expressed as SD/mean
#'   \item \code{WD} = numeric, mean water depth (m)
#'   \item \code{WDE} = numeric, error around water depth, expressed as SD/mean
#'   \item \code{D2B} = numeric, mean distance to bottom (m)
#'   \item \code{D2BE} = numeric, error around distance to bottom, expressed as SD/mean
#' }
#' @param TotNFish		A numeric scalar indicating the target number of fish to put in the lake.  
#' The actual number of fish in the population will likely be smaller than \code{TotNFish}, because the process used to populate the lake with fish
#' ends up with some fish out of water (beyond the boundaries of the artificial lake), which are then removed from the population.
#' Memory on your computer limits the size of \code{TotNFish} (see Details).
#' @param TSRange		A numeric vector of length 2, the range of target strengths to use for the fish (in db), default c(-65, -20).
#' @param PlotsPdf		A character scalar, name of pdf file to store the diagnostic plots in.  If NA, the default, 
#' plots are displayed on the screen instead.  If FALSE, no plots are created.
#' @param Seed			An integer scalar, starting seed for stochasticity incorporated in fish population generation.  
#' Use \code{Seed} to ensure the same population is generated with each call to \code{SimFish}.  
#' Otherwise, if set to NULL, the default, a random seed is used, resulting in a different population with each call to \code{SimFish}.  
#'
#' @return 				A list with 6 elements:
#' \itemize{
#' 	\item \code{Truth}, a data frame with the total number and weight 
#' of each fish group in the population; 
#' 	\item \code{LakeInfo}, a list with the lake inputs supplied as arguments to \code{SimFish} as well as a few additional objects 
#' which are used by \code{\link{SampFish}} in surveying the population,
#' 		\itemize{
#'   		\item \code{ints} = a numeric vector of length 2, the intercepts of the angled lake beds along the west and east shores of the lake (in m)
#'   		\item \code{slopes} = a numeric vector of length 2, the slopes of the angled lake beds along the west and east shores of the lake (unitless)
#'   		\item \code{d2shr.we} = a numeric vector of length 2, distance (in the west-east direction) from west and east shores excluded from lake (in m);
#' 		}
#' 	\item \code{FishInfo}, a list with the fish inputs supplied as arguments to \code{SimFish};
#' 	\item \code{FishParam}, the data frame supplied as an argument to \code{SimFish};
#' 	\item \code{FishPop}, a data frame in which each row is a fish, and 10 columns describe the fish group (\code{G}), location 
#' (easting \code{f.east}, northing \code{f.north}, distance to shore \code{f.d2sh}, bottom depth \code{f.botdep}, 
#' water depth \code{f.wdep}, distance to bottom \code{f.d2bot}, all in m), and 
#' fish size (total length in mm \code{len}, weight in g \code{wt}, and target strength in db \code{ts}); and
#' 	\item \code{PropExcluded}, a numeric vector showing the proportion of the requested number of fish, \code{TotNFish}, that were eliminated from
#' the population based on their size (\code{len}, \code{wt}, \code{ts}) or 
#' location (\code{f.east}, \code{f.north}, \code{f.d2sh}, \code{f.botdep}, \code{f.wdep}).
#' If you end up with far fewer fish than requested, this can be useful in troubleshooting where the problem might lie.
#' }
#'
#' @details
#'
#' The artificial lake can be imagined as a rectangular subset of a "real" lake.  
#' The east and west boundaries of the artificial lake do not reach the shoreline of the "real" lake, unless \code{BotDepMin} is set to zero.
#' The north and south boundaries of the artificial lake do not ascend to a shoreline, 
#' instead the bottom depth remains constant in the south-north direction (i.e., for a given easting).
#' The angle of the western lake bed is twice as steep as the angle of the eastern lake bed.
#' View the top and side views of the artificial lake in this diagram \href{https://raw.githubusercontent.com/JVAdams/artiFISHal/master/images/LakeFigures.JPG}{[link]}.
#'
#' You may wish to cap the total number of fish at 5 million if your computer has a memory of about 2 GB (2047 MB).  
#' This limit can be increased if you have more memory available in R.
#' You can check the memory available with \code{\link{memory.limit}}.
#'
#' The diagnostic plots produced, if \code{PlotsPdf} is not FALSE, include 
#' scatterplots of 1,000 fish randomly selected from the population, scatterplots of 250 fish randomly selected from each group, 
#' and histograms of the size and spatial distribution of all the fish in the lake.
#'
#' @export
#' @import 				MASS
#' @seealso \code{\link{SampFish}}
#' @references Yule, DL, JV Adams, DM Warner, TR Hrabik, PM Kocovsky, BC Weidel, LG Rudstam, and PJ Sullivan.  2013.  
#' Evaluating analytical approaches for estimating pelagic fish biomass using simulated fish communities. 
#' Canadian Journal of Fisheries and Aquatic Sciences 70:1845-1857.  
#' \emph{http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0072#.U1KYxPldXTQ}
#'
#' @examples
#' \dontrun{
#'
#' # parameters for small (a) and large (A) alewife as input to the simulator
#' fishp <- data.frame(
#' 	G = c("a", "A", "A"), 
#' 	Z = c(50, 140, 140), ZE = c(0.25, 0.2, 0.2), 
#' 	LWC1 = 0.000014, LWC2 = 2.8638, LWCE = 0.18, 
#' 	TSC1 = -64.2, TSC2 = 20.5, TSCE = c(0.02, 0.07, 0.07), 
#' 	PropN = c(0.55, 0.25, 0.20), 
#' 	E = c(NA, 900, 2800), EE = c(NA, 4.5, 0.3), 
#' 	N = NA, NE = NA, 
#' 	WD = c(5, 15, 15), WDE = c(0.5, 0.7, 0.7), 
#' 	D2B = NA, D2BE = NA)
#' 
#' # simulate the fish population
#' res <- SimFish(LakeName="Clear Lake", LkWidth=3000, LkLength=2000, 
#'	BotDepMin=20, BotDepMax=100, FishParam=fishp, TotNFish=1000, Seed=667)
#'
#' # look at the results
#' res$Truth
#' res$LakeInfo
#' res$FishInfo
#' head(res$FishParam)
#' head(res$Fish)
#' res$PropExcluded
#' 
#' }

SimFish <- function(LakeName, LkWidth, LkLength, BotDepMin, BotDepMax, BotDepVertex=2*BotDepMax, 
	FishParam, TotNFish, TSRange=c(-65, -20), PlotsPdf=NA, Seed=NULL) {

	if(!is.na(PlotsPdf) & PlotsPdf!=FALSE) pdf(PlotsPdf, width=9, height=6.5, title="Diagnostics", paper="USr")

	# total number of fish
	FishParam$Nfish <- floor(FishParam$PropN * TotNFish / sum(FishParam$PropN))
	FishParam <- FishParam[FishParam$Nfish > 0, ]

	# minimum bottom depth, maximum bottom depth, and "vertex" bottom depth in m (z direction)
	botdepr <- c(BotDepMin, BotDepMax, BotDepVertex)

	# easting (m), range, (x direction)
	eastr <- c(0, LkWidth)

	# northing (m), range, (y direction)
	northr <- c(0, LkLength)

	# slopes and intercepts of west and east shores
	slopes <- c( (botdepr[3]-botdepr[1])/(eastr[2]/3), (botdepr[1]-botdepr[3])/(2*eastr[2]/3) )
	ints <- c(botdepr[1], -eastr[2]*slopes[2]+botdepr[1])

	# distance from west and east shore excluded from lake in m (x direction)
	d2shr.we <- c(0+ints[1]/slopes[1], -ints[2]/slopes[2]-eastr[2])

	# distance to shore (m), range, (x direction)
	mid.d <- ((eastr[1]-d2shr.we[1]) + (eastr[2]+d2shr.we[2]) ) / 2
	d2shr <- c(min(d2shr.we), d2shr.we[1]+mid.d)


	###  FISH  (targets)  ###

	# check that exactly ONE of these two variables (WD, D2B) were set to NA
	look <- FishParam[, c("WD", "D2B")]
	not1na <- apply(is.na(look), 1, sum) != 1
	if(sum(not1na)>0) stop(paste("Rows ", paste(seq_along(not1na)[not1na], collapse=", "), 
		".  Either a mean water depth or a mean distance to bottom MUST be specified, but NOT BOTH!", sep=""))

	nrowz <- dim(FishParam)[1]
	start.i <- (c(0, cumsum(FishParam$Nfish))+1)[1:nrowz]
	end.i <- cumsum(FishParam$Nfish)
	 
	totfish <- sum(FishParam$Nfish)

	fish <- data.frame(G=rep(FishParam$G, FishParam$Nfish), f.east=NA, f.north=NA, f.d2sh=NA, f.botdep=NA, f.wdep=NA, f.d2bot=NA, 
		len=NA, wt=NA, ts=NA)
	if(!is.null(Seed)) set.seed(Seed)
	for(i in seq(nrowz)) {

		# easting available? no, then yes
		if(is.na(FishParam$E[i])) {
			f.east <- runif(FishParam$Nfish[i], eastr[1], eastr[2])
			f.d2sh <- dfromx(x=f.east, d2shr.we=d2shr.we, eastr=eastr)
			f.botdep <- zfromx(x=f.east, maxz=BotDepMax, eastr=eastr, ints=ints, slopes=slopes)
			} else {
			f.east <- rnorm(FishParam$Nfish[i], FishParam$E[i], FishParam$EE[i]*FishParam$E[i])
			f.d2sh <- dfromx(x=f.east, d2shr.we=d2shr.we, eastr=eastr)
			f.botdep <- zfromx(x=f.east, maxz=BotDepMax, eastr=eastr, ints=ints, slopes=slopes)
			}

		# northing available? no, then yes
		if(is.na(FishParam$N[i])) {
			f.north <- runif(FishParam$Nfish[i], northr[1], northr[2])
			} else {
			f.north <- rnorm(FishParam$Nfish[i], FishParam$N[i], FishParam$NE[i]*FishParam$N[i])
			}

		# distance to bottom available? no, then yes (if not, water depth is)
		if(is.na(FishParam$D2B[i])) {
			f.wdep <- rnorm(FishParam$Nfish[i], FishParam$WD[i], FishParam$WDE[i]*FishParam$WD[i])
			f.d2bot <- f.botdep - f.wdep
			} else {
			f.d2bot <- rnorm(FishParam$Nfish[i], FishParam$D2B[i], FishParam$D2BE[i]*FishParam$D2B[i])
			f.wdep <- f.botdep - f.d2bot
			}

		# generate lengths from gamma distribution
		shape <- 1/FishParam$ZE[i]^2
		scale <- FishParam$Z[i]/shape
		len <- rgamma(FishParam$Nfish[i], shape=shape, scale=scale)
		# predict weight from length-weight regression coefficients and add error
		wt. <- FishParam$LWC1[i]*len^FishParam$LWC2[i]
		wt <- wt. + rnorm(FishParam$Nfish[i], 0, FishParam$LWCE[i]*wt.)
		# predict target strength from length-ts regression coefficients and add error
		ts. <- FishParam$TSC1[i] + FishParam$TSC2[i]*log10(len/10)
		ts <- ts. + rnorm(FishParam$Nfish[i], 0, FishParam$TSCE[i]*abs(ts.))

		fish[start.i[i]:end.i[i], -1] <- cbind(f.east, f.north, f.d2sh, f.botdep, f.wdep, f.d2bot, len, wt, ts)
		}

	rm(i, start.i, end.i, f.east, f.north, f.d2sh, f.botdep, f.wdep, f.d2bot, shape, scale, len, wt., wt, ts., ts)

	bad.lenwt <- fish$len < 0 | fish$wt < 0
	bad.ts <- fish$ts < TSRange[1] | fish$ts > TSRange[2]
	bad.east <- fish$f.east < eastr[1] | fish$f.east > eastr[2]
	bad.north <- fish$f.north < northr[1] | fish$f.north > northr[2]
	bad.wdep <- fish$f.wdep < 0 | fish$f.wdep > fish$f.botdep
	bad.botdep <- fish$f.botdep < 0 | fish$f.botdep > BotDepMax
	bad.d2sh <- fish$f.d2sh < d2shr[1] | fish$f.d2sh > d2shr[2] | is.na(fish$f.d2sh)

	n <- dim(fish)[1]

	# if you are losing a lot of fish, you can use this to try and determine why
	PropExcluded <- c(LengthWeight=sum(bad.lenwt), TS=sum(bad.ts), Easting=sum(bad.east), Northing=sum(bad.north), WaterDepth=sum(bad.wdep), 
		BottomDepth=sum(bad.botdep), D2Shore=sum(bad.d2sh))/n

	# get rid of rows that are beyond the bounds we set or that have other problems
	bad <- bad.d2sh | bad.east | bad.north | bad.botdep | bad.wdep | bad.lenwt | bad.ts
	fish <- fish[!bad, ]

	
	
	###  diagnostic plots  ###
	if(is.na(PlotsPdf) | PlotsPdf!=FALSE) {

		# a random selection of 1,000 fish (total)
		n <- dim(fish)[1]
		pick <- if(n<1001) 1:n else sample(1:n, 1000)
		
		fpick <- fish[pick, ]

		sug <- sort(unique(fpick$G))

		if(is.na(PlotsPdf)) dev.new(w=9, h=6.5)
		par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
		plotblank(eastr/1000, c(-BotDepMax, 0), xlab="Easting  (km)", ylab="Water depth  (m)", main=paste(LakeName, "- Side View"))
		lines(c(0, xfromz(z=rep(BotDepMax-0.01, 2), maxz=BotDepMax, ints=ints, slopes=slopes, shore=0:1), eastr[c(2, 2, 1, 1)])/1000, 
			-c(botdepr[1], rep(BotDepMax, 2), botdepr[1], 0, 0, botdepr[1]))

		for(i in seq(along=sug)) {
			sel <- fpick$G==sug[i]
			text(fpick$f.east[sel]/1000, -fpick$f.wdep[sel], fpick$G[sel], col=i)
			}

		if(is.na(PlotsPdf)) dev.new(w=9, h=6.5)
		par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
		plotblank(eastr/1000, northr/1000, xlab="Easting  (km)", ylab="Northing  (km)", main=paste(LakeName, "- Top View"))
		for(i in seq(along=sug)) {
			sel <- fpick$G==sug[i]
			text(fpick$f.east[sel]/1000, fpick$f.north[sel]/1000, fpick$G[sel], col=i)
			}

		
		
		# a random selection of 250 fish from each group

		rows.g <- split(seq(along=fish$G), fish$G)
		pick <- unlist(lapply(rows.g, function(x) sample(x, min(250, length(x)))))
		
		fpick <- fish[pick, ]

		sug <- sort(unique(fpick$G))

		if(is.na(PlotsPdf)) dev.new(w=9, h=6.5)
		par(mfrow=n2mfrow(length(sug)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
		for(i in seq(along=sug)) {
			sel <- fpick$G==sug[i]
			plotblank(eastr/1000, c(-BotDepMax, 0))
			lines(c(0, xfromz(z=rep(BotDepMax-0.01, 2), maxz=BotDepMax, ints=ints, slopes=slopes, shore=0:1), eastr[c(2, 2, 1, 1)])/1000, 
				-c(botdepr[1], rep(BotDepMax, 2), botdepr[1], 0, 0, botdepr[1]))
			text(fpick$f.east[sel]/1000, -fpick$f.wdep[sel], fpick$G[sel], col=i)
			}
		mtext("Easting  (km)", side=1, outer=TRUE)
		mtext("Water depth  (m)", side=2, outer=TRUE)
		mtext(paste(LakeName, "- Side View"), side=3, outer=TRUE, font=2)

		if(is.na(PlotsPdf)) dev.new(w=9, h=6.5)
		par(mfrow=n2mfrow(length(sug)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
		for(i in seq(along=sug)) {
			sel <- fpick$G==sug[i]
			plotblank(eastr/1000, northr/1000)
			text(fpick$f.east[sel]/1000, fpick$f.north[sel]/1000, fpick$G[sel], col=i)
			}
		mtext("Easting  (km)", side=1, outer=TRUE)
		mtext("Northing  (km)", side=2, outer=TRUE)
		mtext(paste(LakeName, "- Top View"), side=3, outer=TRUE, font=2)

		if(is.na(PlotsPdf)) dev.new(w=9, h=6.5)
		par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
		plotblank(fpick$len, -fpick$f.wdep, xlab="Fish length  (mm)", ylab="Water depth  (m)", main=paste(LakeName, "- Size at Depth"))
		for(i in seq(along=sug)) {
			sel <- fpick$G==sug[i]
			text(fpick$len[sel], -fpick$f.wdep[sel], fpick$G[sel], col=i)
			}

		if(is.na(PlotsPdf)) dev.new(w=9, h=6.5)
		par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
		plotblank(fpick$ts, -fpick$f.wdep, xlab="Target strength  (dB)", ylab="Water depth  (m)", main=paste(LakeName, "- Size at Depth"))
		for(i in seq(along=sug)) {
			sel <- fpick$G==sug[i]
			text(fpick$ts[sel], -fpick$f.wdep[sel], fpick$G[sel], col=i)
			}

		if(is.na(PlotsPdf)) dev.new(w=9, h=6.5)
		par(mfrow=n2mfrow(length(sug)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
		for(i in seq(along=sug)) {
			sel <- fpick$G==sug[i]
			plot(fpick$len[sel], fpick$wt[sel], las=1, xlab="", ylab="")
			mtext(sug[i], side=3, adj=0.1, line=-2, font=2)
			}
		mtext("Length  (mm)", side=1, outer=TRUE)
		mtext("Weight  (mm)", side=2, outer=TRUE)
		mtext(paste(LakeName, "- Length-Weight Relation"), side=3, outer=TRUE, font=2)


		
		# histograms of all fish

		fishhist <- function(x, xlab, title, ...) {
			if(is.na(PlotsPdf)) dev.new(w=9, h=6.5)
			par(mfrow=n2mfrow(length(sug)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
			for(i in seq(along=sug)) {
				sel <- fish$G==sug[i]
				hist(x[sel], nclass=25, col="gray", las=1, xlab="", ylab="", main="", ...)
				box()
				mtext(sug[i], side=3, adj=0.9, line=-2, font=2)
				}
			mtext(xlab, side=1, outer=TRUE)
			mtext("Frequency", side=2, outer=TRUE)
			mtext(paste(LakeName, "-", title), side=3, outer=TRUE, font=2)
		}

		sug <- sort(unique(fish$G))
		fishhist(fish$len, "Length  (mm)", "Length Distribution")
		fishhist(fish$wt, "Weight  (g)", "Weight Distribution")
		fishhist(fish$ts, "Target strength  (dB)", "TS Distribution")
		fishhist(fish$f.east/1000, "Easting  (km)", "Easting Distribution", xlim=eastr/1000)
		fishhist(fish$f.north/1000, "Northing  (km)", "Northing Distribution", xlim=northr/1000)
		fishhist(fish$f.d2sh, "Distance to Shore  (m)", "Distance to Shore Distribution", xlim=d2shr)
		fishhist(fish$f.wdep, "Water Depth  (m)", "Water Depth Distribution", xlim=c(0, BotDepMax))
		fishhist(fish$f.d2bot, "Distance to Bottom  (m)", "Distance to Bottom Distribution", xlim=c(0, BotDepMax))
		fishhist(fish$f.botdep, "Bottom Depth  (m)", "Bottom Depth Distribution", xlim=c(0, BotDepMax))

		if(!is.na(PlotsPdf)) graphics.off()

		}

	# total number and weight of each species in population
	truth <- cbind(n=table(fish$G), kg=tapply(fish$wt, fish$G, sum)/1000)
	truth <- as.data.frame(rbind(truth, Total=apply(truth, 2, sum, na.rm=TRUE)))

	LkArea <- LkWidth * LkLength / 10000
	truth$nperha <- truth$n/LkArea
	truth$kgperha <- truth$kg/LkArea

	print(truth)

	# output selected objects for use in sampling programs
	list(
		Truth = truth, 
		LakeInfo = list(LakeName=LakeName, LkWidth=LkWidth, LkLength=LkLength, 
			BotDepMin=BotDepMin, BotDepMax=BotDepMax, BotDepVertex=BotDepVertex,
			ints=ints, slopes=slopes, d2shr.we=d2shr.we), 
		FishInfo = list(TotNFish=TotNFish, TSRange=TSRange, Seed=Seed), 
		FishParam = FishParam,
		Fish = fish,
		PropExcluded = PropExcluded
		)

}
