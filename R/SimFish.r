#' Simulate a Fish Population
#'
#' Create a simulated population of pelagic fish in an artificial lake.
#' @param LakeName 		A character scalar, full name for artificial lake to be used in plot titles.
#' @param LakeNick 		A character scalar, nickname for the artificial lake to be used in file names.
#' @param LkWidth		A numeric scalar, the width of the lake in the west-east direction (in m), default 30,000.
#' @param LkLength		A numeric scalar, the length of the lake in the south-north direction (in m), default 20,000.
#' @param BotDepMin		A numeric scalar, the minimum bottom depth of the lake, at both the west and east shorelines (in m), default 20.
#' @param BotDepMax		A numeric scalar, the maximum bottom depth of the lake (in m).
#' @param BotDepVertex	A numeric scalar, the vertical distance from the surface to the "vertex" of the lake bottom (in m), default \code{2*BotDepMax}.  
#' The "vertex" of the lake bottom is the point at which the angled lake beds along the west and east shores would intersect,
#' were they not cut off first by the specified \code{BotDepMax}.  
#' View this \href{https://raw.githubusercontent.com/JVAdams/artiFISHal/master/LakeFigures.JPG}{figure} for a diagram of the artificial lake.
#' @param FishParam		A data frame with 19 columns in which each row describes a sub-population of fish to be placed in the artificial lake.
#' The first 11 columns must be completely filled in (no missing values).
#' The last 8 columns may have some missing values.  
#' However, in each row, \strong{either} water depth (WD and WDE) \strong{or} distance to bottom (D2B and D2BE) \strong{must} be filled in, 
#' but \strong{not both}.
#' Column names and descriptions:
#' \itemize{
#'   \item \code{G} = character, a one-letter nickname for the group (e.g., fish species and lifestage) that will be used in plotting
#'   \item \code{Group} = character, full name of the group (e.g., fish species and lifestage)
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
#' Memory on your computer will limit the size of \code{TotNFish}.
#' @param TSRange		A numeric vector of length 2, the range of target strengths to use for the fish (in db), default c(-65, -20).
#' @param plots.pdf		A character scalar, name of pdf file to store the diagnostic plots in.  If NA, the default, 
#' plots will be displayed on the screen instead.  If FALSE, no plots will be created.
#' @param Seed			An integer scalar, starting seed for stochasticity incorporated in fish population generation.  
#' Provide an integer to ensure the same population is generated with each call to \code{\link{SimFish}()}.  
#' Otherwise, if set to NULL, the default, a random seed will be used, resulting in a different population with each call to \code{\link{SimFish}()}.  
#' @return 				A list with five elements, (1) \code{Truth}, a data frame (or matrix?) with the total number and weight 
#' of each fish group in the population, 
#' (2) \code{LakeInfo}, a list with the lake inputs supplied as arguments to \code{SimFish} as well as a few additional objects 
#' which will be useful in surveying the population,
#' \itemize{
#'   \item \code{ints} = a numeric vector of length 2, the intercepts of the angled lake beds along the west and east shores of the lake (in m)
#'   \item \code{slopes} = a numeric vector of length 2, the slopes of the angled lake beds along the west and east shores of the lake (unitless)
#'   \item \code{d2shr.we} = a numeric vector of length 2, distance (in the "x" direction) from west and east shores excluded from lake (in m)
#' }
#' (3) \code{FishInfo}, a list with the fish inputs supplied as arguments to \code{SimFish}, 
#' (4) \code{FishParam}, the data frame supplied as an argument to \code{SimFish}, and
#' (5) \code{FishPop}, a data frame in which each row is a fish, and the 10 columns describe the fish group (\code{sp}), location 
#' (easting \code{f.east}, northing \code{f.north}, distance to shore \code{f.d2sh}, bottom depth \code{f.botdep}, 
#' water depth \code{f.fdep}, distance to bottom \code{f.d2bot}, all in m), and 
#' fish size (total length in mm \code{len}, weight in g \code{wt}, and target strength in db, \code{ts}).
#' @details
#'
#' The artificial lake can be imagined as a rectangular subset of a "real" lake.  
#' The east and west boundaries of the artifical lake do not reach the shoreline of the "real" lake, unless \code{BotDepMin} is set to zero.
#' The north and south boundaries of the artificial lake do not ascend to a shoreline, 
#' instead the bottom depth remains constant in the south-north direction (i.e., for a given easting).
#' The angle of the western lake bed is twice as steep as the angle of the eastern lake bed.
#' View the top and side views of the artificial lake in this \href{https://raw.githubusercontent.com/JVAdams/artiFISHal/master/LakeFigures.JPG}{diagram}.
#'
#' You may wish to cap the total number of fish at 5 million if your computer has a memory of about 2 GB (2047 MB).  
#' This limit can be increased if you have more memory available in R.
#' You can check the memory available with \code{\link{memory.limit}()}.
#'
#' @export
#' @import 				MASS

SimFish <- function(LakeName, LakeNick, LkWidth=30000, LkLength=20000, BotDepMin=20, BotDepMax, BotDepVertex=2*BotDepMax, 
	FishParam, TotNFish, TSRange=c(-65, -20), plots.pdf=NA, Seed=NULL) {

	if(!is.na(plots.pdf) & plots.pdf!=FALSE) pdf(plots.pdf, width=9, height=6.5, title="Diagnostics", paper="USr")

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
		".  Either a mean fishing depth or a mean distance to bottom MUST be specified, but NOT BOTH!", sep=""))
#	rm(look, not1na)

	nrowz <- dim(FishParam)[1]
	start.i <- (c(0, cumsum(FishParam$Nfish))+1)[1:nrowz]
	end.i <- cumsum(FishParam$Nfish)
	 
	totfish <- sum(FishParam$Nfish)

	attach(FishParam)

	fish <- data.frame(sp=rep(G, Nfish), f.east=NA, f.north=NA, f.d2sh=NA, f.botdep=NA, f.fdep=NA, f.d2bot=NA, len=NA, wt=NA, ts=NA)
	if(!is.null(Seed)) set.seed(Seed)
	for(i in seq(nrowz)) {

		# easting available? no, then yes
		if(is.na(E[i])) {
			f.east <- runif(Nfish[i], eastr[1], eastr[2])
			f.d2sh <- dfromx(x=f.east, d2shr.we=d2shr.we, eastr=eastr)
			f.botdep <- zfromx(x=f.east, maxz=BotDepMax, eastr=eastr, ints=ints, slopes=slopes)
			} else {
			f.east <- rnorm(Nfish[i], E[i], EE[i]*E[i])
			f.d2sh <- dfromx(x=f.east, d2shr.we=d2shr.we, eastr=eastr)
			f.botdep <- zfromx(x=f.east, maxz=BotDepMax, eastr=eastr, ints=ints, slopes=slopes)
			}

		# northing available? no, then yes
		if(is.na(N[i])) {
			f.north <- runif(Nfish[i], northr[1], northr[2])
			} else {
			f.north <- rnorm(Nfish[i], N[i], NE[i]*N[i])
			}

		# distance to bottom available? no, then yes (if not, fishing depth is)
		if(is.na(D2B[i])) {
			f.fdep <- rnorm(Nfish[i], WD[i], WDE[i]*WD[i])
			f.d2bot <- f.botdep - f.fdep
			} else {
			f.d2bot <- rnorm(Nfish[i], D2B[i], D2BE[i]*D2B[i])
			f.fdep <- f.botdep - f.d2bot
			}

		# generate lengths from gamma distribution
		shape <- 1/ZE[i]^2
		scale <- Z[i]/shape
		len <- rgamma(Nfish[i], shape=shape, scale=scale)
		# predict weight from length-weight regression coefficients and add error
		wt. <- LWC1[i]*len^LWC2[i]
		wt <- wt. + rnorm(Nfish[i], 0, LWCE[i]*wt.)
		# predict target strength from length-ts regression coefficients and add error
		ts. <- TSC1[i] + TSC2[i]*log10(len/10)
		ts <- ts. + rnorm(Nfish[i], 0, TSCE[i]*abs(ts.))

		fish[start.i[i]:end.i[i], -1] <- cbind(f.east, f.north, f.d2sh, f.botdep, f.fdep, f.d2bot, len, wt, ts)
		}


#	rm(i, start.i, end.i, f.east, f.north, f.d2sh, f.botdep, f.fdep, f.d2bot, shape, scale, len, wt., wt, ts., ts)

	bad.d2sh <- fish$f.d2sh < d2shr[1] | fish$f.d2sh > d2shr[2] | is.na(fish$f.d2sh)
	bad.east <- fish$f.east < eastr[1] | fish$f.east > eastr[2]

	bad.north <- fish$f.north < northr[1] | fish$f.north > northr[2]
	bad.botdep <- fish$f.botdep < 0 | fish$f.botdep > BotDepMax
	bad.fdep <- fish$f.fdep < 0 | fish$f.fdep > fish$f.botdep
	bad.lenwt <- fish$len < 0 | fish$wt < 0
	bad.ts <- fish$ts < TSRange[1] | fish$ts > TSRange[2]
	n <- dim(fish)[1]

	# if you are losing a lot of fish, you can use this to try and determine why
	if(FALSE) {
	sum(bad.d2sh)/n
	sum(bad.east)/n
	sum(bad.north)/n
	sum(bad.botdep)/n
	sum(bad.fdep)/n
	sum(bad.lenwt)/n
	sum(bad.ts)/n
	}

	# get rid of rows that are beyond the bounds we set or that have other problems
	bad <- bad.d2sh | bad.east | bad.north | bad.botdep | bad.fdep | bad.lenwt | bad.ts
	fish <- fish[!bad, ]

#	rm(bad.d2sh, bad.east, bad.north, bad.botdep, bad.fdep, bad.lenwt, bad.ts, bad)


	detach(FishParam)



	###  diagnostic plots  ###
	if(is.na(plots.pdf) | plots.pdf!=FALSE) {


	# a random selection of 1,000 fish (total)
	n <- dim(fish)[1]
	pick <- if(n<1001) 1:n else sample(1:n, 1000)
	attach(fish[pick, ])
	sus <- sort(unique(sp))

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
	plot(f.east/1000, -f.fdep, type="n", xlim=eastr/1000, ylim=c(-BotDepMax, 0), 
		xlab="Easting  (km)", ylab="Water depth  (m)", main=paste(LakeName, "- Side View"))
	lines(c(0, xfromz(z=rep(BotDepMax-0.01, 2), maxz=BotDepMax, ints=ints, slopes=slopes, shore=0:1), eastr[c(2, 2, 1, 1)])/1000, 
		-c(botdepr[1], rep(BotDepMax, 2), botdepr[1], 0, 0, botdepr[1]))

	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		text(f.east[sel]/1000, -f.fdep[sel], sp[sel], col=i)
		}

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
	plot(f.east/1000, f.north/1000, type="n", xlim=eastr/1000, ylim=northr/1000, 
		xlab="Easting  (km)", ylab="Northing  (km)", main=paste(LakeName, "- Top View"))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		text(f.east[sel]/1000, f.north[sel]/1000, sp[sel], col=i)
		}

	detach(fish[pick, ])



	# a random selection of 250 fish FROM EACH SPECIES

	rows.sp <- split(seq(along=fish$sp), fish$sp)
	pick <- unlist(lapply(rows.sp, function(x) sample(x, min(250, length(x)))))
	attach(fish[pick, ])
	sus <- sort(unique(sp))

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		plot(f.east/1000, -f.fdep, type="n", xlim=eastr/1000, ylim=c(-BotDepMax, 0), xlab="", ylab="")
		lines(c(0, xfromz(z=rep(BotDepMax-0.01, 2), maxz=BotDepMax, ints=ints, slopes=slopes, shore=0:1), eastr[c(2, 2, 1, 1)])/1000, 
			-c(botdepr[1], rep(BotDepMax, 2), botdepr[1], 0, 0, botdepr[1]))
		text(f.east[sel]/1000, -f.fdep[sel], sp[sel], col=i)
		}
	mtext("Easting  (km)", side=1, outer=TRUE)
	mtext("Water depth  (m)", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Side View"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		plot(f.east/1000, f.north/1000, type="n", xlim=eastr/1000, ylim=northr/1000, xlab="", ylab="")
		text(f.east[sel]/1000, f.north[sel]/1000, sp[sel], col=i)
		}
	mtext("Easting  (km)", side=1, outer=TRUE)
	mtext("Northing  (km)", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Top View"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
	plot(len, -f.fdep, type="n", xlab="Fish length  (mm)", ylab="Water depth  (m)", 
		main=paste(LakeName, "- Size at Depth"))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		text(len[sel], -f.fdep[sel], sp[sel], col=i)
		}

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
	plot(ts, -f.fdep, type="n", xlab="Target strength  (dB)", ylab="Water depth  (m)", 
		main=paste(LakeName, "- Size at Depth"))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		text(ts[sel], -f.fdep[sel], sp[sel], col=i)
		}

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		plot(len[sel], wt[sel], xlab="", ylab="")
		mtext(sus[i], side=3, adj=0.1, line=-2, font=2)
		}
	mtext("Length  (mm)", side=1, outer=TRUE)
	mtext("Weight  (mm)", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Length-Weight Relation"), side=3, outer=TRUE, font=2)

	detach(fish[pick, ])



	# histograms of all fish

	attach(fish)

	sus <- sort(unique(sp))

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(len[sel], nclass=25, col="gray", xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Length  (mm)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Length Distribution"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(wt[sel], nclass=25, col="gray", xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Weight  (g)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Weight Distribution"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(ts[sel], nclass=25, col="gray", xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Target strength  (dB)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- TS Distribution"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(f.east[sel]/1000, nclass=25, col="gray", xlim=eastr/1000, xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Easting  (km)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Easting Distribution"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(f.north[sel]/1000, nclass=25, col="gray", xlim=northr/1000, xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Northing  (km)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Northing Distribution"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(f.d2sh[sel], nclass=25, col="gray", xlim=d2shr, xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Distance to Shore  (m)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Distance to Shore Distribution"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(f.fdep[sel], nclass=25, col="gray", xlim=c(0, BotDepMax), xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Fishing Depth  (m)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Fishing Depth Distribution"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(f.d2bot[sel], nclass=25, col="gray", xlim=c(0, BotDepMax), xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Distance to Bottom  (m)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Distance to Bottom Distribution"), side=3, outer=TRUE, font=2)

	if(is.na(plots.pdf)) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(i in seq(along=sus)) {
		sel <- sp==sus[i]
		hist(f.botdep[sel], nclass=25, col="gray", xlim=c(0, BotDepMax), xlab="", ylab="", main="")
		box()
		mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
		}
	mtext("Bottom Depth  (m)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(LakeName, "- Bottom Depth Distribution"), side=3, outer=TRUE, font=2)

	detach(fish)



	if(!is.na(plots.pdf)) graphics.off()
	}


	# total number and weight of each species in population
	truth <- cbind(n=table(fish$sp), kg=tapply(fish$wt, fish$sp, sum)/1000)
	truth <- as.data.frame(rbind(truth, Total=apply(truth, 2, sum, na.rm=TRUE)))

	LkArea <- LkWidth * LkLength / 10000
	truth$dens <- truth$n/LkArea
	truth$bio <- truth$kg/LkArea

	print(truth)

	# output selected objects for use in sampling programs
	list(
		truth = truth, 
		lkinfo = list(LakeName=LakeName, LakeNick=LakeNick, LkWidth=LkWidth, LkLength=LkLength, 
			BotDepMin=BotDepMin, BotDepMax=BotDepMax, BotDepVertex=BotDepVertex,
			ints=ints, slopes=slopes, d2shr.we=d2shr.we), 
		fishinfo = list(TotNFish=TotNFish, TSRange=TSRange, Seed=Seed), 
		FishParam=FishParam,
		fish = fish
		)

}
