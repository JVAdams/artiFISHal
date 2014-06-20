#' Summarize Acoustic Survey Data
#'
#' Summarize acoustic survey data by interval and layer.
#'
#' @param AcTarg	 	A data frame list with information on fish detected in an acoustics survey, specifically the \code{Targets} data frame
#' output from \code{\link{SampFish}}.  Each row represents a single target, columns describe the specific location of fish in the lake and their target strengths.  
#' @param LakeInfo		A list with the lake inputs supplied as arguments to \code{\link{SimFish}} as well as a few additional objects.
#' @param SurvParam		A named vector with the survey inputs supplied as arguments to \code{\link{SampFish}}.
#'
#' @return 				A list with 2 elements.  
#' \itemize{
#'   \item \code{AcCell} is a data frame with information on the acoustic survey summarized by acoustic interval and layer.
#'   \item \code{AcColumn} is a data frame with information on the acoustic survey summarized by acoustic interval.
#' }
#'
#' @details
#'
#' Acoustic intervals are identified by the easting (in m) of their midpoint.
#' Acoustic layers are indentified by the water depth (in m) of their midpoint.
#'
#' A weighting variable, range weight (Yule 2000), is used to account for different volumes of water
#' sampled in the acoustic survey as a function of the distance from the transducer (in m) and the transducer half angle (0.5 * \code{AcAngle}).
#' The sum of the range weights is reported as \code{sum.rw} in the \code{AcCell} and \code{AcColumn} data frames.
#'
#' @export
#' @import 				MASS
#' @seealso \code{\link{SampFish}}
#' @references 
#' 
#' Yule, DL.  2000.  
#' \href{http://www.tandfonline.com/doi/abs/10.1577/1548-8675(2000)020%3C0759%3ACOHAAP%3E2.3.CO%3B2}{Comparison 
#' of horizontal acoustic and purse-seine estimates of salmonid densities and sizes in eleven Wyoming waters}. 
#' North American Journal of Fisheries Management 20:759-775.
#'
#' @examples
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
#' res <- SimFish(LakeName="Clear Lake", LkWidth=3000, LkLength=2000, BotDepMin=20, BotDepMax=100, 
#' 	FishParam=fishp, TotNFish=50000, Seed=667)
#'
#' # survey the population
#' surv <- SampFish(SimPop=res, NumEvents=2, AcNum=5, AcInterval=3000, AcLayer=10, AcAngle=7, MtNum=25, MtHt=10, MtWd=10, MtLen=200, Seed=545)
#'
#' AcSmry(AcTarg=surv$Targets, LakeInfo=res$LakeInfo, SurvParam=surv$SurvParam)
#'
AcSmry <- function(AcTarg, LakeInfo, SurvParam) {

	# AcTarg=surv$Targets
	# LakeInfo=res$LakeInfo
	# SurvParam=surv$SurvParam

	# add "range weight" to the AcTarg data (Yule 2000)
	# a weighting variable to account for different volumes sampled as a function of range (dist. from ducer in m) and ducer half angle
	AcTarg$rng.wt <- 1/(2*AcTarg$f.wdep*tan(SurvParam["AcAngle"]))

	# summarize acoustic data by interval and layer
	AC2 <- aggregate(AcTarg[, "rng.wt"], AcTarg[, c("Event", "ACid", "ACnorth", "interval", "layer")], sum)
	names(AC2)[names(AC2)=="x"] <- "sum.rw"
	AC2$nperha <- 10000 * AC2$sum.rw/SurvParam["AcInterval"]

	# create a matrix of all possible interval-by-layer combinations for each ACid
	# this will be used to ensure that interval-by-layers with no fish are included in the summaries

	eastr <- c(0, LakeInfo$LkWidth)
	intbrks <- seq(eastr[1], eastr[2] + SurvParam["AcInterval"] - 1, SurvParam["AcInterval"])
	intmids <- intbrks[-1] - SurvParam["AcInterval"]/2

	laybrks <- seq(0, LakeInfo$BotDepMax + SurvParam["AcLayer"] - 1, SurvParam["AcLayer"])
	laymids <- laybrks[-1] - SurvParam["AcLayer"]/2

	# each AC transect goes over the same depth profile, so we can determine the max depth for each interval
	# then see whether the max layer included

	easts <- seq(eastr[1], eastr[2], length=1000)
	depth.contour <- zfromx(x=easts, maxz=LakeInfo$BotDepMax, eastr=eastr, ints=LakeInfo$ints, slopes=LakeInfo$slopes)
	depth.cont.int <- intmids[cut(easts, include.lowest=TRUE, breaks=intbrks, labels=FALSE)]
	all.maxes <- tapply(depth.contour, depth.cont.int, max)
	max.lays <- laymids[cut(all.maxes, include.lowest=TRUE, breaks=laybrks, labels=FALSE)]

	# full matrix of all ACids, all intervals, and all layers
	sua <- sort(unique(AcTarg$ACid))
	full.mat <- expand.grid(layer=seq(laymids[1], max(max.lays), SurvParam["AcLayer"]), interval=intmids, ACid=1:(SurvParam["NumEvents"]*SurvParam["AcNum"]))
	full.mat$Event <- recode(full.mat$ACid, sua, AcTarg$Event[match(sua, AcTarg$ACid)])

	sub.mat <- merge(data.frame(interval=intmids, max.layer=max.lays), full.mat, all=TRUE)
	sub.mat <- sub.mat[sub.mat$layer <= (sub.mat$max.layer + 0.001), c("Event", "ACid", "interval", "layer")]

	ACsmryIL <- merge(AC2, sub.mat, all=TRUE)
	ACsmryIL$nperha[is.na(ACsmryIL$nperha)] <- 0
	ACsmryIL$sum.rw[is.na(ACsmryIL$sum.rw)] <- 0
	ACsmryIL$ACnorth <- recode(ACsmryIL$ACid, sua, AcTarg$ACnorth[match(sua, AcTarg$ACid)])

	rm(easts, depth.contour, depth.cont.int, all.maxes, max.lays, full.mat)

	d2shr.we <- c(0 + LakeInfo$ints[1]/LakeInfo$slopes[1], -LakeInfo$ints[2]/LakeInfo$slopes[2] - eastr[2])

	ACsmryIL$botdep <- zfromx(x=ACsmryIL$interval, maxz=LakeInfo$BotDepMax, eastr=eastr, ints=LakeInfo$ints, slopes=LakeInfo$slopes)
	ACsmryIL$d2bot <- ACsmryIL$botdep - ACsmryIL$layer
	ACsmryIL$d2sh <- dfromx(x=ACsmryIL$interval, d2shr.we=d2shr.we, eastr=eastr)

	# ACsmryIL <- merge(ACsmryIL, sub.mat, all=TRUE)
	# ACsmryIL$nperha[is.na(ACsmryIL$nperha)] <- 0

	# summarize acoustic data by interval only
	ACsmryI <- aggregate(ACsmryIL[, c("sum.rw", "nperha")], 
		ACsmryIL[, c("Event", "ACid", "ACnorth", "interval", "botdep", "d2sh")], sum)
	ACsmryI <- ACsmryI[order(ACsmryI$Event, ACsmryI$ACid, ACsmryI$interval), ]

	acsumil.out.vars <- c("Event", "ACid", "ACnorth", "interval", "layer", "d2sh", "botdep", "d2bot", "sum.rw", "nperha")

	acsumi.out.vars <-  c("Event", "ACid", "ACnorth", "interval",          "d2sh", "botdep",          "sum.rw", "nperha")

	list(AcCell=ACsmryIL[, acsumil.out.vars], AcColumn=ACsmryI[, acsumi.out.vars] )

	}
