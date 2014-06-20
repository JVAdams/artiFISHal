#' Combine Acoustic and Midwater Trawl Survey Data
#'
#' Combine survey data from acoustic transects and midwater trawl tows (created by \code{\link{SampFish}}).  
#' Apply availability to the acoustic data and catchability (availability and selectivity) to the midwater trawl catch.
#'
#' @param SimPop	 	A list with elements \code{LakeInfo}, \code{FishInfo}, \code{FishParam}, \code{FishPop}, typically output from \code{\link{SimFish}}.
#' See \code{\link{SimFish}} for details on the list elements.
#' @param AcMtSurv	 	A list with elements \code{Targets}, \code{AcSummaryCell}, \code{AcSummaryColumn}, \code{MtCatch}, 
#' typically output from \code{\link{SampFish}}.
#' @param SelecParam	A data frame with 6 columns in which each row provides the midwater trawl selectivity parameters for
#' a given fish group and mesh panel zone.
#' All columns must be completely filled in (no missing values).
#' Selectivity is assumed to be 100% for any group-zone combination not represented as a row in the data frame.
#' For 100% selectivity of small fish, use MtL50Small = -Inf and any slope.  
#' For 100% selectivity of large fish, use MtL50Large = Inf and any slope. 
#' Column names and descriptions:
#' \itemize{
#'   \item \code{G} = character, a one-letter nickname for the group (e.g., fish species and lifestage) used in plotting
#'   \item \code{Zone} = character, mesh panel zone, one of "mouth", "middle", "aft", or "cod"
#'   \item \code{MtL50Small} = the length (in mm) at which small fish have a 50% probability of being captured by the trawl
#'   \item \code{MtSlopeSmall} = the (inverse) slope at which small fish probability of capture increases with length, smaller values are steeper
#'   \item \code{MtL50Large} = the length (in mm) at which large fish have a 50% probability of being captured by the trawl
#'   \item \code{MtSlopeLarge} = the (absolute value of the inverse) slope at which large fish probability of capture decreases with length, smaller values are steeper
#' }
#' @param PanelProps 	A numeric vector of length 4, size of the different mesh panel zones of the midwater trawl, 
#' mouth (outermost), middle, aft, and cod (inner), default c(0.4, 0.3, 0.2, 0.1).
#' Sizes are expressed as proportions of the distance from the outer edge of the trawl to the trawl center in both the vertical and horizontal directions,
#' and they should add up to 1.  Use \code{\link{ViewZones}} to visualize the mesh panel zones.
#' @param AcExcl		A numeric vector of length 2, depth of acoustic "dead" zones at the surface and at the bottom (in m), 
#' default of c(0, 0) represents 100% acoustic availability of fish.
#' @param MtExcl		A numeric vector of length 2, depth of zones unfishable with the midwater trawl at the surface and at the bottom (in m), 
#' default of c(0, 0) represents 100% midwater trawl availability of fish.
#' @param Seed			An integer scalar, starting seed for stochasticity incorporated in acoustic and midwater trawl catchability.  
#' Use \code{Seed} to ensure the same individual fish are included in the surveys with each call to \code{CatchComb}.  
#' Otherwise, if set to NULL, the default, a random seed is used, resulting in a different fish selection with each call to \code{CatchComb}.  

#'
#' @return				A data frame with estimated fish density (in number per ha) and biomass (in kg per ha) for each sampling event and group (species, lifestage).
#'
#' @details
#'
#' A classification tree is used to relate the catch composition of the midwater trawl to the location of the trawl in the lake 
#' (e.g., MTReast, ACnorth, MTRd2sh, MTRbdep).  This tree is then used to assign a single midwater trawl catch to each acoustic cell
#' (interval x layer), such that the estimated acoustic densities can be assigned to specific fish groups (species, life stages).  See,
#' for example, Yule et al. (2013).
#'
#' @export
#' @import				rpart
#' @seealso \code{\link{SimFish}}, \code{\link{SampFish}}, \code{\link{ViewZones}}, \code{\link{TuneSelec}}.
#' @references 
#' 
#' Yule, DL, JV Adams, DM Warner, TR Hrabik, PM Kocovsky, BC Weidel, LG Rudstam, and PJ Sullivan.  2013.  
#' \href{http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0072#.U1KYxPldXTQ}{Evaluating 
#' analytical approaches for estimating pelagic fish biomass using simulated fish communities}. 
#' Canadian Journal of Fisheries and Aquatic Sciences 70:1845-1857.
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
#' 	FishParam=fishp, TotNFish=50000)
#'
#' # survey the population
#' surv <- SampFish(SimPop=res, NumEvents=2, AcNum=5, AcInterval=3000, AcLayer=10, AcAngle=7, MtNum=25, MtHt=10, MtWd=10, MtLen=200)
#'
#' selec <- data.frame(
#' 	G = c("A", "a", "A", "a", "A", "a"), 
#' 	Zone = c("mouth", "mouth", "middle", "middle", "aft", "aft"), 
#' 	MtL50Small = c(100, 100, 60, 60, 30, 30), 
#' 	MtSlopeSmall = c(40, 40, 30, 30, 20, 20), 
#' 	MtL50Large = c(180, 180, Inf, Inf, Inf, Inf), 
#' 	MtSlopeLarge = c(20, 20, 100, 100, 100, 100))
#'
#' AcMtEst(SimPop=res, AcMtSurv=surv, Seed=927)
#' AcMtEst(SimPop=res, AcMtSurv=surv, SelecParam=selec, AcExcl=c(5, 10), MtExcl=c(2, 2), Seed=204)
#'
AcMtEst <- function(SimPop, AcMtSurv, PanelProps=c(0.4, 0.3, 0.2, 0.1), SelecParam=NULL, AcExcl=c(0, 0), MtExcl=c(0, 0), Seed=NULL) {

# SimPop=res
# AcMtSurv=surv
# SelecParam=selec
# PanelProps=c(0.4, 0.3, 0.2, 0.1)
# AcExcl=c(5, 10)
# MtExcl=c(0, 0)
# Seed=245

	if(!is.null(Seed)) set.seed(Seed)
	mtr <- surv$MtCatch
	srv <- surv$SurvParam
	ac <- surv$Targets

	# check validity of the trawl zone proportions that were input
	names(PanelProps) <- c("mouth", "middle", "aft", "cod")
	panelprops <- rev(PanelProps)
	if(round(sum(panelprops), 0.0000001) != 1) stop("cod proportions should sum to 1")

	# availability function, = 0 at surface and bottom and = 1 in the middle
	dblcut <- function(wdep, surfacecut, d2bot, bottomcut) {
		as.numeric(wdep > surfacecut & d2bot > bottomcut)
		}


	# acoustic availability
	ac$keep <- with(ac, dblcut(wdep=f.wdep, surfacecut=AcExcl[1], d2bot=f.d2bot, bottomcut=AcExcl[2]))

	# summarize by cell (interval x layer)
	acs <- AcSmry(AcTarg=ac[ac$keep==1, ], LakeInfo=SimPop$LakeInfo, SurvParam=AcMtSurv$SurvParam)$AcCell


	# midwater trawl availability

	if(is.null(SelecParam)) {
		SelecParam <- data.frame(
			G = character(),
			Zone = character(),
			MtL50Small = numeric(0),
			MtSlopeSmall = numeric(0),
			MtL50Large = numeric(0),
			MtSlopeLarge = numeric(0))
		}

	# check validity of zones
	suz <- c("mouth", "middle", "aft", "cod")
	uz <- unique(SelecParam$Zone)
	badzones <- setdiff(uz, suz)
	if(length(badzones) > 0) stop('Zones must be one of "mouth", "middle", "aft", or "cod".')

	# check for missings
	missings <- sum(is.na(SelecParam))
	if(missings > 0) stop("SelectParam data frame may not have any missing values.")

	# fill in 100% selectivities for group-zones with no parameters
	sug <- sort(unique(AcMtSurv$MtCatch$G))
	full <- expand.grid(G=sug, Zone=suz)
	selec2 <- merge(SelecParam, full, all=TRUE)
	sel100 <- is.na(selec2$MtL50Small)
	selec2$MtL50Small[sel100] <- -Inf
	selec2$MtSlopeSmall[sel100] <- 100
	selec2$MtL50Large[sel100] <- Inf
	selec2$MtSlopeLarge[sel100] <- 200

	# for each fish, determine its maximum vertical or horizontal distance from the center of the trawl
	# as a proportion of the trawl dimensions
	mtr$maxdist <- with(mtr, pmax(abs(f.wdep - MTRwdep)/srv["MtHt"], abs(f.north - ACnorth)/srv["MtWd"]))
	# use this distance to assign each fish to a zone of the trawl
	mtr$Zone <- cut(mtr$maxdist, breaks=c(0, cumsum(panelprops)), include.lowest=TRUE, labels=names(panelprops))
	mtrsel <- merge(mtr, selec2, all.x=TRUE)

	# Think about trawl availability ... when we are cutting off trawls near the surface or the bottom, shouldn't this constraint happen
	# during the survey itself, where trawls that encompass those "dead" zones can be eliminated?

	mtrsel$p.avail <- with(mtrsel, dblcut(wdep=f.wdep, surfacecut=MtExcl[1], d2bot=f.d2bot, bottomcut=MtExcl[2]))
	mtrsel$p.selec <- with(mtrsel, logit2(x=len, x50a=MtL50Small, slopea=MtSlopeSmall, x50b=MtL50Large, slopeb=-MtSlopeLarge))
	mtrsel$p.catch <- mtrsel$p.avail*mtrsel$p.selec

	# apply catchability (selectivity AND availability) functions to "perfect" MTR catch
	mtrsel$keep <- sapply(mtrsel$p.catch, function(p) sample(0:1, size=1, replace=TRUE, prob=c(1-p, p)))


	sue <- sort(unique(acs$Event))
	results <- expand.grid(G=sug, Event=sue, nperha=NA, kgperha=NA)

	for(k in sue) {

		# subset the MT data
		mtk <- mtrsel[mtrsel$Event==k & mtrsel$keep==1, ]

		# only do these calculations if there were "keep" fish in the midwater trawl
		# without "keep" fish, no species-specific density and biomass can be estimated
		if(dim(mtk)[1] > 0) {

			# subset the AC data
			ack <- acs[acs$Event==k, ]

			# make the variable names in mtk the same as in those in ack for tree prediction
			names(mtk)[match(c("MTReast", "MTRd2sh", "MTRbdep", "MTRwdep", "MTRd2bot"), names(mtk))] <- c("interval", "d2sh", "botdep", "layer", "d2bot")

			# fit a classification tree to the MTR data for Event k
			treek <- rpart(as.factor(G) ~ interval + ACnorth + d2sh + botdep + layer + d2bot, data=mtk, control=list(cp=0.05, minsplit=10, minbucket=5))

			# mean density for each species
			# suffixes:  e = event, a = ac transect, i = interval, l = layer
			# use the fitted tree to predict species composition of each AC interval/layer
			pred.props <- predict(treek, newdata=ack)
			dens.eail <- ack$nperha*pred.props
			dens.eai <- aggregate(dens.eail, ack[, c("ACid", "interval")], sum)
			dens.e <- apply(dens.eai[, sug], 2, mean)

			# mean biomass for each species
			# calculate the mean weight of each group at each node of the fitted tree
			mwt <- tapply(mtk$wt, list(treek$where, mtk$G), mean)
			mwt[is.na(mwt)] <- 0

			# proportions corresponding to each node
			pred.node <- prednode(treek, newdata=ack)

			# mean weigth expanded to full dimensions of ack, by matching node number
			mwt.eail <- mwt[match(pred.node, row.names(mwt)), ]
			bio.eail <- mwt.eail * dens.eail
			bio.eai <- aggregate(bio.eail, ack[, c("ACid", "interval")], sum)
			bio.e <- apply(bio.eai[, sug], 2, mean)

			for(g in seq(sug)) {
				sel <- results$Event==k & results$G==sug[g]
				results$nperha[sel] <- dens.e[sug[g]]
				results$kgperha[sel] <- bio.e[sug[g]]/1000
				}
			}
		}

	results$nperha[is.na(results$nperha)] <- 0
	results$kgperha[is.na(results$kgperha)] <- 0
	results[, c("Event", "G", "nperha", "kgperha")]
#	write.csv(results, paste0(subdir, "/ResultsSummary-lake", sel.lk, "-run", sel.run, ".csv"), row.names=FALSE)
	}
