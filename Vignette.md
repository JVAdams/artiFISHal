artiFISHal Vignette
===================

The following is an example using the functions in the **artiFISHal** package.  
See [README.md](https://github.com/JVAdams/artiFISHal/blob/master/README.md) for instructions on installing the package.

Load the package.

	library(artiFISHal)

Create a data frame with 18 columns in which each row describes a sub-population of fish to be placed in the artificial lake. 
The first 11 columns must be completely filled in (no missing values). The last 8 columns may have some missing values. 
However, in each row, either water depth (WD and WDE) or distance to bottom (D2B and D2BE) must be filled in, but not both. 
See the help file for Simyfish for more information on the column names and descriptions.

	?Simyfish

You may find that the easiest way to create this data frame is by entering the information in an external file, 
e.g., a comma delimited file.  In that case, you would read in the data frame using the `read.csv` function.  
The data frame should be all character and numeric variables, so if you use `read.csv`, you should specify `as.is=TRUE`.

	ExInputs <- read.csv("C:/temp/ExInputs.csv", as.is=TRUE)

For the purposes of this vignette, I will create the data frame using code.

	ExInputs <- data.frame(
		Description = c("alewife small", "alewife small", "alewife large", "alewife large", 
			"alewife large", "alewife large", "rainbow smelt small", "rainbow smelt small", "rainbow smelt large", 
			"rainbow smelt large", "bloater small", "bloater small", "bloater large", "bloater large"), 
		G = c("a", "a", "A", "A", "A", "A", "s", "s", "S", "S", "b", "b", "B", "B"), 
		Z = c(65, 65, 127, 127, 127, 127, 55, 55, 127, 127, 70, 70, 190, 190), 
		ZE = c(0.15, 0.15, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.15, 0.15, 0.3, 0.3), 
		LWC1 = 1e-8*c(1400, 1400, 1400, 1400, 1400, 1400, 481, 481, 481, 481, 102, 102, 102, 102), 
		LWC2 = c(2.8638, 2.8638, 2.8638, 2.8638, 2.8638, 2.8638, 3.0331, 3.0331, 3.0331, 3.0331, 3.3812, 3.3812, 3.3812, 3.3812), 
		LWCE = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05), 
		TSC1 = c(-64.2, -64.2, -64.2, -64.2, -64.2, -64.2, -67.8, -67.8, -67.8, -67.8, -70.88, -70.88, -70.88, -70.88), 
		TSC2 = c(20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 19.9, 19.9, 19.9, 19.9, 25.54, 25.54, 25.54, 25.54), 
		TSCE = c(0.03, 0.03, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.08, 0.08, 0.03, 0.03, 0.04, 0.04), 
		PropN = c(0.3, 0.3, 0.035, 0.035, 0.015, 0.015, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375), 
		E = 100*c(93, 170, 6, 289, 90, 179, 93, 170, 6, 289, 93, 170, 34, 233), 
		EE = c(0.7, 0.4, 2.5, 0.1, 0.9, 0.6, 0.7, 0.4, 2.2, 0.1, 0.7, 0.4, 0.3, 0.1), 
		N = 100*c(100, 100, 100, 100, 100, 100, 100, 100, 170, 170, 100, 100, 75, 75), 
		NE = c(1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 1, 1, 0.5, 0.5), 
		WD = c(14, 14, 11, 11, 11, 11, 16, 16, 29, 29, 20, 20, NA, NA), 
		WDE = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.2, 0.2, NA, NA), 
		D2B = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 15, 15), 
		D2BE = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.3, 0.3))

Now you can use the `Simyfish` function to simulate the fish population based on these parameters.

	myfish <- Simyfish(LakeName="Lake M", LkWidth=30000, LkLength=20000, BotDepMin=20, BotDepMax=100, 
		FishParam=ExInputs, TotNFish=100000, Seed=545)

Look at the results.  The true total number and weight of each fish group in the population.

	myfish$Truth

The fish population that you just simulated.  Each row represents a single fish.  
Columns describe the fish group, location, and fish size.

	head(myfish$Fish)

Simulating a fish population that looks the way you intended may take some trial and error.
In the process of creating the population, `SimFish` will eliminate fish based on their size (LengthWeight, TS) or 
location (Easting, Northing, WaterDepth, BottomDepth, D2Shore).
You will likely end up with fewer fish than you requested.
But, if you end up with *far* fewer fish than expected, you should look at the proportion of fish excluded
for various reasons to troubleshoot where the problem might lie in the inputs provided.

	myfish$PropExcluded

The output from `SimFish` also keeps the inputs that were supplied as arguments for use in later sampling and summarization.

	myfish$LakeInfo
	myfish$FishInfo
	head(myfish$FishParam)

Once you have the population the way you want it, you can sample the population with a virtual acoustic and midwater trawl survey.

	mysurv <- SampFish(SimPop=myfish, NumEvents=10, AcNum=15, AcInterval=30000, AcLayer=10, AcAngle=7, 
		MtNum=30, MtHt=10, MtWd=10, MtLen=2000, Seed=341)

This will yeild two sample data frames, one for the acoustic targets detected and one for the midwater trawl catch.

	head(mysurv$Targets)
	head(mysurv$MtCatch)

The output from `SampFish` also keeps the inputs supplied as arguments for use in later sampling and summarization.

	mysurv$SurvParam

Estimate the number and biomass of fish per unit area sampled using the `AcMtEst` function
which combines the acoustic and midwater trawl data.

	perf <- AcMtEst(SimPop=myfish, AcMtSurv=Msurv, Seed=927)
	head(perf)

The virtual acoustic and midwater trawl surveys carried out by `SampFish` are assumed to be *perfect* with no issues of fish availability, 
acoustic dead zones, or trawl selectivity.
If you wish to introduce availability or selectivity functions, you may do so in the `AcMtEst` function.

For the target data, fish are excluded from the sample if they are located in zones unavailable to the acoustic transducer
specified with the argument `AcExcl`.  

For the catch data, fish are excluded from the sample if they are located in zones unavailable to the midwater trawl
specified with the argument `MtExcl`.  
In addition, selectivity curves can be specified for different panels of the midwater trawl using the arguments 
`PanelProps` and `SelecParam`.  
The `PanelProps` argument is used to denote the relative sizes of four different panels in the trawl,
from the mouth (outermost panel) to the cod end (innermost panel).
You can get a fish's eye view of the panel sizes denoted using the `ViewZones` function.

	ViewZones(PanelProp = c(0.4, 0.3, 0.2, 0.1))

The `SelecParam` argument is used to denote the mesh selectivity for each fish group (species, lifestage) and trawl panel zone.
A double logistic regression function is used to specify selectivity as a function of fish total length
using two slopes and two inflection points.
If you have some notion of what the selectivity curve should look like, 
you can use the `TuneSelec` function to discover the corresponding parameters using interactive *sliders*.

	TuneSelec()

Once you have chosen the selectivity parameters you want to use, put them in a data frame.

	selec <- data.frame(
		G = c("A", "a", "A", "a", "A", "a"),
		Zone = c("mouth", "mouth", "middle", "middle", "aft", "aft"),
		MtL50Small = c(100, 90, 60, 50, 30, 2),
		MtSlopeSmall = c(40, 40, 30, 30, 20, 20),
		MtL50Large = c(180, 180, Inf, Inf, Inf, Inf),
		MtSlopeLarge = c(20, 20, 100, 100, 100, 100))

Then you can visually compare the selectivity curves for the different fish groups using the `ViewSelec` function.

	ViewSelec(selec)

If all looks well, you can now derive estimates from the survey, incorporating these availabilities and selectivities.

	imperf <- AcMtEst(SimPop=myfish, AcMtSurv=Msurv, SelecParam=selec, AcExcl=c(5, 10), MtExcl=c(2, 2), Seed=204)
	head(imperf)
