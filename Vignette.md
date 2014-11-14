artiFISHal Vignette
===================

The following is an example using the 
functions in the **artiFISHal** package.  See [README.md](https://github.com/JVAdams/artiFISHal/blob/master/README.md) 
for instructions on installing the package.

Start by loading the package.

	library(artiFISHal)

Create a data frame with 18 columns in which each row describes a sub-population of fish to be placed in the artificial lake. 
The first 11 columns must be completely filled in (no missing values). 

- **G** = character, a one-letter nickname for the group (e.g., fish species and lifestage) used in plotting
- **Z** = numeric, mean length (in mm)
- **ZE** = numeric, error around mean length, expressed as SD/mean
- **LWC1**, **LWC2** = numeric, length-weight regression coefficients, where wt = LWC1*len^LWC2, (wt in g, len in mm)
- **LWCE** = numeric, error around weight estimate, expressed as SD(estimate)/estimate
- **TSC1**, **TSC2** = numeric, target strength and length relation coefficients, ts = TSC1 + TSC2*log10(len/10), (ts in db, len in mm)
- **TSCE** = numeric, error around target strength estimate, expressed as SD(estimate)/estimate
- **PropN** = numeric, approximate proportion of population that the row represents (automatically adjusted to ensure they sum to 1)

The last 8 columns may have some missing values. 
However, in each row, either water depth (WD and WDE) or distance to bottom (D2B and D2BE) must be filled in, but not both. 

- **E** = numeric, mean easting (m)
- **EE** = numeric, error around easting, expressed as SD/mean
- **N** = numeric, mean northing (m)
- **NE** = numeric, error around northing, expressed as SD/mean
- **WD** = numeric, mean water depth (m)
- **WDE** = numeric, error around water depth, expressed as SD/mean
- **D2B** = numeric, mean distance to bottom (m)
- **D2BE** = numeric, error around distance to bottom, expressed as SD/mean

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

This data frame contains information on six groups of fish in 14 rows of data.
To get some idea of how this works, have a look at rows 3-6 which includes specifications for group A, alewife large.
The size of the fish being simulated are the same in all four rows (Z, ZE, LWC1, LWC3, LWCE, TSC1, TSC2, and TSCE).
But, the location of the fish are different: 3.5% are at easting 600 (with SD 2.5\*600), 1.5% are at easting 9000 (SD 0.9\*9000),
1.5% are at easting 17900 (with SD 0.6\*17900), and 3.5% are at easting 28900 (SD 0.1\*28900).
All of the groups of fish are located using water depth (WD and WDE), except for group B, bloater large,
which are located using distance to bottom (D2B and D2BE).
		
			   Description G   Z   ZE       LWC1   LWC2 LWCE   TSC1  TSC2 TSCE  PropN     E  EE     N  NE WD WDE D2B D2BE
	1        alewife small a  65 0.15 0.00001400 2.8638 0.05 -64.20 20.50 0.03 0.3000  9300 0.7 10000 1.0 14 0.2  NA   NA
	2        alewife small a  65 0.15 0.00001400 2.8638 0.05 -64.20 20.50 0.03 0.3000 17000 0.4 10000 1.0 14 0.2  NA   NA
	3        alewife large A 127 0.30 0.00001400 2.8638 0.05 -64.20 20.50 0.05 0.0350   600 2.5 10000 1.0 11 0.2  NA   NA
	4        alewife large A 127 0.30 0.00001400 2.8638 0.05 -64.20 20.50 0.05 0.0350 28900 0.1 10000 1.0 11 0.2  NA   NA
	5        alewife large A 127 0.30 0.00001400 2.8638 0.05 -64.20 20.50 0.05 0.0150  9000 0.9 10000 1.0 11 0.2  NA   NA
	6        alewife large A 127 0.30 0.00001400 2.8638 0.05 -64.20 20.50 0.05 0.0150 17900 0.6 10000 1.0 11 0.2  NA   NA
	7  rainbow smelt small s  55 0.30 0.00000481 3.0331 0.05 -67.80 19.90 0.05 0.0375  9300 0.7 10000 1.0 16 0.2  NA   NA
	8  rainbow smelt small s  55 0.30 0.00000481 3.0331 0.05 -67.80 19.90 0.05 0.0375 17000 0.4 10000 1.0 16 0.2  NA   NA
	9  rainbow smelt large S 127 0.30 0.00000481 3.0331 0.05 -67.80 19.90 0.08 0.0375   600 2.2 17000 0.1 29 0.3  NA   NA
	10 rainbow smelt large S 127 0.30 0.00000481 3.0331 0.05 -67.80 19.90 0.08 0.0375 28900 0.1 17000 0.1 29 0.3  NA   NA
	11       bloater small b  70 0.15 0.00000102 3.3812 0.05 -70.88 25.54 0.03 0.0375  9300 0.7 10000 1.0 20 0.2  NA   NA
	12       bloater small b  70 0.15 0.00000102 3.3812 0.05 -70.88 25.54 0.03 0.0375 17000 0.4 10000 1.0 20 0.2  NA   NA
	13       bloater large B 190 0.30 0.00000102 3.3812 0.05 -70.88 25.54 0.04 0.0375  3400 0.3  7500 0.5 NA  NA  15  0.3
	14       bloater large B 190 0.30 0.00000102 3.3812 0.05 -70.88 25.54 0.04 0.0375 23300 0.1  7500 0.5 NA  NA  15  0.3

Now you can use the `SimFish` function to simulate the fish population based on these parameters. 

	myfish <- SimFish(LakeName="My Lake", LkWidth=30000, LkLength=20000, BotDepMin=20, BotDepMax=100, 
		FishParam=ExInputs, TotNFish=100000, Seed=545)

Look at the results.  Look at the true total number and weight of each fish group in the population. 

	myfish$Truth
	
				  n         kg     nperha       kgperha
	a     38521  88.639197 0.64201667 0.00147731995
	A      4791  88.455773 0.07985000 0.00147426288
	b      4896   9.576029 0.08160000 0.00015960048
	B      7328 527.063943 0.12213333 0.00878439905
	s      4756   5.584219 0.07926667 0.00009307032
	S      3648  54.411041 0.06080000 0.00090685069
	Total 63940 773.730202 1.06566667 0.01289550337

Look at the fish population that you just simulated.  Each row represents a single fish. 
Columns describe the fish group, location, and fish size. 

	head(myfish$Fish)
	
	  G    f.east   f.north    f.d2sh  f.botdep   f.wdep  f.d2bot      len       wt        ts
	1 a 12355.212 16823.528 13466.323 100.00000 13.60617 86.39383 59.81102 1.595275 -49.49159
	2 a 12559.997 12749.377 13671.108 100.00000 12.48335 87.51665 64.93574 2.106270 -44.95919
	3 a  4237.675 14277.050  5348.786  96.27815 20.25294 76.02521 80.84223 4.227741 -46.51572
	4 a  2185.455  5228.477  3296.566  59.33820 16.64036 42.69784 52.59334 1.191047 -48.87585
	5 a  4965.597 17912.420  6076.708 100.00000 12.44079 87.55921 78.76110 3.766953 -46.85894
	9 a  9106.080 16052.540 10217.191 100.00000 12.18392 87.81608 64.67095 2.049004 -48.81689

Simulating a fish population that looks the way you intended may take some trial and error. 

![Dose-effect relation](https://github.com/JVAdams/LW1949/blob/master/Capture.PNG)

\href{https://raw.githubusercontent.com/JVAdams/artiFISHal/master/images/LakeFigures.JPG}{diagram}.

In the process of creating the population, `SimFish` will eliminate fish based on their size 
(Length or Weight less than zero or TS outside the allowable TS range) or 
location (Easting, Northing, WaterDepth, BottomDepth, D2Shore outside their allowable ranges). 
You will likely end up with fewer fish than you requested. 
But, if you end up with *far* fewer fish than expected, you should look at the proportion of fish excluded
to troubleshoot where in the inputs the problem might lie. 

	myfish$PropExcluded
	
	LengthWeight           TS      Easting     Northing   WaterDepth  BottomDepth      D2Shore 
		 0.00000      0.00020      0.09728      0.27628      0.09472      0.04808      0.08345 

The output from `SimFish` also keeps the inputs that were supplied as arguments for use in later sampling and summarization. 

	myfish$LakeInfo
	
	$LakeName
	[1] "My Lake"

	$LkWidth
	[1] 30000

	$LkLength
	[1] 20000

	$BotDepMin
	[1] 20

	$BotDepMax
	[1] 100

	$BotDepVertex
	[1] 200

	$ints
	[1]  20 290

	$slopes
	[1]  0.018 -0.009

	$d2shr.we
	[1] 1111.111 2222.222

	
	myfish$FishInfo
	
	$TotNFish
	[1] 100000

	$TSRange
	[1] -65 -20

	$Seed
	[1] 545


	head(myfish$FishParam)
	
	    Description G   Z   ZE     LWC1   LWC2 LWCE  TSC1 TSC2 TSCE PropN     E  EE     N NE WD WDE D2B D2BE Nfish
	1 alewife small a  65 0.15 0.000014 2.8638 0.05 -64.2 20.5 0.03 0.300  9300 0.7 10000  1 14 0.2  NA   NA 30000
	2 alewife small a  65 0.15 0.000014 2.8638 0.05 -64.2 20.5 0.03 0.300 17000 0.4 10000  1 14 0.2  NA   NA 30000
	3 alewife large A 127 0.30 0.000014 2.8638 0.05 -64.2 20.5 0.05 0.035   600 2.5 10000  1 11 0.2  NA   NA  3500
	4 alewife large A 127 0.30 0.000014 2.8638 0.05 -64.2 20.5 0.05 0.035 28900 0.1 10000  1 11 0.2  NA   NA  3500
	5 alewife large A 127 0.30 0.000014 2.8638 0.05 -64.2 20.5 0.05 0.015  9000 0.9 10000  1 11 0.2  NA   NA  1500
	6 alewife large A 127 0.30 0.000014 2.8638 0.05 -64.2 20.5 0.05 0.015 17900 0.6 10000  1 11 0.2  NA   NA  1500

Once you have the population the way you want it, you can sample the population with a virtual acoustic and midwater trawl survey. 

	mysurv <- SampFish(SimPop=myfish, NumEvents=10, AcNum=15, AcInterval=30000, AcLayer=10, AcAngle=7, 
		MtNum=30, MtHt=10, MtWd=10, MtLen=2000, Seed=341)

This will yield two sample data frames, one for the acoustic targets detected and one for the midwater trawl catch. 

	head(mysurv$Targets)
	head(mysurv$MtCatch)

The output from `SampFish` also keeps the inputs supplied as arguments for use in later sampling and summarization. 

	mysurv$SurvParam

Estimate the number and biomass of fish per unit area sampled using the `AcMtEst` function
which combines the acoustic and midwater trawl data. 

	perf <- AcMtEst(SimPop=myfish, AcMtSurv=mysurv, Seed=927)
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

The `SelecParam` argument is used to denote the mesh selectivity for each fish group (species, life stage) and trawl panel zone. 
You can use the `MeshPass` function to determine the largest height (or depth) of a fish that can pass through a 
single diamond-shaped mesh of a net, provided the size of the mesh `BarMesh`, 
the height to width ratio of the mesh `H2WRatio`,
and the length to height ratio of the fish `L2HRatio`.

	MeshPass(BarMesh=2, H2WRatio=0.3, L2HRatio=4)

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

	imperf <- AcMtEst(SimPop=myfish, AcMtSurv=mysurv, SelecParam=selec, AcExcl=c(5, 10), MtExcl=c(2, 2), Seed=204)
	head(imperf)
