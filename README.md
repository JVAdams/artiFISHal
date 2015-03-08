artiFISHal
==========

**artiFISHal** is a package of functions for the [R programming language](http://www.r-project.org/).  **artiFISHal** is a pelagic fish community simulator, which can be used to create artificial lakes and populate them with known numbers of fish (identified in species-size groups) to mimic pelagic fish communities. **artiFISHal** can then be used to sample these fish with virtual acoustic and midwater trawl surveys.

[Yule et al. (2013)](http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0072#.U1KYxPldXTQ) used **artiFISHal** to evaluate several different approaches for estimating the biomass of pelagic fish species in the Great Lakes. An example of how to use the functions in **artiFISHal** is given in this [vignette](https://github.com/JVAdams/artiFISHal/blob/master/Vignette.md).

- - -

You should be able to access the functions by installing them directly from within R.

	library("devtools")
	devtools::install_github("JVAdams/artiFISHal")
	library(artiFISHal)

If you don't already have `Rtools` and `devtools`, you will need to download and install (as administrator, if using a PC) `Rtools` from [CRAN](http://cran.r-project.org/bin/windows/Rtools/) then run the following lines of code before submitting the code above:

	find_rtools()
	install.packages("devtools")

An alternative approach for Windows users is to download this [zip file](https://github.com/JVAdams/artiFISHal/raw/master/a and install the package from the R menu:
- Packages
- Install package(s) from local zip files...
	
- - -

*NOTE:  If you are looking for files that were used in the Great Lakes Acoustic Users Group 2014 Workshop on Trawl Performance, they are now located at the [GLAUG2014](https://github.com/JVAdams/GLAUG2014) repository instead.

- - -

_Thanks to Hilary Parker whose blog post [Writing an R package from scratch](http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/) encouraged me to create my first R package._

- - -

_U.S. Geological Survey_ (USGS) Computer Program **artiFISHal** version 1.0.0. 
Written by Jean V. Adams, [USGS - Great Lakes Science Center](http://www.glsc.usgs.gov/), Ann Arbor, Michigan, USA. 
Written in programming language R (R Core Team, 2014, www.R-project.org), version 3.1.0 (2014-04-10). 
Run on a PC with Intel(R) Core(TM) I7-4600m CPU, 2.90 GHz processor, 16.0 GB RAM, and Microsoft Windows 7 Enterprise operating system 2009 Service Pack 1. 
Source code is available from Jean V. Adams on [GitHub](https://github.com/JVAdams/artiFISHal), _jvadams (at) usgs (dot) gov_.

_Disclaimer:_ Although this program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the United States Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.
